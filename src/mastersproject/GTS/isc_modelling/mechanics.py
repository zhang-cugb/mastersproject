import porepy as pp
import numpy as np
from porepy.models.contact_mechanics_model import ContactMechanics
from porepy.models.abstract_model import AbstractModel


import GTS as gts
from GTS.prototype_1.mechanics.isotropic_setup import IsotropicSetup


class ContactMechanicsISC(ContactMechanics):
    """ Implementation of ContactMechanics for ISC"""

        # Time
        self._set_time_parameters()

        # Scaling coefficients
        self.scalar_scale = scales['scalar_scale']
        self.length_scale = scales['length_scale']

        # Grid
        self.gb = None
        self.Nd = None

        # --- BOUNDARY, INITIAL, SOURCE CONDITIONS ---
        self.source_scalar_borehole_shearzone = source_scalar_borehole_shearzone

        self.shearzone_names = shearzone_names
        self.n_frac = len(self.shearzone_names)

        self.mesh_args = mesh_args
        self.box = bounding_box

        self.isc_data_path = isc_data_path
        self.isc = gts.ISCData(path=self.isc_data_path)

    def create_grid(self, overwrite_grid=False):
        """ Create a GridBucket of a 3D domain with fractures
        defined by the ISC dataset.

        The method requires the following attribute:
            mesh_args (dict): Containing the mesh sizes.

        The method assigns the following attributes to self:
            gb (pp.GridBucket): The produced grid bucket.
            box (dict): The bounding box of the domain, defined through minimum and
                maximum values in each dimension.
            Nd (int): The dimension of the matrix, i.e., the highest dimension in the
                grid bucket.

        """
        if (self.gb is None) or overwrite_grid:
            network = gts.fracture_network(
                shearzone_names=self.shearzone_names,
                export=True,
                path=self.isc_data_path,
                domain=self.box,
            )
            path = f"{self.viz_folder_name}/gmsh_frac_file"
            self.gb = network.mesh(mesh_args=self.mesh_args, file_name=path)
            pp.contact_conditions.set_projections(self.gb)
            self.Nd = self.gb.dim_max()

            # TODO: Make this procedure "safe".
            #   E.g. assign names by comparing normal vector and centroid.
            #   Currently, we assume that fracture order is preserved in creation process.
            #   This may be untrue if fractures are (completely) split in the process.
            # Set fracture grid names:
            self.gb.add_node_props(keys="name")  # Add 'name' as node prop to all grids.
            fracture_grids = self.gb.get_grids(lambda g: g.dim == 2)
            for i, sz_name in enumerate(self.shearzone_names):
                self.gb.set_node_prop(fracture_grids[i], key="name", val=sz_name)
                # Note: Use self.gb.node_props(g, 'name') to get value.
        else:
            assert (self.Nd is not None), \
                "Attribute Nd must be set in an existing grid."

            # We require that 2D grids have a name.
            g = self.gb.get_grids(lambda g: g.dim == 2)
            for i, sz in enumerate(self.shearzone_names):
                assert (self.gb.node_props(g[i], "name") is not None), \
                    "All 2D grids must have a name."

    def faces_to_fix(self, g):
        """ Fix some boundary faces to dirichlet to ensure unique solution to problem.

        Identify three boundary faces to fix (u=0). This should allow us to assign
        Neumann "background stress" conditions on the rest of the boundary faces.

        Credits: Keilegavlen et al (2019) - Source code.
        """
        all_bf, *_ = self.domain_boundary_sides(g)
        point = np.array(
            [
                [(self.box["xmin"] + self.box["xmax"]) / 2],
                [(self.box["ymin"] + self.box["ymax"]) / 2],
                [self.box["zmin"]],
            ]
        )
        distances = pp.distances.point_pointset(point, g.face_centers[:, all_bf])
        indexes = np.argsort(distances)
        faces = all_bf[indexes[: self.Nd]]
        return faces

    def bc_type(self, g):
        """
        We set Neumann values on all but the fracture faces,
        which are fixed to ensure a unique solution.
        """
        all_bf, *_ = self.domain_boundary_sides(g)
        faces = self.faces_to_fix(g)
        bc = pp.BoundaryConditionVectorial(g, faces, ["dir"]*g.num_faces)
        fracture_faces = g.tags["fracture_faces"]
        bc.is_neu[:, fracture_faces] = False
        bc.is_dir[:, fracture_faces] = True
        return bc

    def stress_tensor(self):
        """ Compute the stress tensor"""

        # Values from Krietsch et al 2019
        stress_value = np.array([13.1, 9.2, 8.7])
        dip_direction = np.array([104.48, 259.05, 3.72])
        dip = np.array([39.21, 47.90, 12.89])

        def r(th, gm):
            """ Compute direction vector of a dip (th) and dip direction (gm)."""
            rad = np.pi / 180
            x = np.cos(th * rad) * np.sin(gm * rad)
            y = np.cos(th * rad) * np.cos(gm * rad)
            z = - np.sin(th * rad)
            return np.array([x, y, z])

        rot = r(th=dip, gm=dip_direction)

        # Orthogonalize the rotation matrix (which is already close to orthogonal)
        rot, _ = np.linalg.qr(rot)

        # Stress tensor in principal coordinate system
        stress = np.diag(stress_value)

        # Stress tensor in euclidean coordinate system
        stress_eucl = np.dot(np.dot(rot, stress), rot.T)
        return stress_eucl

    def bc_values(self, g):
        """ Stress values as ISC
        Credits: PorePy paper"""
        # Retrieve the boundaries where values are assigned
        all_bf, east, west, north, south, top, bottom = self.domain_boundary_sides(g)
        A = g.face_areas
        # Domain centred at 480 m below surface

        # Gravity acceleration
        gravity = (
                pp.GRAVITY_ACCELERATION
                * self.rock.DENSITY
                * self._depth(g.face_centers)
                / self.scalar_scale
        )
        bc_values = np.zeros((g.dim, g.num_faces))
        # TODO: Compute the actual (unperturbed) stress tensor
        we, sn, bt = 9.2 * pp.MEGA * pp.PASCAL, 8.7 * pp.MEGA * pp.PASCAL, 13.1 * pp.MEGA * pp.PASCAL

        # we = sn = bt = 9
        bc_values[0, west] = (we * gravity[west]) * A[west]
        bc_values[0, east] = -(we * gravity[east]) * A[east]
        bc_values[1, south] = (sn * gravity[south]) * A[south]
        bc_values[1, north] = -(sn * gravity[north]) * A[north]
        if self.Nd > 2:
            bc_values[2, bottom] = (bt * gravity[bottom]) * A[bottom]
            bc_values[2, top] = -(bt * gravity[top]) * A[top]

        faces = self.faces_to_fix(g)
        bc_values[:, faces] = 0

        return bc_values.ravel("F")


class ContactMechanicsIsotropicISC(IsotropicSetup):
    def __init__(self, **kwargs):
        """ Initialize the mechanics class for an ISC geometry."""
        super().__init__()

        # Mesh size arguments
        default_mesh_args = {"mesh_size_frac": 10, "mesh_size_min": 10}
        self.mesh_args = kwargs.get("mesh_args", default_mesh_args)

        # Bounding box of the domain
        default_box = {
            "xmin": -6,
            "xmax": 80,
            "ymin": 55,
            "ymax": 150,
            "zmin": 0,
            "zmax": 50,
        }
        self.box = kwargs.get("box", default_box)

    def create_grid(self):
        """ Create a GridBucket of a 3D domain with fractures
        defined by the ISC dataset.

        The method requires the following attribute:
            mesh_args (dict): Containing the mesh sizes.

        The method assigns the following attributes to self:
            gb (pp.GridBucket): The produced grid bucket.
            box (dict): The bounding box of the domain, defined through minimum and
                maximum values in each dimension.
            Nd (int): The dimension of the matrix, i.e., the highest dimension in the
                grid bucket.

        """
        network = gts.fracture_network(
            shearzone_names=None, export=True, path="linux", domain=self.box
        )
        self.gb = network.mesh(mesh_args=self.mesh_args)
        pp.contact_conditions.set_projections(self.gb)
        self.Nd = self.gb.dim_max()


def run_model(setup: AbstractModel):
    """
    Set up and run a model.

    Parameters
    setup : pp.AbstractModel
        A model that (ultimately) inherits from AbstractModel.

    """
    params = {"folder_name": "GTS/isc_modelling/results/isotropic_setup_viz"}
    model = setup(params=params)
    model.prepare_simulation()
    model.init_viz()

    tol = 1e-10

    # Get fracture grid(s):
    # Set zero values there to facilitate export.
    for i in [1, 2]:
        g_list = model.gb.grids_of_dimension(i)
        for g in g_list:
            data = model.gb.node_props(g)
            data[pp.STATE]["u_"] = np.zeros((3, g.num_cells))

    # Get the 3D data.
    g3 = model.gb.grids_of_dimension(3)[0]
    d3 = model.gb.node_props(g3)
    # Solve problem
    x = model.assemble_and_solve_linear_system(tol)
    model.after_newton_convergence(x, None, None)  # Distribute solution

    # Get the state, transform it, and save to another state variable
    sol3 = d3[pp.STATE][model.displacement_variable]
    trsol3 = np.reshape(np.copy(sol3), newshape=(g3.dim, g3.num_cells), order="F")
    d3[pp.STATE][model.displacement_variable + "_"] = trsol3

    # Export solution
    model.export_step()

    return model
