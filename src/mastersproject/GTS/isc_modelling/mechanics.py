import logging
from typing import (  # noqa
    Any,
    Coroutine,
    Generator,
    Generic,
    Iterable,
    List,
    Mapping,
    Optional,
    Set,
    Tuple,
    Type,
    TypeVar,
    Union,
)

import porepy as pp
import numpy as np
from porepy.models.contact_mechanics_model import ContactMechanics
from porepy.models.abstract_model import AbstractModel

import GTS as gts
from GTS.prototype_1.mechanics.isotropic_setup import IsotropicSetup

logger = logging.getLogger(__name__)


class ContactMechanicsISC(ContactMechanics):
    """ Implementation of ContactMechanics for ISC

    Run a Contact Mechanics model from porepy on the geometry
    defined by the In-Situ Stimulation and Circulation (ISC)
    project at the Grimsel Test Site (GTS).
    """

    def __init__(
            self,
            viz_folder_name: str,
            result_file_name: str,
            isc_data_path: str,
            mesh_args: Mapping[str, int],
            bounding_box: Mapping[str, int],
            shearzone_names: List[str],
            # source_scalar_borehole_shearzone: Mapping[str, str],
            scales: Mapping[str, float],
            stress: np.ndarray,
            solver: str,
    ):
        """ Initialize a Contact Mechanics model for GTS-ISC geometry.

        Parameters
        ----------
        viz_folder_name : str
            Absolute path to folder where grid and results will be stored
        result_file_name : str
            Root name for simulation result files
        isc_data_path : str
            Path to isc data: path/to/GTS/01BasicInputData
            Alternatively 'linux' or 'windows' for certain default paths (only applies to haakon's computers).
        --- SIMULATION RELATED PARAMETERS ---
        mesh_args : Mapping[str, int]
            Arguments for meshing of domain.
            Required keys: 'mesh_size_frac', 'mesh_size_min, 'mesh_size_bound'
        bounding_box : Mapping[str, int]
            Bounding box of domain
            Required keys: 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'.
        shearzone_names : List[str]
            Which shear-zones to include in simulation
        # source_scalar_borehole_shearzone : Mapping[str, str]
        #     Which borehole and shear-zone intersection to do injection in.
        #     Required keys: 'shearzone', 'borehole'
        scales : Mapping[str, float]
            Length scale and scalar variable scale.
            Required keys: 'scalar_scale', 'length_scale'
        solver : str, {'direct', 'pyamg'}
            Which solver to use
        --- PHYSICAL PARAMETERS ---
        stress : np.ndarray
            Stress tensor for boundary conditions
        """

        self.name = "contact mechanics on ISC dataset"
        logger.info(f"Running: {self.name}")

        params = {
            'folder_name': viz_folder_name,  # saved in self.viz_folder_name
            'linear_solver': solver,
        }
        super().__init__(params=params)

        # Root name of solution files
        self.file_name = result_file_name

        # Scaling coefficients
        self.scalar_scale = scales['scalar_scale']
        self.length_scale = scales['length_scale']

        # --- PHYSICAL PARAMETERS ---
        self.stress = stress
        self.set_rock()

        # --- BOUNDARY, INITIAL, SOURCE CONDITIONS ---
        # self.source_scalar_borehole_shearzone = source_scalar_borehole_shearzone

        # --- FRACTURES ---
        self.shearzone_names = shearzone_names
        self.n_frac = len(self.shearzone_names)
        # Initialize data storage for normal and tangential jumps
        self.u_jumps_tangential = np.empty((1, self.n_frac))
        self.u_jumps_normal = np.empty((1, self.n_frac))

        # --- COMPUTATIONAL MESH ---
        self.mesh_args = mesh_args
        self.box = bounding_box
        self.gb = None
        self.Nd = None

        # --- GTS-ISC DATA ---
        self.isc_data_path = isc_data_path
        self.isc = gts.ISCData(path=self.isc_data_path)

    def create_grid(self, overwrite_grid=False):
        """ Create a GridBucket of a 3D domain with fractures
        defined by the ISC dataset.

        Parameters
        ----------
        overwrite_grid : bool
            Overwrite an existing grid.

        The method requires the following attribute:
            mesh_args (dict): Containing the mesh sizes.

        Returns
        -------
        None

        Attributes
        ----------
        The method assigns the following attributes to self:
            gb : pp.GridBucket
                The produced grid bucket.
            box : dict
                The bounding box of the domain, defined through
                minimum and maximum values in each dimension.
            Nd : int
                The dimension of the matrix, i.e., the highest
                dimension in the grid bucket.

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
        We set Neumann values on all but a few boundary faces. Fracture faces also set to Dirichlet.

        Three boundary faces (see method faces_to_fix(self, g)) are set to 0 displacement (Dirichlet).
        This ensures a unique solution to the problem.
        Furthermore, the fracture faces are set to 0 displacement (Dirichlet).
        """

        all_bf, *_ = self.domain_boundary_sides(g)
        faces = self.faces_to_fix(g)
        bc = pp.BoundaryConditionVectorial(g, faces, ["dir"]*len(faces))
        fracture_faces = g.tags["fracture_faces"]
        bc.is_neu[:, fracture_faces] = False
        bc.is_dir[:, fracture_faces] = True
        return bc

    def bc_values(self, g):
        """ Mechanical stress values as ISC
        """
        # Retrieve the domain boundary
        all_bf, *_ = self.domain_boundary_sides(g)

        # Get outward facing normal vectors for domain boundary, weighted for face area
        # 1. Get normal vectors on the boundary
        bf_normals = g.face_normals[:, all_bf]
        # 2. Adjust direction so they face outwards
        flip_normal_to_outwards = np.where(g.cell_face_as_dense()[0, all_bf] >= 0, 1, -1)
        outward_normals = bf_normals * flip_normal_to_outwards

        # Boundary values
        bc_values = np.zeros((g.dim, g.num_faces))

        bf_stress = np.dot(self.stress, outward_normals)
        bc_values[:, all_bf] = bf_stress

        faces = self.faces_to_fix(g)
        bc_values[:, faces] = 0

        return bc_values.ravel("F")

    def source(self, g):
        """
        """
        return np.zeros(self.Nd * g.num_cells)

    def set_rock(self):
        """ Set rock properties of the ISC rock.
        """

        class GrimselGranodiorite(pp.UnitRock):
            def __init__(self):
                super().__init__()

                self.PERMEABILITY = 1
                self.THERMAL_EXPANSION = 1
                self.DENSITY = 2700 * pp.KILOGRAM / (pp.METER ** 3)
                self.POROSITY = 1

                # Lamé parameters
                self.YOUNG_MODULUS = 49 * pp.GIGA * pp.PASCAL  # Krietsch et al 2018 (Data Descriptor) - Dynamic E
                self.POISSON_RATIO = 0.32  # Krietsch et al 2018 (Data Descriptor) - Dynamic Poisson
                self.LAMBDA, self.MU = pp.params.rock.lame_from_young_poisson(
                    self.YOUNG_MODULUS, self.POISSON_RATIO
                )

                self.FRICTION_COEFFICIENT = 0.8
                self.POROSITY = 0.7 / 100

        self.rock = GrimselGranodiorite()

    def set_mu(self, g):
        """ Set mu

        Set mu in linear elasticity stress-strain relation.
        stress = mu * trace(eps) + 2 * lam * eps

        Assumes self.set_rock() is called
        """
        return np.ones(g.num_cells) * self.rock.MU

    def set_lam(self, g):
        """ Set lambda

        Set lambda in linear elasticity stress-strain relation.
        stress = mu * trace(eps) + 2 * lam * eps

        Assumes self.set_rock() is called
        """
        return np.ones(g.num_cells) * self.rock.LAMBDA

    def _set_friction_coefficient(self, g):
        """ The friction coefficient is uniform, and equal to 1.

        Assumes self.set_rock() is called
        """
        return np.ones(g.num_cells) * self.rock.FRICTION_COEFFICIENT

    def set_parameters(self):
        """
        Set the parameters for the simulation.
        """
        gb = self.gb

        for g, d in gb:
            if g.dim == self.Nd:
                # Rock parameters
                lam = self.set_lam(g)
                mu = self.set_mu(g)
                C = pp.FourthOrderTensor(mu, lam)

                # BC and source values
                bc = self.bc_type(g)
                bc_val = self.bc_values(g)
                source_val = self.source(g)

                pp.initialize_data(
                    g,
                    d,
                    self.mechanics_parameter_key,
                    {
                        "bc": bc,
                        "bc_values": bc_val,
                        "source": source_val,
                        "fourth_order_tensor": C,
                        # "max_memory": 7e7,
                        # "inverter": python,
                    },
                )

            elif g.dim == self.Nd - 1:
                friction = self._set_friction_coefficient(g)
                pp.initialize_data(
                    g,
                    d,
                    self.mechanics_parameter_key,
                    {"friction_coefficient": friction},
                )
        for _, d in gb.edges():
            mg = d["mortar_grid"]
            pp.initialize_data(mg, d, self.mechanics_parameter_key)

    def set_viz(self):
        """ Set exporter for visualization """
        self.viz = pp.Exporter(self.gb, name=self.file_name, folder=self.viz_folder_name)
        # list of time steps to export with visualization.

        self.u_exp = 'u_exp'
        self.traction_exp = 'traction_exp'
        self.export_fields = [
            self.u_exp,
            self.traction_exp,
        ]

    def prepare_simulation(self):
        """ Is run prior to a time-stepping scheme. Use this to initialize
        discretizations, linear solvers etc.
        """

        self.create_grid()
        self.Nd = self.gb.dim_max()
        self.set_parameters()
        self.assign_variables()
        self.assign_discretizations()
        self.initial_condition()
        self.discretize()
        self.initialize_linear_solver()

        self.set_viz()

    def export_step(self):
        """ Export a step

        Inspired by Keilegavlen 2019 (code)
        """
        gb = self.gb
        Nd = self.Nd

        for g, d in gb:
            if g.dim == Nd:  # On matrix
                u = d[pp.STATE][self.displacement_variable].reshape((Nd, -1), order='F').copy() * self.length_scale

                if g.dim != 3:  # Only called if solving a 2D problem
                    u = np.vstack(u, np.zeros(u.shape[1]))

                d[pp.STATE][self.u_exp] = u

                d[pp.STATE][self.traction_exp] = np.zeros(d[pp.STATE][self.u_exp].shape)

            else:  # In fractures or intersection of fractures (etc.)
                g_h = gb.node_neighbors(g, only_higher=True)[0]  # Get the higher-dimensional neighbor
                if g_h.dim == Nd:  # In a fracture
                    data_edge = gb.edge_props((g, g_h))
                    u_mortar_local = self.reconstruct_local_displacement_jump(
                        data_edge=data_edge, from_iterate=True).copy()
                    u_mortar_local = u_mortar_local * self.length_scale

                    traction = d[pp.STATE][self.contact_traction_variable].reshape((Nd, -1), order="F")

                    if g.dim == 2:
                        d[pp.STATE][self.u_exp] = u_mortar_local
                        d[pp.STATE][self.traction_exp] = traction
                    # TODO: Check when this statement is actually called
                    else:  # Only called if solving a 2D problem (i.e. this is a 0D fracture intersection)
                        d[pp.STATE][self.u_exp] = np.vstack(u_mortar_local, np.zeros(u_mortar_local.shape[1]))
                else:  # In a fracture intersection
                    d[pp.STATE][self.u_exp] = np.zeros((Nd, g.num_cells))
                    d[pp.STATE][self.traction_exp] = np.zeros((Nd, g.num_cells))
        self.viz.write_vtk(data=self.export_fields, time_dependent=False)  # Write visualization
        self.save_data()

    def save_data(self):
        """ Save normal and tangential jumps to a class attribute
        Inspired by Keilegavlen 2019 (code)
        """
        gb = self.gb
        Nd = self.Nd
        n = self.n_frac

        tangential_u_jumps = np.zeros((1, n))
        normal_u_jumps = np.zeros((1, n))

        for frac_num, frac_name in enumerate(self.shearzone_names):
            g_lst = gb.get_grids(lambda _g: gb.node_props(_g)['name'] == frac_name)
            assert len(g_lst) == 1  # Currently assume each fracture is uniquely named.

            g = g_lst[0]
            g_h = gb.node_neighbors(g, only_higher=True)[0]  # Get higher-dimensional neighbor
            assert g_h.dim == Nd  # We only operate on fractures of dim Nd-1.

            data_edge = gb.edge_props((g, g_h))
            u_mortar_local = self.reconstruct_local_displacement_jump(
                data_edge=data_edge, from_iterate=True).copy() * self.length_scale

            # Jump distances
            tangential_jump = np.linalg.norm(u_mortar_local[:-1, :], axis=0)
            normal_jump = np.linalg.norm(u_mortar_local[-1, :])

            # Ad-hoc average normal and tangential jump "estimates"
            # TODO: Find a proper way to express the "total" displacement of a fracture
            avg_tangential_jump = np.sum(tangential_jump * g.cell_volumes) / np.sum(g.cell_volumes)
            avg_normal_jump = np.sum(normal_jump * g.cell_volumes) / np.sum(g.cell_volumes)

            tangential_u_jumps[0, frac_num] = avg_tangential_jump
            normal_u_jumps[0, frac_num] = avg_normal_jump

        self.u_jumps_tangential = np.concatenate((self.u_jumps_tangential, tangential_u_jumps))
        self.u_jumps_normal = np.concatenate((self.u_jumps_normal, normal_u_jumps))

    def after_newton_iteration(self, solution_vector):
        """
        Extract parts of the solution for current iterate.

        The iterate solutions in d[pp.STATE]["previous_iterate"] are updated for the
        mortar displacements and contact traction are updated.
        Method is a tailored copy from assembler.distribute_variable.

        OVERWRITES parent to remove writing to vtk.

        Parameters:
            assembler (pp.Assembler): assembler for self.gb.
            solution_vector (np.array): solution vector for the current iterate.

        Returns:
            (np.array): displacement solution vector for the Nd grid.

        """
        self.update_state(solution_vector)

    def after_newton_convergence(self, solution, errors, iteration_counter):
        """ What to do at the end of a step."""
        self.assembler.distribute_variable(solution)
        self.export_step()


class ContactMechanicsISCWithGrid(ContactMechanicsISC):
    """ Solve contact mechanics with a pre-existing grid.
    """

    def __init__(
            self,
            viz_folder_name: str,
            result_file_name: str,
            isc_data_path: str,
            mesh_args: Mapping[str, int],
            bounding_box: Mapping[str, int],
            shearzone_names: List[str],
            # source_scalar_borehole_shearzone: Mapping[str, str],
            scales: Mapping[str, float],
            stress: np.ndarray,
            solver: str,
            gb: pp.GridBucket
    ):

        super().__init__(viz_folder_name, result_file_name, isc_data_path, mesh_args, bounding_box, shearzone_names,
                         scales, stress, solver)

        self.gb = gb
        self.Nd = gb.dim_max()

    def create_grid(self, overwrite_grid=False):
        # Overwrite method to ensure no new grid is created.
        assert self.gb is not None
        return


# --- OLD FILES --- DEPRECATED ---

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
