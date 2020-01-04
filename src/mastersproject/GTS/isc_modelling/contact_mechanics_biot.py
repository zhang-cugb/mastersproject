import os
import logging

import porepy as pp
import numpy as np
import scipy.sparse as sps
from porepy.models.contact_mechanics_biot_model import ContactMechanicsBiot

import GTS as gts


class ContactMechanicsBiotISC(ContactMechanicsBiot):
    def __init__(self, params=None, **kwargs):
        """ Initialize the Contact Mechanics Biot

        Parameters
        params : dict : Optional
            folder_name : str
                name of storage folder, relative to root.
            root : str
                root to folder name (path must exist)
                Default: /home/haakon/mastersproject/src/mastersproject/GTS/isc_modelling/results/
        kwargs
            time_step : float : Default = 1
                size of a time step (post-scaled to self.length_scale ** 2)
            num_steps : int : Default = 2
                Total number of time steps
            data_path : str
                path to isc_data: path/to/GTS/01BasicInputData
        """
        self.name = "contact mechanics biot on ISC dataset"
        logging.info(f"Running: {self.name}")

        # Specify absolute visualization storage path
        if params is None:
            params = {}
        assert isinstance(params, dict), "Params should be a dictionary."
        _root = (
            "/home/haakon/mastersproject/src/mastersproject/GTS/isc_modelling/results/"
        )
        root = params.get("root", _root)
        assert os.path.isdir(root)  # Root must be an existing path
        self.path = root + params.get("folder_name", "biot_contact_mechanics_viz")
        params["folder_name"] = self.path
        logging.info(f"Visualization folder path: {self.path}")

        super().__init__(params)
        self.viz = None

        # Time
        num_steps = kwargs.get("num_steps", 2)
        self.time_step = kwargs.get("time_step", 1) * self.length_scale ** 2
        self.end_time = self.time_step * (num_steps - 1)
        self.time_steps_array = np.linspace(start=0, stop=self.end_time, num=num_steps)
        self.step_count = np.arange(len(self.time_steps_array))
        self.current_step = self.step_count[0]

        # Grid
        self.gb = None
        self.Nd = None

        # Boundary conditions, initial conditions, source conditions:
        # Scalar source
        self.source_scalar_borehole_shearzone = {
            "shearzone": "S1_1",
            "borehole": "INJ1",
        }

        # Fractures are created in the order of self.shearzone_names.
        # This is effectively an index of the shearzone at hand.
        default_shearzone_set = ["S1_1", "S1_2", "S1_3", "S3_1", "S3_2"]
        self.shearzone_names = kwargs.get("shearzone_names", default_shearzone_set)

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

        # TODO: Think of a good way to include ISCData in this class
        self.isc = gts.ISCData(path=kwargs.get("data_path", "linux"))

    def create_grid(self, overwrite_grid=False):
        """ Create a GridBucket of a 3D domain with fractures
        defined by the ISC dataset.

        Parameters
        overwrite_grid : bool
            Overwrite an existing grid.

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
                path="linux",
                domain=self.box,
            )
            path = f"{self.path}/gmsh_frac_file"
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
            # Use self.gb.node_props(g, 'name') to get value.
        else:
            assert self.Nd is not None

            # We require that 2D grids have a name.
            g = self.gb.get_grids(lambda g: g.dim == 2)
            for i, sz in enumerate(self.shearzone_names):
                assert self.gb.node_props(g[i], "name") is not None

    def bc_type_mechanics(self, g):
        # TODO: Custom mechanics boundary conditions (type).
        return super().bc_type_mechanics(g)

    def bc_values_mechanics(self, g):
        # TODO: Customer mechanics boundary conditions (values).
        return super().bc_values_mechanics(g)

    def bc_values_scalar(self, g):
        """ Set boundary values to 1 (Neumann) on top face.
        0 (Dirichlet) on bottom face.
        0 (Neumann) otherwise.
        """
        # TODO: Hydrostatic scalar BC's (values).
        all_bf, east, west, north, south, top, bottom = self.domain_boundary_sides(g)
        top_face = np.nonzero(top)[0]
        bc_val = np.zeros(g.num_faces)
        bc_val[top_face] = 1
        return bc_val

    def bc_type_scalar(self, g):
        """ Set boundary conditions dirichlet on bottom face.
        Neumann otherwise.
        """
        # TODO: Hydrostatic scalar BC's (type).
        # Define boundary regions
        all_bf, east, west, north, south, top, bottom = self.domain_boundary_sides(g)
        bottom_face = np.nonzero(bottom)[0]
        # Define boundary condition on faces
        return pp.BoundaryCondition(g, bottom_face, "dir")

    def source_scalar(self, g: pp.Grid):
        """ Well-bore source

        This is an example implementation of a borehole-fracture source.
        """
        # Borehole-shearzone intersection of interest
        bh_sz = self.source_scalar_borehole_shearzone

        # Get name of grid
        grid_name = self.gb.node_props(g, "name")

        # 0-values for all non-2D grids
        if grid_name is None:
            return np.zeros(g.num_cells)

        # Also 0-values for all 2D-grids except the 'S1_1' shear-zone.
        if grid_name != bh_sz["shearzone"]:
            return np.zeros(g.num_cells)

        logging.info(f"Grid of name: {grid_name}, and dimension {g.dim}")
        logging.info(f"Setting non-zero source for scalar variable")

        # Get necessary data
        df = self.isc.borehole_plane_intersection()

        _mask = (df.shearzone == bh_sz["shearzone"]) & (
            df.borehole == bh_sz["borehole"]
        )
        result = df.loc[_mask, ("x_sz", "y_sz", "z_sz")]
        if result.empty:
            raise ValueError("No intersection found.")

        pts = result.to_numpy().T
        assert pts.shape[1] == 1, "Should only be one intersection"

        # Find cell nearest the desired point.
        ids, dsts = g.closest_cell(pts, return_distance=True)
        logging.info(f"Closest cell found has distance: {dsts[0]:4f}")

        # Set the source term.
        values = np.zeros(g.num_cells)
        values[ids] = 10
        return values

    def set_mu(self, g):
        """ Set mu

        Set mu in linear elasticity stress-strain relation.
        stress = mu * trace(eps) + 2 * lam * eps
        """
        # TODO: Custom mu
        return np.ones(g.num_cells)

    def set_lam(self, g):
        """ Set lambda

        Set lambda in linear elasticity stress-strain relation.
        stress = mu * trace(eps) + 2 * lam * eps
        """
        # TODO: Custom lambda
        return np.ones(g.num_cells)

    def set_mechanics_parameters(self):
        """
        Set the parameters for the simulation.
        """
        gb = self.gb

        for g, d in gb:
            if g.dim == self.Nd:
                # Rock parameters
                lam = self.set_lam(g) / self.scalar_scale
                mu = self.set_mu(g) / self.scalar_scale
                C = pp.FourthOrderTensor(mu, lam)

                # Define boundary condition
                bc = self.bc_type_mechanics(g)
                # BC and source values
                bc_val = self.bc_values_mechanics(g)
                source_val = self.source_mechanics(g)

                pp.initialize_data(
                    g,
                    d,
                    self.mechanics_parameter_key,
                    {
                        "bc": bc,
                        "bc_values": bc_val,
                        "source": source_val,
                        "fourth_order_tensor": C,
                        "time_step": self.time_step,
                        "biot_alpha": self.biot_alpha(g),
                    },
                )

            elif g.dim == self.Nd - 1:
                friction = self._set_friction_coefficient(g)
                pp.initialize_data(
                    g,
                    d,
                    self.mechanics_parameter_key,
                    {"friction_coefficient": friction, "time_step": self.time_step},
                )

        for _, d in gb.edges():
            mg = d["mortar_grid"]
            pp.initialize_data(mg, d, self.mechanics_parameter_key)

    def init_viz(self, file_name="test_biot", overwrite=False):
        """ Initialize visualization.
        Will only create a new object if none exists.
        Alternatively, the existing exporter can be overwritten using 'overwrite'.

        """
        if (self.viz is None) or overwrite:
            # g3 = self.gb.grids_of_dimension(self.gb.dim_max())[0]
            self.viz = pp.Exporter(self.gb, name=file_name, folder=self.viz_folder_name)

    def export_step(self):
        """ Implementation of export step"""
        export_fields = [self.displacement_variable + "_"]  # self.scalar_variable
        # Test out: Export a single grid.
        # g3 = self.gb.grids_of_dimension(self.gb.dim_max())[0]
        # data = self.gb.node_props(g3)
        # export_data = data[pp.STATE][self.displacement_variable]
        # export_fields = {self.displacement_variable: export_data}
        self.viz.write_vtk(export_fields, time_step=self.current_step)

    def export_pvd(self):
        """ Implementation of export pvd"""
        num_steps = np.arange(len(self.time_steps_array))
        self.viz.write_pvd(num_steps)


def run_model(model: ContactMechanicsBiotISC = None, viz_folder_name="biot"):
    """ Set up and run the biot model.

    Parameters
    model : ContactMechanicsBiotISC, Optional
        input model
    viz_folder_name : str
        Path to storage folder.
    """
    if model is None:
        params = {"folder_name": viz_folder_name}
        model = ContactMechanicsBiotISC(params=params)
    model.prepare_simulation()
    model.init_viz(
        overwrite=True
    )  # Overwrite the viz created in pp.contact_mechanics_biot at prepare_simulation()
    time_steps = model.time_steps_array

    # breakpoint()
    print("Starting simulation...")
    tol = 1e-10

    # Get fracture grid(s):
    # Set zero values there to facilitate export.
    frac_dims = [1, 2]
    for dim in frac_dims:
        gd_list = model.gb.grids_of_dimension(dim)
        for g in gd_list:
            data = model.gb.node_props(g)
            data[pp.STATE]["u_"] = np.zeros((3, g.num_cells))

    # Get the 3D data.
    g3 = model.gb.grids_of_dimension(3)[0]
    d3 = model.gb.node_props(g3)
    # Solve problem
    for curr_step, step in enumerate(time_steps):
        model.time = step
        model.current_step = curr_step
        x = model.assemble_and_solve_linear_system(tol)  # Solve time step
        # TODO: Overwrite method and save errors and iteration counter.
        model.after_newton_convergence(x, None, None)  # Distribute solution

        # Get the state, transform it, and save to another state variable
        sol3 = d3[pp.STATE][model.displacement_variable]
        trsol3 = np.reshape(np.copy(sol3), newshape=(g3.dim, g3.num_cells), order="F")
        d3[pp.STATE][model.displacement_variable + "_"] = trsol3

        # breakpoint()

        model.export_step()

    print("Successful simulation.")
    model.export_pvd()

    return model
