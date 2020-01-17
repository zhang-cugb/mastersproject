import os
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
import scipy.sparse as sps
from porepy.models.contact_mechanics_biot_model import ContactMechanicsBiot

import GTS as gts


class ContactMechanicsBiotISC(ContactMechanicsBiot):
    """

    Attributes
    ----------
    name : str
        descriptive name of this class
    file_name : str
        Root name of solution files
    scalar_scale : float
        Scale of scalar variable
    length_scale : float
        Scale of lengths

    """
    def __init__(
            self,
            viz_folder_name: str,
            result_file_name: str,
            isc_data_path: str,
            mesh_args: Mapping[str, int],
            bounding_box: Mapping[str, int],
            shearzone_names: List[str],
            source_scalar_borehole_shearzone: Mapping[str, str],
            scales: Mapping[str, float],
    ):
        """ Initialize the Contact Mechanics Biot

        Parameters
        ----------
        viz_folder_name : str
            Absolute path to folder where grid and results will be stored
        result_file_name : str
            Root name for simulation result files
        isc_data_path : str
            Path to isc data: path/to/GTS/01BasicInputData
            Alternatively 'linux' or 'windows' for certain default paths (only applies to haakon's computers).
        mesh_args : Mapping[str, int]
            Arguments for meshing of domain.
            Required keys: 'mesh_size_frac', 'mesh_size_min, 'mesh_size_bound'
        bounding_box : Mapping[str, int]
            Bounding box of domain
            Required keys: 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'.
        shearzone_names : List[str]
            Which shear-zones to include in simulation
        source_scalar_borehole_shearzone : Mapping[str, str]
            Which borehole and shear-zone intersection to do injection in.
            Required keys: 'shearzone', 'borehole'
        scales : Mapping[str, float]
            Length scale and scalar variable scale.
            Required keys: 'scalar_scale', 'length_scale'
        """

        self.name = "contact mechanics biot on ISC dataset"
        logging.info(f"Running: {self.name}")

        if not os.path.exists(folder_name):
            os.makedirs(folder_name, exist_ok=True)
        logging.info(f"Visualization folder path: {folder_name}")

        params = {'folder_name': folder_name}

        super().__init__(params)

        # Root name of solution files
        self.file_name = result_file_name

        # Time
        self._set_time_parameters()

        # Scaling coefficients
        # self.scalar_scale = 10 * pp.GIGA
        # self.length_scale = 50

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
        self.mesh_args = mesh_args

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
        self.n_frac = 5

        # Basic parameters
        self.set_rock_and_fluid()

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
            # Use self.gb.node_props(g, 'name') to get value.
        else:
            assert self.Nd is not None

            # We require that 2D grids have a name.
            g = self.gb.get_grids(lambda g: g.dim == 2)
            for i, sz in enumerate(self.shearzone_names):
                assert self.gb.node_props(g[i], "name") is not None

    def faces_to_fix(self, g):
        """
        Identify three boundary faces to fix (u=0). This should allow us to assign
        Neumann "background stress" conditions on the rest of the boundary faces.
        Credits: PorePy paper
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

    def bc_type_mechanics(self, g):
        """
        We set Neumann values imitating an anisotropic background stress regime on all
        but the fracture faces, which are fixed to ensure a unique solution.
        credits: porepy article
        """
        all_bf, *_ = self.domain_boundary_sides(g)
        faces = self.faces_to_fix(g)
        bc = pp.BoundaryConditionVectorial(g, faces, "dir")
        frac_face = g.tags["fracture_faces"]
        bc.is_neu[:, frac_face] = False
        bc.is_dir[:, frac_face] = True
        return bc

    def bc_values_mechanics(self, g):
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
        # we, sn, bt = 9.2 * pp.MEGA * pp.PASCAL, 8.7 * pp.MEGA * pp.PASCAL, 13.1 * pp.MEGA * pp.PASCAL
        # we, sn, bt = 7 / 8, 5 / 4, 1
        # we = sn = bt = 9 * pp.MEGA * pp.PASCAL
        we = sn = bt = 9
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

    def bc_values_scalar(self, g):
        """ Hydrostatic flow values
        credit: porepy paper
        """
        # TODO: Hydrostatic scalar BC's (values).
        all_bf, *_ = self.domain_boundary_sides(g)
        bc_values = np.zeros(g.num_faces)
        depth = self._depth(g.face_centers[:, all_bf])
        bc_values[all_bf] = self.fluid.hydrostatic_pressure(depth) / self.scalar_scale
        # return bc_values
        return np.zeros(g.num_faces)

    def bc_type_scalar(self, g):
        """ Known boundary conditions (Dirichlet)
        """
        # TODO: Hydrostatic scalar BC's (type).
        # Define boundary regions
        all_bf, *_ = self.domain_boundary_sides(g)
        # Define boundary condition on faces
        return pp.BoundaryCondition(g, all_bf, "dir")

    def source_flow_rate(self):
        """
        Rate given in l / min = m^3/s * 1e-3 / 60.
        Length scaling needed to convert from the scaled length to m.
        """
        liters = 10
        return liters * pp.MILLI * (pp.METER / self.length_scale) ** self.Nd / pp.MINUTE

    def well_cells(self):
        """
        Tag well cells with unity values, positive for injection cells and
        negative for production cells.
        """
        df = self.isc.borehole_plane_intersection()
        # Borehole-shearzone intersection of interest
        bh_sz = self.source_scalar_borehole_shearzone

        _mask = (df.shearzone == bh_sz["shearzone"]) & (
                df.borehole == bh_sz["borehole"]
        )
        result = df.loc[_mask, ("x_sz", "y_sz", "z_sz")]
        if result.empty:
            raise ValueError("No intersection found.")

        pts = result.to_numpy().T
        assert pts.shape[1] == 1, "Should only be one intersection"

        for g, d in self.gb:
            tags = np.zeros(g.num_cells)

            # Get name of grid
            grid_name = self.gb.node_props(g, "name")

            # We only tag cells in the desired fracture
            if grid_name == bh_sz['shearzone']:
                logging.info(f"Grid of name: {grid_name}, and dimension {g.dim}")
                logging.info(f"Setting non-zero source for scalar variable")

                ids, dsts = g.closest_cell(pts, return_distance=True)
                logging.info(f"Closest cell found has distance: {dsts[0]:4f}")

                # Tag the injection cell
                tags[ids] = 1

            g.tags["well_cells"] = tags
            pp.set_state(d, {"well": tags.copy()})

    def source_scalar(self, g: pp.Grid):
        """ Well-bore source

        This is an example implementation of a borehole-fracture source.
        """
        flow_rate = self.source_flow_rate()

        # TODO: Ask if scalar source must be multiplied by time_step.
        values = flow_rate * g.tags["well_cells"] * self.time_step
        return values

    def source_mechanics(self, g):
        """
        Gravity term.
        Credits: PorePy paper
        """
        values = np.zeros((self.Nd, g.num_cells))
        values[2] = (
                pp.GRAVITY_ACCELERATION
                * self.rock.DENSITY
                * g.cell_volumes
                * self.length_scale
                / self.scalar_scale
        )
        return values.ravel("F")

    def compute_aperture(self, g):
        """ Compute aperture"""
        # TODO: Set aperture from transmissibilities.
        apertures = np.ones(g.num_cells)
        shearzone = self.gb.node_props(g, 'name')
        # The mean measured aperture per shear-zone.
        mean_apertures = {'S1_1': 1026 * pp.MILLI * pp.METER,
                          'S1_2': 909.5 * pp.MILLI * pp.METER,
                          'S1_3': 1634 * pp.MILLI * pp.METER,
                          'S3_1': 180 * pp.MILLI * pp.METER,
                          'S3_2': 159 * pp.MILLI * pp.METER,
                          None: 1,  # 3D matrix
                          }
        apertures *= mean_apertures[shearzone] / self.length_scale
        return apertures

    def set_permeability_from_aperture(self):
        """
        Cubic law in fractures, rock permeability in the matrix.
        Credits: PorePy paper
        """
        # TODO: Find a valid alternative to cubic law.
        viscosity = self.fluid.dynamic_viscosity() / self.scalar_scale
        gb = self.gb
        for g, d in gb:
            if g.dim < self.Nd:
                # Use cubic law in fractures
                apertures = self.compute_aperture(g)
                apertures_unscaled = apertures * self.length_scale
                k = np.power(apertures_unscaled, 2) / 12
                # Multiply with the cross-sectional area, which equals the apertures
                # for 2d fractures in 3d
                k *= apertures
                kxx = k / viscosity / self.length_scale ** 2
            else:
                # Use the rock permeability in the matrix
                kxx = (
                    self.rock.PERMEABILITY
                    / viscosity
                    * np.ones(g.num_cells)
                    / self.length_scale ** 2
                )

            K = pp.SecondOrderTensor(kxx)
            d[pp.PARAMETERS][self.scalar_parameter_key]["second_order_tensor"] = K

        # Normal permeability inherited from the neighboring fracture g_l
        for e, d in gb.edges():
            mg = d["mortar_grid"]
            g_l, _ = gb.nodes_of_edge(e)
            data_l = gb.node_props(g_l)
            a = self.compute_aperture(g_l)
            # We assume isotropic permeability in the fracture, i.e. the normal
            # permeability equals the tangential one
            k_s = data_l[pp.PARAMETERS][self.scalar_parameter_key][
                "second_order_tensor"
            ].values[0, 0]
            # Division through half the aperture represents taking the (normal) gradient
            kn = mg.slave_to_mortar_int() * np.divide(k_s, a / 2)
            pp.initialize_data(
                mg, d, self.scalar_parameter_key, {"normal_diffusivity": kn}
            )

    def set_rock_and_fluid(self):
        """
        Set rock and fluid properties to those of granite and water.
        We ignore all temperature effects.
        Credits: PorePy paper
        """
        self.rock = pp.Granite()

        # Lame parameters
        self.rock.YOUNG_MODULUS = 20.0 * pp.GIGA * pp.PASCAL
        self.rock.POISSON_RATIO = 0.33

        def lam_from(E, v):
            return E * v / ((1 + v) * (1 - 2 * v))

        def mu_from(E, v):
            return E / (2 * (1 + v))

        self.rock.LAMBDA = lam_from(self.rock.YOUNG_MODULUS, self.rock.POISSON_RATIO)
        self.rock.MU = mu_from(self.rock.YOUNG_MODULUS, self.rock.POISSON_RATIO)

        self.rock.FRICTION_COEFFICIENT = 0.2
        self.rock.POROSITY = 0.7 / 100

        self.fluid = pp.Water()
        # The permeability is for the intact rock
        self.rock.PERMEABILITY = 5e-19
        # Initial hydraulic aperture in m
        self.initial_aperture = 1e-3 / self.length_scale

    def set_mu(self, g):
        """ Set mu

        Set mu in linear elasticity stress-strain relation.
        stress = mu * trace(eps) + 2 * lam * eps
        """
        # TODO: Custom mu
        return np.ones(g.num_cells) * self.rock.MU

    def set_lam(self, g):
        """ Set lambda

        Set lambda in linear elasticity stress-strain relation.
        stress = mu * trace(eps) + 2 * lam * eps
        """
        # TODO: Custom lambda
        return np.ones(g.num_cells) * self.rock.LAMBDA

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

    def set_viz(self):
        """ Set exporter for visualization """
        self.viz = pp.Exporter(self.gb, name=self.file_name, folder=self.viz_folder_name)
        # list of time steps to export with visualization.
        self.export_times = []

        self.u_exp = 'u_exp'
        self.p_exp = 'p'
        self.traction_exp = 'traction_exp'
        self.export_fields = [
            self.u_exp,
            self.p_exp,
            # self.traction_exp,
        ]

    def _export_step(self):
        """
        export_step also serves as a hack to update parameters without changing the biot
        run method, since it is the only method of the setup class which is called at
        (the end of) each time step.
        """
        if "exporter" not in self.__dict__:
            self.set_viz()
        for g, d in self.gb:
            if g.dim == self.Nd:
                u = d[pp.STATE][self.displacement_variable].reshape(
                    (self.Nd, -1), order="F"
                )
                if g.dim == 3:
                    d[pp.STATE][self.u_exp] = u * self.length_scale

                else:
                    d[pp.STATE][self.u_exp] = np.vstack(
                        (u * self.length_scale, np.zeros(u.shape[1]))
                    )
                d[pp.STATE]["traction_exp"] = np.zeros(d[pp.STATE][self.u_exp].shape)
            else:
                g_h = self.gb.node_neighbors(g)[0]
                if g_h.dim == self.Nd:
                    data_edge = self.gb.edge_props((g, g_h))
                    u_mortar_local = self.reconstruct_local_displacement_jump(data_edge)
                    traction = d[pp.STATE][self.contact_traction_variable].reshape(
                        (self.Nd, -1), order="F"
                    )

                    if g.dim == 2:
                        d[pp.STATE][self.u_exp] = u_mortar_local * self.length_scale
                        d[pp.STATE][self.traction_exp] = traction
                    else:
                        d[pp.STATE][self.u_exp] = np.vstack(
                            (
                                u_mortar_local * self.length_scale,
                                np.zeros(u_mortar_local.shape[1]),
                            )
                        )
            d[pp.STATE][self.p_exp] = d[pp.STATE][self.scalar_variable] * self.scalar_scale
        self.viz.write_vtk(self.export_fields, time_step=self.time)
        self.export_times.append(self.time)
        self.save_data()
        # self.adjust_time_step()
        self.set_parameters()

    def save_data(self):
        n = self.n_frac
        if "u_jumps_tangential" not in self.__dict__:
            self.u_jumps_tangential = np.empty((1, n))
            self.u_jumps_normal = np.empty((1, n))
        tangential_u_jumps = np.zeros((1, n))
        normal_u_jumps = np.zeros((1, n))
        for g, d in self.gb:
            if g.dim < self.Nd:
                g_h = self.gb.node_neighbors(g)[0]
                # assert g_h.dim == self.Nd
                if not g_h.dim == self.Nd:
                    continue
                data_edge = self.gb.edge_props((g, g_h))
                u_mortar_local = self.reconstruct_local_displacement_jump(data_edge)
                tangential_jump = np.linalg.norm(
                    u_mortar_local[:-1] * self.length_scale, axis=0
                )
                normal_jump = u_mortar_local[-1] * self.length_scale
                vol = np.sum(g.cell_volumes)
                tangential_jump_norm = np.sqrt(np.sum(tangential_jump ** 2 * g.cell_volumes)) / vol
                normal_jump_norm = (
                        np.sqrt(np.sum(normal_jump ** 2 * g.cell_volumes)) / vol
                )
                tangential_u_jumps[0, g.frac_num] = tangential_jump_norm
                normal_u_jumps[0, g.frac_num] = normal_jump_norm

        self.u_jumps_tangential = np.concatenate((self.u_jumps_tangential, tangential_u_jumps))
        self.u_jumps_normal = np.concatenate((self.u_jumps_normal, normal_u_jumps))

    def export_step(self):
        """ Implementation of export step"""

        # Get fracture grids:
        # Set zero values there to facilitate export.
        frac_dims = [1, 2]
        for dim in frac_dims:
            gd_list = self.gb.grids_of_dimension(dim)
            for g in gd_list:
                data = self.gb.node_props(g)
                data[pp.STATE][self.u_exp] = np.zeros((3, g.num_cells))

        # Get the 3D data.
        g3 = self.gb.grids_of_dimension(3)[0]
        d3 = self.gb.node_props(g3)

        # Get the state, transform it, and save to another state variable
        u_exp = d3[pp.STATE][self.displacement_variable]
        d3[pp.STATE][self.u_exp] = np.reshape(np.copy(u_exp), newshape=(g3.dim, g3.num_cells), order="F")

        # Write step
        self.viz.write_vtk(self.export_fields, time_step=self.time)
        self.export_times.append(self.time)

        self.save_data()

    def export_pvd(self):
        """ Implementation of export pvd"""
        self.viz.write_pvd(self.export_times)

    def before_newton_loop(self):
        """ Will be run before entering a Newton loop.
        E.g.
           Discretize time-dependent quantities etc.
           Update time-dependent parameters (captured by assembly).
        """
        self.set_parameters()
        # The following is expensive, as it includes Biot. Consider making a custom  method
        # discretizing only the term you need!

        # TODO: Discretize only the terms you need.
        self.discretize()

    def after_newton_convergence(self, solution, errors, iteration_counter):
        """ Overwrite from parent to export solution steps."""
        self.assembler.distribute_variable(solution)
        self.save_mechanical_bc_values()
        self.export_step()

    def _set_time_parameters(self):
        """
        Set time parameters

        """
        self.time = 0
        self.time_step = 1 * pp.DAY
        self.end_time = 4 * pp.DAY
        # Set initial time step
        self.initial_time_step = self.time_step

        # num_steps = 2
        # self.time_step = 1 * self.length_scale ** 2
        # self.end_time = self.time_step * (num_steps - 1)
        # self.time_steps_array = np.linspace(start=0, stop=self.end_time, num=num_steps)
        # self.step_count = np.arange(len(self.time_steps_array))
        # self.current_step = self.step_count[0]

    def prepare_simulation(self):
        """ Is run prior to a time-stepping scheme. Use this to initialize
        discretizations, linear solvers etc.


        ONLY CHANGE FROM PARENT:
        - Set self.viz with custome method.
        """
        self.create_grid()
        self.Nd = self.gb.dim_max()
        self.well_cells()  # Tag the well cells
        self.set_parameters()
        self.assign_variables()
        self.assign_discretizations()
        self.initial_condition()
        self.discretize()
        self.initialize_linear_solver()

        self.set_viz()

    def _depth(self, coords):
        """
        Unscaled depth. We center the domain at 480m below the surface.
        (See Krietsch et al, 2018a)
        """
        return 480.0 * pp.METER - self.length_scale * coords[2]

    def _is_nonlinear_problem(self):
        """
        OVERWRITTEN FROM PARENT:
        This problem is non-linear, so return True.
        """
        return True


def main(
        viz_folder_name: str = None
):
    """Prepare the ContactMechanicsBiotISC solver for the porepy run_model method.

     Parameters
     viz_folder_name : str
        Absolute path to storage folder.
    """
    if viz_folder_name is None:
        viz_folder_name = (
            "/home/haakon/mastersproject/src/mastersproject/GTS/isc_modelling/results/cm_biot_1"
        )

    # Define mesh sizes for grid generation.
    mesh_size = 10  # .36
    # mesh_args = {
    #     "mesh_size_frac": mesh_size,
    #     "mesh_size_min": 0.1 * mesh_size,
    #     "mesh_size_bound": 6 * mesh_size,
    # }
    mesh_args = {
        "mesh_size_frac": mesh_size,
        "mesh_size_min": mesh_size,
        "mesh_size_bound": mesh_size,
    }

    setup = ContactMechanicsBiotISC(
        mesh_args=mesh_args,
        folder_name=viz_folder_name
    )

    # SOLVE THE PROBLEM
    default_options = {  # Parameters for Newton solver.
        "max_iterations": 2,
        "convergence_tol": 1,
        "divergence_tol": 1e5,
    }
    pp.run_time_dependent_model(setup=setup, params=default_options)

    setup.export_pvd()
    logging.info("Simulation done.")
    logging.info(f"All files stored to: \n {viz_folder_name}")

    return setup


def run_model(
        model: ContactMechanicsBiotISC = None,
        viz_folder_name: str = None,
        file_name: str = "test_biot"):
    """ Set up and run the biot model.

    Parameters
    model : ContactMechanicsBiotISC, Optional
        input model
    viz_folder_name : str
        Absolute path to storage folder.
    file_name : str
        root name of output files
    """
    if viz_folder_name is None:
        viz_folder_name = (
            "/home/haakon/mastersproject/src/mastersproject/GTS/isc_modelling/results/cm_biot_1"
        )

    # Define mesh sizes for grid generation.
    mesh_size = 5  # .36
    mesh_args = {
        "mesh_size_frac": mesh_size,
        "mesh_size_min": 0.1 * mesh_size,
        "mesh_size_bound": 6 * mesh_size,
    }

    if model is None:
        model = ContactMechanicsBiotISC(
            mesh_args=mesh_args,
            folder_name=viz_folder_name
        )

    model.prepare_simulation()
    model.set_viz()  # Overwrite the viz created in pp.contact_mechanics_biot at prepare_simulation()

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
    errors = []
    t_end = model.end_time
    k = 0
    while model.time < t_end:
        model.time += model.time_step
        k += 1
        logging.debug(
            f"\n Time step {k} at time {model.time:.1e} of {t_end:.1e} with time step {model.time_step:.1e}"
        )

        # Prepare for Newton

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
