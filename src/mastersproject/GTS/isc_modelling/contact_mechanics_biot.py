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
from GTS.isc_modelling.mechanics import ContactMechanicsISC

import GTS as gts

# --- LOGGING UTIL ---
from util.logging_util import timer, trace

logger = logging.getLogger(__name__)


# TODO: Re-introduce all scaling when it is properly understood.
class ContactMechanicsBiotISC(ContactMechanicsISC, ContactMechanicsBiot):
    """
    TODO: Write class description
    """

    def __init__(self, params: dict):
        """ Initialize the Contact Mechanics Biot

        Parameters
        ----------
        params : dict
            Should contain the following key-value pairs:
                viz_folder_name : str
                    Absolute path to folder where grid and results will be stored
                mesh_args : dict[str, int]
                    Arguments for meshing of domain.
                    Required keys: 'mesh_size_frac', 'mesh_size_min, 'mesh_size_bound'
                bounding_box : d[str, int]
                    Bounding box of domain
                    Required keys: 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'.
                shearzone_names : List[str]
                    Which shear-zones to include in simulation
                source_scalar_borehole_shearzone : dict[str, str]
                    Which borehole and shear-zone intersection to do injection in.
                    Required keys: 'shearzone', 'borehole'
                length_scale, scalar_scale : float : Optional (Default: 100, pp.GIGA, respectively)
                    Length scale and scalar variable scale.
        """

        logger.info(f"Initializing contact mechanics biot on ISC dataset")

        # --- BOUNDARY, INITIAL, SOURCE CONDITIONS ---
        self.source_scalar_borehole_shearzone = params.get('source_scalar_borehole_shearzone')

        super().__init__(params=params)

        # Set file name of the pre-run first.
        self.file_name = 'initialize_run'

        # Time
        self.prepare_initial_run()

        # Initialize phase and injection rate
        self.current_phase = 0
        self.current_injection_rate = 0

        # --- PHYSICAL PARAMETERS ---
        self.set_rock_and_fluid()

        self.transmissivity = {  # Unscaled
            'S1_1': 1e-12 * pp.METER ** 2 / pp.SECOND,
            'S1_2': 1e-12 * pp.METER ** 2 / pp.SECOND,
            'S1_3': 1e-12 * pp.METER ** 2 / pp.SECOND,
            'S3_1': 1e-6 * pp.METER ** 2 / pp.SECOND,
            'S3_2': 1e-6 * pp.METER ** 2 / pp.SECOND,
            None: 1e-14 * pp.METER ** 2 / pp.SECOND,  # 3D matrix
        }
        self.aquifer_thickness = {  # Unscaled
            'S1_1': 1000 * pp.MILLI * pp.METER,
            'S1_2': 1000 * pp.MILLI * pp.METER,
            'S1_3': 1000 * pp.MILLI * pp.METER,
            'S3_1': 170 * pp.MILLI * pp.METER,
            'S3_2': 170 * pp.MILLI * pp.METER,
            None: 1,  # 3D matrix
        }

    def bc_type_mechanics(self, g):
        """
        We set Neumann values on all but a few boundary faces. Fracture faces also set to Dirichlet.

        Three boundary faces (see method faces_to_fix(self, g)) are set to 0 displacement (Dirichlet).
        This ensures a unique solution to the problem.
        Furthermore, the fracture faces are set to 0 displacement (Dirichlet).
        """
        return super().bc_type(g)

    def bc_values_mechanics(self, g):
        """ Scaled mechanical stress values as ISC
        """
        return super().bc_values(g)

    def bc_values_scalar(self, g):
        """ Hydrostatic flow values
        credit: porepy paper
        """
        # TODO: Hydrostatic scalar BC's (values).
        all_bf, *_ = self.domain_boundary_sides(g)
        bc_values = np.zeros(g.num_faces)
        depth = self._depth(g.face_centers[:, all_bf])

        # DIRICHLET
        bc_values[all_bf] = self.fluid.hydrostatic_pressure(depth) / self.scalar_scale
        return bc_values
        # return np.zeros(g.num_faces)

    def bc_type_scalar(self, g):
        """ Known boundary conditions (Dirichlet)
        """
        # Define boundary regions
        all_bf, *_ = self.domain_boundary_sides(g)
        # Define boundary condition on faces
        return pp.BoundaryCondition(g, all_bf, ["dir"] * all_bf.size)

    def source_flow_rate(self):
        """ Scaled source flow rate

        [From Grimsel Experiment Description]:

        We simulate one part of the injection procedure (Phase 3).
        In this phase, they inject by flow rate.
        Four injection steps of 10, 15, 20 and 25 l/min, each for 10 minutes.

        Afterwards, the system is shut-in for 40 minutes.
        """
        self.simulation_protocol()
        injection_rate = self.current_injection_rate
        return injection_rate * pp.MILLI * (pp.METER / self.length_scale) ** self.Nd

    def well_cells(self):
        """
        Tag well cells with unity values, positive for injection cells and
        negative for production cells.
        """
        # TODO: Use unscaled grid to find result.
        df = self.isc.borehole_plane_intersection()
        # Borehole-shearzone intersection of interest
        bh_sz = self.source_scalar_borehole_shearzone

        _mask = (df.shearzone == bh_sz["shearzone"]) & (
                df.borehole == bh_sz["borehole"]
        )
        result = df.loc[_mask, ("x_sz", "y_sz", "z_sz")]
        if result.empty:
            raise ValueError("No intersection found.")

        pts = result.to_numpy().T / self.length_scale
        assert pts.shape[1] == 1, "Should only be one intersection"
        tagged = False

        for g, d in self.gb:
            tags = np.zeros(g.num_cells)

            # Get name of grid
            grid_name = self.gb.node_props(g, "name")

            # We only tag cells in the desired fracture
            if grid_name == bh_sz['shearzone']:
                logger.info(f"Tagging grid of name: {grid_name}, and dimension {g.dim}")
                logger.info(f"Setting non-zero source value for pressure")

                ids, dsts = g.closest_cell(pts, return_distance=True)
                logger.info(f"Closest cell found has distance: {dsts[0]:4f}")

                # Tag the injection cell
                tags[ids] = 1
                tagged = True

            g.tags["well_cells"] = tags
            pp.set_state(d, {"well": tags.copy()})

        if not tagged:
            logger.warning("No injection cell was tagged.")

    def source_scalar(self, g: pp.Grid):
        """ Well-bore source

        This is an example implementation of a borehole-fracture source.
        """
        flow_rate = self.source_flow_rate()  # Already scaled by self.length_scale
        values = flow_rate * g.tags["well_cells"] * self.time_step
        # TODO: Hydrostatic pressure
        return values

    def source_mechanics(self, g):
        """ Scaled gravity term. """
        return super().source(g)

    def _permeability_from_transmissivity(self, T, b, theta=None):
        """ Compute permeability [m2] from transmissivity [m2/s]

        We can relate permeability, k [m2] with transmissivity, T [m2/s]
        through the relation
        k = T * mu / (rho * g * b)
        where mu is dynamic viscosity [Pa s], rho is density [kg/m3],
        g in gravitational acceleration [m/s2] and b is aquifer thickness [m]

        Assumes that self.fluid is set
        """
        mu = self.fluid.dynamic_viscosity(theta=theta)
        rho = self.fluid.density(theta=theta)
        g = pp.GRAVITY_ACCELERATION
        k = T * mu / (rho * g * b)
        return k

    def grid_permeability_from_transmissivity(self, g, theta=None):
        """ Wrapper for permeability_from_transmissivity
        Returns the uniform permeability over the grid g (np.array of size g.num_cells)
        """
        shearzone = self.gb.node_props(g, 'name')
        T = self.transmissivity[shearzone]  # Unscaled
        b = self.aquifer_thickness[shearzone]  # Unscaled
        permeability = self._permeability_from_transmissivity(T, b, theta=theta)  # Unscaled
        return permeability * np.ones(g.num_cells)

    def set_permeability_from_transmissivity(self):
        """ Set permeability in fracture and matrix from transmissivity

        The formula,
        k = T * mu / (rho * g * b)
        where b is aquifer thickness, and T is transmissivity (see above)
        is used.

        In shear zones, we use its thickness as b.
        In the 3D rock matrix, we use 1.
        """

        viscosity = self.fluid.dynamic_viscosity() / self.scalar_scale
        gb = self.gb
        for g, d in gb:
            k = self.grid_permeability_from_transmissivity(g) / (self.length_scale ** 2)

            kxx = k / viscosity

            k_tensor = pp.SecondOrderTensor(kxx)
            d[pp.PARAMETERS][self.scalar_parameter_key]["second_order_tensor"] = k_tensor

        # TODO: Understand how permeability works on the mortar grid.
        # Normal permeability inherited from the neighboring fracture g_l
        for e, d in gb.edges():
            mg = d["mortar_grid"]
            g_l, _ = gb.nodes_of_edge(e)
            data_l = gb.node_props(g_l)

            a = self.grid_aperture_from_transmissivity(g_l)  # Unscaled

            # We assume isotropic permeability in the fracture, i.e. the normal
            # permeability equals the tangential one
            k_s = data_l[pp.PARAMETERS][self.scalar_parameter_key][
                "second_order_tensor"
            ].values[0, 0]
            # Division through half the aperture represents taking the (normal) gradient
            # TODO: Check scaling for 'kn'
            kn = mg.slave_to_mortar_int() * np.divide(k_s, a / 2) * self.scalar_scale
            d = pp.initialize_data(
                e, d, self.scalar_parameter_key, {"normal_diffusivity": kn}
            )  # TODO: Check if it should be mg

    def _aperture_from_transmissivity(self, T, b, theta=None):
        """ Compute hydraulic aperture [m] from transmissivity [m2/s]

        We use the following relation (derived from cubic law):
        a = sqrt( 12 * mu * T / (rho * g * b) )
        where mu is dynamic viscosity [Pa s], rho is density [kg/m3],
        g is gravitational acceleration [m/s2], and b is aquifer thickness [m]

        Assumes that self.fluid is set
        """
        mu = self.fluid.dynamic_viscosity(theta=theta)
        rho = self.fluid.density(theta=theta)
        g = pp.GRAVITY_ACCELERATION
        hydraulic_aperture = np.sqrt(12 * mu * T / (rho * g * b))
        return hydraulic_aperture

    def grid_aperture_from_transmissivity(self, g: pp.Grid, theta=None):
        """ Grid wrapper for aperture_from_transmissivity
        Returns the uniform aperture over the grid g (np.array of size g.num_cells)
        """
        shearzone = self.gb.node_props(g, 'name')
        T = self.transmissivity[shearzone]
        b = self.aquifer_thickness[shearzone]
        aperture = self._aperture_from_transmissivity(T, b, theta=theta)  # Unscaled
        return aperture * g.num_cells

    def set_rock_and_fluid(self):
        """
        Set rock and fluid properties to those of granite and water.
        We ignore all temperature effects.
        Credits: PorePy paper
        """

        super().set_rock()

        # Fluid. Temperature at ISC is 11 degrees average.
        self.fluid = pp.Water(theta_ref=11)

    def set_parameters(self):
        """ Set biot parameters
        """
        self.set_mechanics_parameters()
        self.set_scalar_parameters()

    def set_mechanics_parameters(self):
        """ Set mechanics parameters for the simulation.
        """
        # TODO Consider calling super().set_parameters(),
        #  then set the remaining parameters here.
        gb = self.gb

        for g, d in gb:
            if g.dim == self.Nd:
                # Rock parameters
                lam = self.rock.LAMBDA * np.ones(g.num_cells) / self.scalar_scale
                mu = self.rock.MU * np.ones(g.num_cells) / self.scalar_scale
                C = pp.FourthOrderTensor(mu, lam)

                # Define boundary condition
                bc = self.bc_type_mechanics(g)
                # BC and source values
                bc_val = self.bc_values_mechanics(g)  # Already scaled
                source_val = self.source_mechanics(g)  # Already scaled

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

    def set_scalar_parameters(self):
        """ Set scalar parameters for the simulation
        """
        gb = self.gb

        compressibility = self.fluid.COMPRESSIBILITY * self.scalar_scale  # Units [1/Pa]
        porosity = self.rock.POROSITY
        for g, d in gb:
            # specific volume
            aperture = self.grid_aperture_from_transmissivity(g)
            specific_volume = np.power(aperture, self.Nd - g.dim)

            # Boundary and source conditions
            bc = self.bc_type_scalar(g)
            bc_values = self.bc_values_scalar(g)  # Already scaled
            source_values = self.source_scalar(g)  # Already scaled

            # Biot alpha
            alpha = self.biot_alpha(g)

            # Initialize data
            pp.initialize_data(
                g,
                d,
                self.scalar_parameter_key,
                {
                    "bc": bc,
                    "bc_values": bc_values,
                    "mass_weight": compressibility * porosity * specific_volume,
                    "biot_alpha": alpha,
                    "source": source_values,
                    "time_step": self.time_step,
                },
            )

        # Set permeability on grid, fracture and mortar grids.
        self.set_permeability_from_transmissivity()

    def set_viz(self):
        """ Set exporter for visualization """
        self.viz = pp.Exporter(self.gb, file_name=self.file_name, folder_name=self.viz_folder_name)
        # list of time steps to export with visualization.
        self.export_times = []

        self.u_exp = 'u_exp'
        self.p_exp = 'p_exp'
        self.traction_exp = 'traction_exp'
        self.normal_frac_u = 'normal_frac_u'
        self.tangential_frac_u = 'tangential_frac_u'

        self.export_fields = [
            self.u_exp,
            self.p_exp,
            # self.traction_exp,
            self.normal_frac_u,
            self.tangential_frac_u,
        ]

    @trace(logger, timeit=False)
    def export_step(self):
        """ Export a step

        Inspired by Keilegavlen 2019 (code)
        """
        # TODO: Check that everything is unscaled correctly

        self.save_frac_jump_data()  # Save fracture jump data to pp.STATE
        gb = self.gb
        Nd = self.Nd
        ss = self.scalar_scale
        ls = self.length_scale

        for g, d in gb:
            # Export pressure variable
            if self.scalar_variable in d[pp.STATE]:
                d[pp.STATE][self.p_exp] = d[pp.STATE][self.scalar_variable].copy() * ss
            else:
                d[pp.STATE][self.p_exp] = np.zeros((Nd, g.num_cells))

            if g.dim != 2:  # We only define tangential jumps in 2D fractures
                d[pp.STATE][self.normal_frac_u] = np.zeros(g.num_cells)
                d[pp.STATE][self.tangential_frac_u] = np.zeros(g.num_cells)

            if g.dim == Nd:  # On matrix
                u = d[pp.STATE][self.displacement_variable].reshape((Nd, -1), order='F').copy() * ls

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
        self.viz.write_vtk(data=self.export_fields, time_step=self.time)  # Write visualization
        self.export_times.append(self.time)

    def export_pvd(self):
        """ Implementation of export pvd"""
        self.viz.write_pvd(self.export_times)

    def initial_condition(self):
        """
        Initial guess for Newton iteration, scalar variable and bc_values (for time
        discretization).

        When stimulation phase is reached, we use displacements of last solution in
        initialize phase as initial condition for the cell displacements.
        """
        super().initial_condition()
        # TODO: Set hydrostatic initial condition

        # TODO: Scale variables
        if self.current_phase > 0:  # Stimulation phase

            for g, d in self.gb:
                if g.dim == self.Nd:
                    initial_displacements = d["initial_cell_displacements"]
                    pp.set_state(d, {self.displacement_variable: initial_displacements})

            for e, d in self.gb.edges():
                if e[0].dim == self.Nd:
                    try:
                        initial_displacements = d["initial_cell_displacements"]
                    except KeyError:
                        logger.warning("We got KeyError on d['initial_cell_displacements'].")
                        mg = d["mortar_grid"]
                        initial_displacements = np.zeros(mg.num_cells * self.Nd)
                    state = {
                        self.mortar_displacement_variable: initial_displacements,
                        "previous_iterate": {
                            self.mortar_displacement_variable: initial_displacements,
                        },
                    }
                    pp.set_state(d, state)

    @trace(logger)
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

    @trace(logger, timeit=False)
    def after_newton_convergence(self, solution, errors, iteration_counter):
        """ Overwrite from parent to export solution steps."""
        self.assembler.distribute_variable(solution)
        self.save_mechanical_bc_values()
        self.export_step()

    def after_simulation(self):
        """ Called after a time-dependent problem
        """
        self.export_pvd()
        logger.info(f"Solution exported to folder \n {self.viz_folder_name}")

    def after_newton_failure(self, solution, errors, iteration_counter):
        """ Instead of raising error on failure, save and return available data.
        """
        logger.error("Newton iterations did not converge")
        self.after_newton_convergence(solution, errors, iteration_counter)

        self.after_simulation()
        return self

    def prepare_initial_run(self):
        """
        Set time parameters for the preparation phase

        First, we run no flow for 6 hours to observe deformation due to mechanics itself.
        Then, [from Grimsel Experiment Description]:
        flow period is 40 minutes, followed by a shut-in period of 40 minutes.
        """

        # For the initialization phase, we use the following
        # start time
        self.time = - 6 * pp.HOUR
        # time step
        self.time_step = 3 * pp.HOUR
        # end time
        self.end_time = 0

        # self.time_step = 30 * pp.MINUTE

    def prepare_main_run(self):
        """ Adjust parameters between initial run and main run

        Total time: 80 minutes.
        Time step: 5 minutes
        """

        # New file name for this run
        self.file_name = 'main_run'
        self.set_viz()

        # We use the following time parameters
        # start time
        self.time = 0
        # time step
        self.time_step = 5 * pp.MINUTE
        # end time
        self.end_time = 40 * pp.MINUTE  # TODO: Change back to 40 minutes.

        # Store initial displacements
        for g, d in self.gb:
            if g.dim == 3:
                u = d[pp.STATE][self.displacement_variable]
                d["initial_cell_displacements"] = u

        for e, d in self.gb.edges():
            if e[0].dim == self.Nd:
                u = d[pp.STATE][self.mortar_displacement_variable]


                d["initial_cell_displacements"] = u

    def simulation_protocol(self):
        """ Adjust time step and other parameters for simulation protocol

                Here, we consider Doetsch et al (2018) [see e.g. p. 78/79 or App. J]
                Hydro Shearing Protocol:
                * Injection Cycle 3:
                    - Four injection steps of 10, 15, 20 and 25 l/min
                    - Each step lasts 10 minutes.
                    - Then, the interval is shut-in and monitored for 40 minutes.
                    - Venting was forseen at 20 minutes

                For this setup, we only consider Injection Cycle 3.

                Attributes set here:
                    current_phase : int
                        phase as a number (0 - 5)
                    current_injection_rate : float
                        fluid injection rate (l/min)
                """
        time_intervals = [
            # Phase 0: 0 l/min
            0,
            # Phase 1: 10 l/min
            10 * pp.MINUTE,
            # Phase 2: 15 l/min
            20 * pp.MINUTE,
            # Phase 3: 20 l/min
            30 * pp.MINUTE,
            # Phase 4: 25 l/min
            40 * pp.MINUTE,
            # Phase 5: 0 l/min
        ]

        injection_amount = [
            0,      # Phase 0
            10,     # Phase 1
            15,     # Phase 2
            20,     # Phase 3
            25,     # Phase 4
            0,      # Phase 5
        ]
        next_phase = np.searchsorted(time_intervals, self.time, side='right')
        if next_phase > self.current_phase:
            logger.info(f"A new phase has started: Phase {next_phase}")

        # Current phase number:
        self.current_phase = next_phase

        # Current injection amount [litres / second]
        self.current_injection_rate = injection_amount[self.current_phase] / pp.MINUTE

    @timer(logger)
    def prepare_simulation(self):
        """ Is run prior to a time-stepping scheme. Use this to initialize
        discretizations, linear solvers etc.


        ONLY CHANGE FROM PARENT:
        - Set self.viz with custom method.
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

    def check_convergence(self, solution, prev_solution, init_solution, nl_params=None):
        g_max = self._nd_grid()

        if not self._is_nonlinear_problem():
            # At least for the default direct solver, scipy.sparse.linalg.spsolve, no
            # error (but a warning) is raised for singular matrices, but a nan solution
            # is returned. We check for this.
            diverged = np.any(np.isnan(solution))
            converged = not diverged
            error = np.nan if diverged else 0
            return error, converged, diverged

        mech_dof = self.assembler.dof_ind(g_max, self.displacement_variable)
        scalar_dof = self.assembler.dof_ind(g_max, self.scalar_variable)

        # Also find indices for the contact variables
        contact_dof = np.array([], dtype=np.int)
        for e, _ in self.gb.edges():
            if e[0].dim == self.Nd:
                contact_dof = np.hstack(
                    (
                        contact_dof,
                        self.assembler.dof_ind(e[1], self.contact_traction_variable),
                    )
                )

        # Pick out the solution from current, previous iterates, as well as the
        # initial guess.
        u_mech_now = solution[mech_dof] * self.length_scale
        u_mech_prev = prev_solution[mech_dof] * self.length_scale
        u_mech_init = init_solution[mech_dof] * self.length_scale

        contact_now = solution[contact_dof] * self.scalar_scale * self.length_scale ** 2
        contact_prev = prev_solution[contact_dof] * self.scalar_scale * self.length_scale ** 2
        contact_init = init_solution[contact_dof] * self.scalar_scale * self.length_scale ** 2

        # Pressure solution
        p_scalar_now = solution[scalar_dof] * self.scalar_scale
        p_scalar_prev = prev_solution[scalar_dof] * self.scalar_scale
        p_scalar_init = init_solution[scalar_dof] * self.scalar_scale

        # Calculate errors

        # Displacement error
        difference_in_iterates_mech = np.sum((u_mech_now - u_mech_prev) ** 2)
        difference_from_init_mech = np.sum((u_mech_now - u_mech_init) ** 2)

        logger.info(f"diff iter u = {difference_in_iterates_mech:.6e}")
        logger.info(f"diff init u = {difference_from_init_mech:.6e}")

        # Contact traction error
        # TODO: Unsure about units of contact traction
        contact_norm = np.sum(contact_now ** 2)
        difference_in_iterates_contact = np.sum((contact_now - contact_prev) ** 2)
        difference_from_init_contact = np.sum((contact_now - contact_init) ** 2)

        logger.info(f"diff iter contact = {difference_in_iterates_contact:.6e}")
        logger.info(f"diff init contact = {difference_from_init_contact:.6e}")

        # Pressure scalar error
        scalar_norm = np.sum(p_scalar_now ** 2)
        difference_in_iterates_scalar = np.sum((p_scalar_now - p_scalar_prev) ** 2)
        difference_from_init_scalar = np.sum((p_scalar_now - p_scalar_init) ** 2)
        logger.info(f"diff iter scalar = {difference_in_iterates_scalar:.6e}")
        logger.info(f"diff init scalar = {difference_from_init_scalar:.6e}")

        tol_convergence = nl_params["nl_convergence_tol"]
        # Not sure how to use the divergence criterion
        # tol_divergence = nl_params["nl_divergence_tol"]

        converged = False
        diverged = False

        # Converge in displacement and pressure on 3D grid
        converged_u = False
        converged_p = False

        # Check absolute convergence criterion
        if difference_in_iterates_mech < tol_convergence:
            # converged = True
            converged_u = True
            error_mech = difference_in_iterates_mech
            logger.info(f"u converged absolutely.")

        else:
            # Check relative convergence criterion
            if (
                difference_in_iterates_mech
                < tol_convergence * difference_from_init_mech
            ):
                # converged = True
                converged_u = True
                logger.info(f"u converged relatively")
            error_mech = difference_in_iterates_mech / difference_from_init_mech

        # The if is intended to avoid division through zero
        if difference_in_iterates_contact < 1e-10:
            # converged = True
            error_contact = difference_in_iterates_contact
            logger.info(f"contact variable converged absolutely")
        else:
            error_contact = (
                difference_in_iterates_contact / difference_from_init_contact
            )

        # -- Scalar solution --
        # The if is intended to avoid division through zero
        if difference_in_iterates_scalar < tol_convergence:
            converged_p = True
            error_scalar = difference_in_iterates_scalar
            logger.info(f"pressure converged absolutely")
        else:
            # Relative convergence criterion:
            if difference_in_iterates_scalar < tol_convergence * difference_from_init_scalar:
                # converged = True
                converged_p = True
                logger.info(f"pressure converged relatively")

            error_scalar = (difference_in_iterates_scalar / difference_from_init_scalar)

        logger.info(f"Error in contact force is {error_contact:.6e}")
        logger.info(f"Error in matrix displacement is {error_mech:.6e}")
        logger.info(f"Error in pressure is {error_scalar:.6e}.")

        converged = converged_p and converged_u

        return error_mech, converged, diverged

    def _depth(self, coords):
        """
        Unscaled depth. We center the domain at 480m below the surface.
        (See Krietsch et al, 2018a)
        """
        return 480.0 * pp.METER - self.length_scale * coords[2]

