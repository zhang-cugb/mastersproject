import os
import logging
from typing import (  # noqa
    Any,
    Callable,
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
from pathlib import Path
from pprint import pformat
import functools
import inspect

import pendulum
import porepy as pp
import numpy as np
import scipy.sparse as sps
from porepy.models.contact_mechanics_biot_model import ContactMechanicsBiot
from porepy.models.contact_mechanics_model import ContactMechanics

import GTS as gts
from refinement import refine_mesh

# --- LOGGING UTIL ---
from util.logging_util import (
    timer,
    trace,
    __setup_logging,
)

logger = logging.getLogger(__name__)


@trace(logger)
def run_biot_model(
        *,
        viz_folder_name: str = None,
        mesh_args: Mapping[str, int] = None,
        bounding_box: Mapping[str, int] = None,
        shearzone_names: List[str] = None,
        source_scalar_borehole_shearzone: Mapping[str, str] = None,
        length_scale: float = None,
        scalar_scale: float = None,
):
    """ Send all initialization parameters to contact mechanics biot class

        Parameters
        ----------
        viz_folder_name : str
            Absolute path to folder where grid and results will be stored
            Default: /home/haakon/mastersproject/src/mastersproject/GTS/isc_modelling/results/default
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
        length_scale, scalar_scale : float : Optional
            Length scale and scalar variable scale.
        """
    params = {
        "folder_name": viz_folder_name,
        "mesh_args": mesh_args,
        "bounding_box": bounding_box,
        "shearzone_names": shearzone_names,
        "source_scalar_borehole_shearzone": source_scalar_borehole_shearzone,
        "length_scale": length_scale,
        "scalar_scale": scalar_scale,
    }

    setup = run_abstract_model(
        model=gts.ContactMechanicsBiotISC,
        run_model_method=pp.run_time_dependent_model,
        params=params,
    )

    return setup


@trace(logger)
def run_biot_gts_model(params):
    """ Set up and run biot model with
    an initialization run and a main run.
    """

    setup = run_abstract_model(
        model=gts.ContactMechanicsBiotISC,
        run_model_method=gts_biot_model,
        params=params,
    )

    return setup


@trace(logger)
def run_mechanics_model(
        *,
        viz_folder_name: str = None,
        mesh_args: Mapping[str, int] = None,
        bounding_box: Mapping[str, int] = None,
        shearzone_names: List[str] = None,
        length_scale: float = None,
        scalar_scale: float = None,
):
    """ Send all initialization parameters to contact mechanics class

    Parameters
    ----------
    viz_folder_name : str
        Absolute path to folder where grid and results will be stored
        Default: /home/haakon/mastersproject/src/mastersproject/GTS/isc_modelling/results/default
    mesh_args : Mapping[str, int]
        Arguments for meshing of domain.
        Required keys: 'mesh_size_frac', 'mesh_size_min, 'mesh_size_bound'
    bounding_box : Mapping[str, int]
        Bounding box of domain
        Required keys: 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'.
    shearzone_names : List[str]
        Which shear-zones to include in simulation
    length_scale, scalar_scale : float : Optional
            Length scale and scalar variable scale.
    """
    params = {
        "folder_name": viz_folder_name,
        "mesh_args": mesh_args,
        "bounding_box": bounding_box,
        "shearzone_names": shearzone_names,
        "length_scale": length_scale,
        "scalar_scale": scalar_scale,
    }

    setup = run_abstract_model(
        model=gts.ContactMechanicsISC,
        run_model_method=pp.run_stationary_model,
        params=params,
    )

    return setup


def gts_biot_model(setup, params):
    """ Setup for time-dependent model run at Grimsel Test Site

    Usually called by run_abstract_model if this method is supplied
    as the argument 'run_model'.
    """

    # Initialization phase
    pp.run_time_dependent_model(setup=setup, params=params)
    logger.info(
        f"Initial simulation complete. Exporting solution. Time: {pendulum.now().to_atom_string()}"
    )
    # Stimulation phase
    logger.info(f"Starting stimulation phase at time: {pendulum.now().to_atom_string()}")
    setup.prepare_main_run()
    logger.info("Setup complete. Starting time-dependent simulation")
    pp.run_time_dependent_model(setup=setup, params=params)


def run_abstract_model(
        model: Type[ContactMechanics],
        run_model_method: Callable,
        params: dict = None,
        newton_params: dict = None,
):
    """ Set up and run an abstract model

    Parameters
    ----------
    model : Type[ContactMechanics]
        Which model to run
        Only tested for subclasses of ContactMechanicsISC
    run_model_method : Callable
        Which method to run model with
        Typically pp.run_stationary_model or pp.run_time_dependent_model
    params : dict (Default: None)
        Any non-default parameters to use
    newton_params : dict (Default: None)
        Any non-default newton solver parameters to use
    """

    setup = _abstract_model_setup(
        model=model,
        params=params,
    )
    # -------------------------
    # --- SOLVE THE PROBLEM ---
    # -------------------------
    default_options = {  # Parameters for Newton solver.
        "max_iterations": 20,
        "nl_convergence_tol": 1e-6,
        "nl_divergence_tol": 1e5,
    }
    newton_options = default_options.update(newton_params)
    logger.info(f"Options for Newton solver: \n {newton_options}")
    logger.info("Setup complete. Starting simulation")

    run_model_method(setup=setup, params=default_options)

    logger.info(f"Simulation complete. Exporting solution. Time: {pendulum.now().to_atom_string()}")

    return setup


def _abstract_model_setup(
        model: Type[ContactMechanics],
        params: dict = None,
):
    """ Helper method to assemble model setup for biot and mechanics.

    Parameters
    ----------
    model : Type[ContactMechanics]
        Which model to run
        Accepted: 'ContactMechanicsISC', 'ContactMechanicsBiotISC'
    params : dict (Default: None)
        Custom parameters to pass to model
        See below of default values
    """
    # ------------------------------------------
    # --- FOLDER AND FILE RELATED PARAMETERS ---
    # ------------------------------------------
    _this_file = Path(os.path.abspath(__file__)).parent
    path_head = "default/default_1"
    _results_path = _this_file / f"results/{path_head}"

    # ------------------------------------
    # --- MODELLING RELATED PARAMETERS ---
    # ------------------------------------
    sz = 10
    mesh_args = {
        'mesh_size_frac': sz,
        'mesh_size_min': 0.1 * sz,
        'mesh_size_bound': 6 * sz,
    }

    bounding_box = {
        'xmin': -20,
        'xmax': 80,
        'ymin': 50,
        'ymax': 150,
        'zmin': -25,
        'zmax': 75
    }

    # Set which borehole / shearzone to inject fluid to
    # This corresponds to setup in HS2 from Doetsch et al 2018
    source_scalar_borehole_shearzone = {
        "shearzone": "S1_2",
        "borehole": "INJ1",
    }  # (If ContactMechanicsISC is run, this is ignored)

    # ---------------------------
    # --- PHYSICAL PARAMETERS ---
    # ---------------------------

    stress = stress_tensor()

    # -----------------------------
    # --- INITIALIZE PARAMETERS ---
    # -----------------------------

    in_params = {
        "folder_name":
            _results_path,
        "mesh_args":
            mesh_args,
        "bounding_box":
            bounding_box,
        "shearzone_names":
            None,
        "length_scale":
            1,
        "scalar_scale":
            1,
        "solver":
            "direct",
        "source_scalar_borehole_shearzone":
            source_scalar_borehole_shearzone,
        "stress":
            stress,
    }

    # Update default parameter set with input parameters
    in_params.update(params)

    # ---------------------
    # --- BACKEND SETUP ---
    # ---------------------

    # Create viz folder path if it does not already exist
    viz_folder_name = in_params["folder_name"]
    Path(viz_folder_name).mkdir(parents=True, exist_ok=True)

    # Set up logging
    __setup_logging(viz_folder_name)
    logger.info(f"Preparing setup for mechanics simulation on {pendulum.now().to_atom_string()}")
    logger.info(f"Simulation parameters:\n {pformat(in_params)}")

    # -------------------
    # --- SETUP MODEL ---
    # -------------------

    setup = model(params=params)
    return setup


def create_isc_domain(
        viz_folder_name: Union[str, Path],
        shearzone_names: List[str],
        bounding_box: dict,
        mesh_args: dict,
        n_refinements: int = 0):
    """ Create a domain (.geo file) for the ISC test site.

    Parameters
    ----------
    viz_folder_name : str or pathlib.Path
        Absolute path to folder to store results in
    shearzone_names : List of str
        Names of shearzones to include
    bounding_box : dict
        Bounding box of domain ('xmin', 'xmax', etc.)
    mesh_args : dict
        Arguments for meshing (of coarsest grid)
    n_refinements : int, Default = 0
        Number of refined grids to produce.
        The grid is refined by splitting.
    """

    # ----------------------------------------
    # --- CREATE FRACTURE NETWORK AND MESH ---
    # ----------------------------------------
    network = gts.fracture_network(
        shearzone_names=shearzone_names,
        export_vtk=True,
        domain=bounding_box,
        network_path=f"{viz_folder_name}/fracture_network.vtu"
    )

    gmsh_file_name = str(viz_folder_name / "gmsh_frac_file")
    gb = network.mesh(mesh_args=mesh_args, file_name=gmsh_file_name)

    gb_list = refine_mesh(
        in_file=f'{gmsh_file_name}.geo',
        out_file=f"{gmsh_file_name}.msh",
        dim=3,
        network=network,
        num_refinements=n_refinements,
    )

    # TODO: Make this procedure "safe".
    #   E.g. assign names by comparing normal vector and centroid.
    #   Currently, we assume that fracture order is preserved in creation process.
    #   This may be untrue if fractures are (completely) split in the process.
    # Assign node prop 'name' to each grid in the grid bucket.
    for _gb in gb_list:
        pp.contact_conditions.set_projections(_gb)
        _gb.add_node_props(keys="name")
        fracture_grids = _gb.get_grids(lambda g: g.dim == _gb.dim_max() - 1)
        if shearzone_names is not None:
            for i, sz_name in enumerate(shearzone_names):
                _gb.set_node_prop(fracture_grids[i], key="name", val=sz_name)
                # Note: Use self.gb.node_props(g, 'name') to get value.

    return gb_list


@timer(logger)
def convergence_study():
    """ Perform a convergence study of a given problem setup.
    """

    # 1. Step: Create n grids by uniform refinement.
    # 2. Step: for grid i in list of n grids:
    # 2. a. Step: Set up the mechanics model.
    # 2. b. Step: Solve the mechanics problem.
    # 2. c. Step: Keep the grid (with solution data)
    # 3. Step: Let the finest grid be the reference solution.
    # 4. Step: For every other grid:
    # 4. a. Step: Map the solution to the fine grid, and compute error.
    # 5. Step: Compute order of convergence, etc.

    # -----------------
    # --- ARGUMENTS ---
    # -----------------
    viz_folder_name = Path(os.path.abspath(__file__)).parent / "results/mech_convergence_2test"
    if not os.path.exists(viz_folder_name):
        os.makedirs(viz_folder_name, exist_ok=True)

    shearzone_names = ["S1_1", "S1_2", "S1_3", "S3_1", "S3_2"]

    mesh_size = 10
    mesh_args = {  # A very coarse grid
        "mesh_size_frac": mesh_size,
        "mesh_size_min": mesh_size,
        "mesh_size_bound": mesh_size,
    }

    bounding_box = {
        "xmin": -6,
        "xmax": 80,
        "ymin": 55,
        "ymax": 150,
        "zmin": 0,
        "zmax": 50,
    }

    # 1. Step: Create n grids by uniform refinement.
    gb_list = create_isc_domain(
        viz_folder_name=viz_folder_name,
        shearzone_names=shearzone_names,
        bounding_box=bounding_box,
        mesh_args=mesh_args,
        n_refinements=1,
    )

    scales = {
        'scalar_scale': 1,
        'length_scale': 1,
    }
    solver = 'direct'

    # ---------------------------
    # --- PHYSICAL PARAMETERS ---
    # ---------------------------

    stress = stress_tensor()

    # ----------------------
    # --- SET UP LOGGING ---
    # ----------------------
    print(viz_folder_name / "results.log")
    logger = __setup_logging(viz_folder_name)

    logger.info(f"Preparing setup for mechanics convergence study on {pendulum.now().to_atom_string()}")
    logger.info(f"Reporting on {len(gb_list)} grid buckets.")
    logger.info(f"Visualization folder path: \n {viz_folder_name}")
    logger.info(f"Mesh arguments for coarsest grid: \n {mesh_args}")
    logger.info(f"Bounding box: \n {bounding_box}")
    logger.info(f"Variable scaling: \n {scales}")
    logger.info(f"Solver type: {solver}")
    logger.info(f"Stress tensor: \n {stress}")

    # -----------------------
    # --- SETUP AND SOLVE ---
    # -----------------------

    newton_options = {  # Parameters for Newton solver.
        "max_iterations": 10,
        "convergence_tol": 1e-10,
        "divergence_tol": 1e5,
    }
    logger.info(f"Options for Newton solver: \n {newton_options}")

    from GTS.isc_modelling.mechanics import ContactMechanicsISCWithGrid
    for gb in gb_list:
        setup = ContactMechanicsISCWithGrid(
            viz_folder_name, 'main_run', 'linux', mesh_args, bounding_box,
            shearzone_names, scales, stress, solver, gb,
        )

        logger.info("Setup complete. Starting simulation")
        pp.run_stationary_model(setup, params=newton_options)
        logger.info("Simulation complete. Exporting solution.")

    return gb_list


def stress_tensor():
    """ Stress at ISC test site

    Values from Krietsch et al 2019
    """

    # Note: Negative side due to compressive stresses
    stress_value = - np.array([13.1, 9.2, 8.7]) * pp.MEGA * pp.PASCAL
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
