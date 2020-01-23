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
import pendulum

import porepy as pp
import numpy as np
import scipy.sparse as sps
from porepy.models.contact_mechanics_biot_model import ContactMechanicsBiot

import GTS as gts

logger = logging.getLogger(__name__)


def run_mechanics_model(
        *,
        viz_folder_name: str = None,
        result_file_name: str = None,
        isc_data_path: str = None,
        mesh_args: Mapping[str, int] = None,
        bounding_box: Mapping[str, int] = None,
        shearzone_names: List[str] = None,
        source_scalar_borehole_shearzone: Mapping[str, str] = None,
        scales: Mapping[str, float] = None,

):
    """ Send all initialization parameters to contact mechanics class

    Parameters
    ----------
    viz_folder_name : str
        Absolute path to folder where grid and results will be stored
        Default: /home/haakon/mastersproject/src/mastersproject/GTS/isc_modelling/results/default
    result_file_name : str
        Root name for simulation result files
        Default: 'main_run'
    isc_data_path : str
        Path to isc data: path/to/GTS/01BasicInputData
        Alternatively 'linux' or 'windows' for certain default paths (only applies to haakon's computers).
        Default: 'linux'
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
        Defaults to 1 for both.
    """
    # ------------------------------------------
    # --- FOLDER AND FILE RELATED PARAMETERS ---
    # ------------------------------------------

    # Set viz folder
    if viz_folder_name is None:
        viz_folder_name = (
            "/home/haakon/mastersproject/src/mastersproject/GTS/isc_modelling/results/default"
        )
    # Create viz folder path if it does not already exist
    if not os.path.exists(viz_folder_name):
        os.makedirs(viz_folder_name, exist_ok=True)

    # Set up logging
    logging.basicConfig(filename=viz_folder_name+"/results.log", level=logging.DEBUG)
    logger.info(f"Preparing setup for mechanics simulation on {pendulum.now().to_atom_string()}")
    logger.info(f"Visualization folder path: \n {viz_folder_name}")

    # Set file name of modelling results files
    if result_file_name is None:
        result_file_name = 'main_run'
    logger.info(f"Root file name of results: {result_file_name}")

    # Get data path to ISC data
    if isc_data_path is None:
        isc_data_path = 'linux'

    # ------------------------------------
    # --- MODELLING RELATED PARAMETERS ---
    # ------------------------------------

    # Set mesh arguments
    if mesh_args is None:
        mesh_size = 10
        mesh_args = {  # A very coarse grid
            "mesh_size_frac": mesh_size,
            "mesh_size_min": mesh_size,
            "mesh_size_bound": mesh_size,
        }
    logger.info(f"Mesh arguments: \n {mesh_args}")

    # Set bounding box
    if bounding_box is None:
        bounding_box = {
            "xmin": -6,
            "xmax": 80,
            "ymin": 55,
            "ymax": 150,
            "zmin": 0,
            "zmax": 50,
        }
    logger.info(f"Bounding box: \n {bounding_box}")

    # Set which shear-zones to include in simulation
    if shearzone_names is None:
        shearzone_names = ["S1_1", "S1_2", "S1_3", "S3_1", "S3_2"]
    logger.info(f"Shear zones in simulation: \n {shearzone_names}")

    # Set length scale and scalar scale
    if scales is None:
        scales = {
            'scalar_scale': 1,
            'length_scale': 1,
        }
    logger.info(f"Variable scaling: \n {scales}")

    # Set solver. 'pyamg' or 'direct'.
    solver = 'direct'
    logger.info(f"Solver type: {solver}")

    # TODO: Set custom time parameters with a class with only one method: _set_time_parameters

    # ---------------------------
    # --- PHYSICAL PARAMETERS ---
    # ---------------------------

    stress = stress_tensor()
    logger.info(f"Stress tensor: \n {stress}")
    # -------------------
    # --- SETUP MODEL ---
    # -------------------

    setup = gts.ContactMechanicsISC(
        viz_folder_name=viz_folder_name,
        result_file_name=result_file_name,
        isc_data_path=isc_data_path,
        mesh_args=mesh_args,
        bounding_box=bounding_box,
        shearzone_names=shearzone_names,
        source_scalar_borehole_shearzone=source_scalar_borehole_shearzone,
        scales=scales,
        stress=stress,
        solver=solver,
    )
    # -------------------------
    # --- SOLVE THE PROBLEM ---
    # -------------------------
    default_options = {  # Parameters for Newton solver.
        "max_iterations": 10,
        "convergence_tol": 1e-10,
        "divergence_tol": 1e5,
    }
    newton_options = default_options
    logger.info(f"Options for Newton solver: \n {newton_options}")

    logger.info("Setup complete. Starting simulation")
    pp.run_stationary_model(setup, params=newton_options)
    logger.info("Simulation complete. Exporting solution.")

    setup.export_step()
    logger.info(f"Solution exported to folder \n {viz_folder_name}")
    logger.info(f"Exits method on {pendulum.now().to_atom_string()}")
    return setup


def stress_tensor():
    """ Stress at ISC test site

    Values from Krietsch et al 2019
    """

    stress_value = np.array([13.1, 9.2, 8.7]) * pp.MEGA * pp.PASCAL
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