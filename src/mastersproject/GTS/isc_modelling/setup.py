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


def setup_model(
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
    """ Send all initialization parameters to contact mechanics or contact mechanics biot class.

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
    logging.info(f"Visualization folder path: {viz_folder_name}")

    # Set file name of modelling results files
    if result_file_name is None:
        result_file_name = 'main_run'

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

    # Set which shear-zones to include in simulation
    if shearzone_names is None:
        shearzone_names = ["S1_1", "S1_2", "S1_3", "S3_1", "S3_2"]

    # Set borehole and shear-zone names for which injection occurs.
    if source_scalar_borehole_shearzone is None:
        source_scalar_borehole_shearzone = {
            "shearzone": "S1_1",
            "borehole": "INJ1",
        }

    # Set length scale and scalar scale
    if scales is None:
        scales = {
            'scalar_scale': 1,
            'length_scale': 1,
        }

    # TODO: Set custom time parameters with a class with only one method: _set_time_parameters

