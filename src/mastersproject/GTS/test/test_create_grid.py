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
import os
from pathlib import Path

import porepy as pp
import numpy as np
from porepy.models.contact_mechanics_model import ContactMechanics

import GTS as gts
from util.logging_util import __setup_logging
logger = logging.getLogger(__name__)


def test_create_grid():
    _this_file = Path(os.path.abspath(__file__)).parent
    _results_path = _this_file / "results"
    _results_path.mkdir(parents=True, exist_ok=True)  # Create path if not exists
    __setup_logging(_results_path)
    logger.info(f"Path to results: {_results_path}")

    # --- DOMAIN ARGUMENTS ---
    mesh_args = {'mesh_size_frac': 10, 'mesh_size_min': 10, 'mesh_size_bound': 10}
    bounding_box = {'xmin': -6, 'xmax': 80, 'ymin': 55, 'ymax': 150, 'zmin': 0, 'zmax': 50}

    shearzone_names = ["S1_1", "S1_2", "S1_3", "S3_1", "S3_2"]

    gb = gts.isc_modelling.create_isc_domain(
        _results_path,
        shearzone_names,
        bounding_box,
        mesh_args,
        n_refinements=0
    )[0]

    # Standard variables for ContactMechanicsISC class
    result_fname = 'test_run'
    isc_data_path = 'linux'

    scales = {
        'scalar_scale': 1,
        'length_scale': 1,
    }
    solver = 'direct'
    stress = gts.isc_modelling.stress_tensor()

    cm = gts.isc_modelling.ContactMechanicsISCWithGrid(
        _results_path,
        'test_run',
        'linux',
        mesh_args,
        bounding_box,
        shearzone_names,
        scales,
        stress,
        solver,
        gb,
    )

    return cm