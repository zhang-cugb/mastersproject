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
from porepy.models.contact_mechanics_biot_model import ContactMechanicsBiot
from GTS.isc_modelling.contact_mechanics_biot import ContactMechanicsBiotISC
import pendulum

import GTS as gts
from util.logging_util import (
    __setup_logging,
    timer,
    trace,
)
from refinement import gb_coarse_fine_cell_mapping
from refinement.convergence import grid_error
import GTS.test.util as test_util


def test_unit_convergence_study():
    """ Unit test for convergence study
    """

    # 1. Prepare parameters
    stress = gts.isc_modelling.stress_tensor()
    hydrostatic = np.mean(np.diag(stress)) * np.ones(stress.shape[0])
    stress = np.diag(hydrostatic)

    params = {
        "mesh_args":
            {
                'mesh_size_frac': 10,
                'mesh_size_min': .1 * 10,
                'mesh_size_bound': 6 * 10,
            },
        "stress":
            stress,
        "shearzone_names":
            None,  # ['S1_1'],
    }

    n_refinements = 2

    # Run ContactMechanicsISC model
    gb_list = gts.isc_modelling.setup.run_models_for_convergence_study(
        model=gts.ContactMechanicsISC,
        params=params,
        n_refinements=n_refinements,
    )

    gb_ref = gb_list[-1]

    errors = []
    for i in range(0, n_refinements):
        gb_i = gb_list[i]
        gb_coarse_fine_cell_mapping(gb=gb_i, gb_ref=gb_ref)

        _error = grid_error(gb=gb_i, gb_ref=gb_ref, variable='u', variable_dof=3)
        errors.append(_error)

    return gb_list, errors







