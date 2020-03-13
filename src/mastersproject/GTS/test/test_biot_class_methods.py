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
from porepy.models.contact_mechanics_biot_model import ContactMechanicsBiot
from GTS.isc_modelling.contact_mechanics_biot import ContactMechanicsBiotISC
import pendulum

import GTS as gts
from src.mastersproject.util.logging_util import (
    __setup_logging,
    timer,
    trace,
)
import GTS.test.util as test_util

logger = logging.getLogger(__name__)


def test_compare_run_mech_and_run_mech_by_filter_term():
    """ This test runs pure mechanics by running the test
    'test_mechanics_class_methods.test_decomposition_of_stress()'
    for the hydrostatic case. Then, it runs the test
    'test_run_mechanics_term_by_filter()'.

    The goal is to compare the output of these two methods and ensure they
    are the same.
    """
    # 1. Prepare parameters
    stress = gts.isc_modelling.stress_tensor()
    # We set up hydrostatic stress
    hydrostatic = np.mean(np.diag(stress)) * np.ones(stress.shape[0])
    stress = np.diag(hydrostatic)
    no_shearzones = None
    gravity = False  # No gravity effects
    params = {
        "stress": stress,
        "shearzone_names": no_shearzones,
        "_gravity_bc": gravity,
        "_gravity_src": gravity,
    }

    # Storage folder
    this_method_name = test_compare_run_mech_and_run_mech_by_filter_term.__name__
    now_as_YYMMDD = pendulum.now().format("YYMMDD")
    _folder_root = f"{this_method_name}/{now_as_YYMMDD}"

    # 1. --- Setup ContactMechanicsISC ---
    setup_mech = test_util.prepare_setup(
        model=gts.ContactMechanicsISC,
        path_head=f"{_folder_root}/test_mech",
        params=params,
        prepare_simulation=False,
        setup_loggers=True,
    )
    setup_mech.create_grid()

    # 2. --- Setup BiotReduceToMechanics ---
    params2 = test_util.prepare_params(
        path_head=f"{_folder_root}/test_biot_reduce_to_mech",
        params=params,
        setup_loggers=False
    )
    setup_biot = BiotReduceToMechanics(params=params2)

    # Recreate the same mesh as for the above setup
    path_to_gb_msh = f"{setup_mech.viz_folder_name}/gmsh_frac_file.msh"
    gb2 = pp.fracture_importer.dfm_from_gmsh(path_to_gb_msh, dim=3, network=setup_mech._network)
    setup_biot.set_grid(gb2)

    # 3. --- Run simulations ---
    nl_params = {}

    # Run ContactMechanicsISC
    pp.run_stationary_model(setup_mech, nl_params)

    # Run BiotReduceToMechanics
    pp.run_time_dependent_model(setup_biot, nl_params)

    # --- Compare results ---
    def get_u(_setup):
        gb = _setup.gb
        g = gb.grids_of_dimension(3)[0]
        d = gb.node_props(g)
        u = d['state']['u'].reshape((3, -1), order='F')
        return u

    u_mech = get_u(setup_mech)
    u_biot = get_u(setup_biot)
    assert np.isclose(np.sum(np.abs(u_mech - u_biot)), 0.0), "Running mechanics or biot (only discretize mechanics " \
                                                             "term should return same result."
    return setup_mech, setup_biot


def test_run_mechanics_term_by_filter():
    """ This test intends to replicate the part of the
     results of 'test_decomposition_of_stress' by only
     discretizing the mechanics term.
    """

    # 1. Prepare parameters
    stress = gts.isc_modelling.stress_tensor()
    # We set up hydrostatic stress
    hydrostatic = np.mean(np.diag(stress)) * np.ones(stress.shape[0])
    stress = np.diag(hydrostatic)

    no_shearzones = None
    gravity = False  # No gravity effects

    params = {
        "stress": stress,
        "shearzone_names": no_shearzones,
        "_gravity_bc": gravity,
        "_gravity_src": gravity,
    }

    # Storage folder
    this_method_name = test_run_mechanics_term_by_filter.__name__
    now_as_YYMMDD = pendulum.now().format("YYMMDD")
    _folder_root = f"{this_method_name}/{now_as_YYMMDD}/test_1"

    #
    # 2. Setup and run test
    params = test_util.prepare_params(path_head=_folder_root, params=params, setup_loggers=True)
    setup = BiotReduceToMechanics(params=params)

    nl_params = {}  # Default Newton Iteration parameters
    pp.run_time_dependent_model(setup, params=nl_params)

    return setup


# --- Overwrite discretize method in ContactMechanics ---
# We overwrite the discretize() method to discretize only
# the desired terms.
class BiotReduceToMechanics(ContactMechanicsBiotISC):
    @trace(logger)
    def discretize(self):
        """ Discretize the mechanics stress term
        """
        if not hasattr(self, "assembler"):
            self.assembler = pp.Assembler(self.gb, active_variables=[self.displacement_variable])

        g_max = self.gb.grids_of_dimension(self.Nd)[0]
        self.assembler.discretize(grid=g_max, variable_filter=self.displacement_variable)


