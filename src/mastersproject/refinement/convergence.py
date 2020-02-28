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

# --- LOGGING UTIL ---
from util.logging_util import timer, trace

logger = logging.getLogger(__name__)


def _grid_error(mappings, gb: pp.GridBucket, gb_ref: pp.GridBucket, variable, variable_dof):
    """ Compute grid errors for a given mapping

    The mapping is the one computed by gb_coarse_fine_cell_mapping()
    """
    # TODO: Fix this method. Something is wrong when passing a mapping to it.
    errors = {}

    for g, g_ref, mapping in mappings:
        # print(i, pair)
        # g = mappings[i][0]
        # g_ref = mappings[i][1]
        # mapping = mappings[i][2]

        assert g.num_cells < g_ref.num_cells

        data = gb.node_props(g)
        data_ref = gb_ref.node_props(g_ref)

        errors[data['node_number']] = {}  # Initialize this dict entry

        states = data[pp.STATE]
        states_ref = data_ref[pp.STATE]

        # TODO: Add some limitation to which keys you want to check,
        #  or how you should compute errors over certain types of keys
        state_keys = set(states.keys())
        state_ref_keys = set(states_ref.keys())
        check_keys = state_keys.intersection(state_ref_keys)

        if variable not in check_keys:
            logger.info(f"{variable} not present on grid number "
                        f"{gb.node_props(g, 'node_number')} of dim {g.dim}.")

        sol = states[variable].reshape((-1, variable_dof))
        mapped_sol = mapping.dot(sol).reshape((-1, 1))
        sol_ref = states_ref[variable]

        absolute_error = np.linalg.norm(mapped_sol - sol_ref)

        norm_ref = np.linalg.norm(sol_ref)
        if norm_ref < 1e-5:
            logger.warning(f"Relative error not reportable. "
                           f"Norm of reference solution is {norm_ref}. "
                           f"Reporting absolute error")
            relative_error = -1

        relative_error = absolute_error / norm_ref

        errors[data['node_number']] = {variable: {'absolute_error': absolute_error,
                                                  'relative_error': relative_error}}

    return errors







