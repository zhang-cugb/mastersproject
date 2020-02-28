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


def grid_error(
        gb: pp.GridBucket,
        gb_ref: pp.GridBucket,
        variable: str,
        variable_dof: int,
) -> dict:
    """ Compute grid errors a grid bucket and refined reference grid bucket

    Assumes that the coarse grid bucket has a node property
    'coarse_fine_cell_mapping' assigned on each grid, which
    maps from coarse to fine cells according to the method
    'coarse_fine_cell_mapping(...)'.

    Parameters
    ----------
    gb, gb_ref : pp.GridBucket
        Coarse and fine grid buckets, respectively
    variable : str
        which variable to compute error over
    variable_dof : int
        Degrees of freedom of 'variable'.

    Returns
    -------
    errors : dict
        Dictionary with top level keys as node_number,
        within which for each variable, the error is
        reported.
    """
    errors = {}

    grids = gb.get_grids()
    grids_ref = gb_ref.get_grids()
    n_grids = len(grids)

    for i in np.arange(n_grids):
        g, g_ref = grids[i], grids_ref[i]
        mapping = gb.node_props(g, "coarse_fine_cell_mapping")

        # Get states
        data = gb.node_props(g)
        data_ref = gb_ref.node_props(g_ref)
        states = data[pp.STATE]
        states_ref = data_ref[pp.STATE]

        # TODO: Compute error over a list of variables and variable_dofs.
        # Check if the variable exists on both the grid and reference grid
        state_keys = set(states.keys())
        state_ref_keys = set(states_ref.keys())
        check_keys = state_keys.intersection(state_ref_keys)
        if variable not in check_keys:
            logger.info(f"{variable} not present on grid number "
                        f"{gb.node_props(g, 'node_number')} of dim {g.dim}.")
            continue

        # Compute errors relative to the reference grid
        # TODO: Should the solution be divided by g.cell_volumes or similar?
        # TODO: If scaling is used, consider that - or use the export-ready variables,
        #   'u_exp', 'p_exp', etc.
        sol = states[variable].reshape((variable_dof, -1), order='F').T  # (num_cells x variable_dof)
        mapped_sol: np.ndarray = mapping.dot(sol).ravel(order='F')
        sol_ref = states_ref[variable]

        absolute_error = np.linalg.norm(mapped_sol - sol_ref)
        norm_ref = np.linalg.norm(sol_ref)

        if norm_ref < 1e-5:
            logger.warning(f"Relative error not reportable. "
                           f"Norm of reference solution is {norm_ref}. "
                           f"Reporting absolute error")
            error = absolute_error
            is_relative = False
        else:
            error = absolute_error / norm_ref
            is_relative = True

        errors[data["node_number"]] = {
            variable: {
                "error":
                    error,
                "is_relative":
                    is_relative,
            }
        }

    return errors



