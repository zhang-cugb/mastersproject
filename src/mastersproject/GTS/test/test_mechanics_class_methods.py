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
from util.logging_util import (
    __setup_logging,
    timer,
    trace,
)

logger = logging.getLogger(__name__)


def test_mechanical_boundary_conditions():
    """ Test that mechanical boundary conditions are correctly set
    """
    # Create setup with no fractures
    params = {'shearzone_names': None}
    setup = _prepare_tests(path_head="test_mechanical_boundary_conditions", params=params)

    gb: pp.GridBucket = setup.gb
    g: pp.Grid = gb.grids_of_dimension(3)[0]
    d = gb.node_props(g)

    mech_params: dict = d[pp.PARAMETERS][setup.mechanics_parameter_key]
    # Get the mechanical bc values
    mech_bc = mech_params['bc_values'].reshape((3, -1), order='F')
    mech_bc_type: pp.BoundaryConditionVectorial = mech_params['bc']

    all_bf, east, west, north, south, top, bottom = setup.domain_boundary_sides(g)
    # --- TESTS ---

    # 1. Test that mechanical BCs have the right type
    #   We want neumann on all faces except setup.faces_to_fix(),
    #   which should be on the bottom of the grid.
    dir_faces = setup.faces_to_fix(g)
    neu_faces = np.setdiff1d(all_bf, dir_faces)
    assert np.all(mech_bc_type.is_dir[:, dir_faces])
    assert np.all(mech_bc_type.is_neu[:, neu_faces])

    # 2.a. Test that Neumann BCs all point into the domain.
    #   The physical conditions at Grimsel Test Site implies
    #   a compressive boundary traction on all faces.

    # We expect the dirichlet faces to be on bottom. Remove them
    # for the sake of testing neumann conditions
    btm = np.setdiff1d(np.where(bottom)[0], dir_faces)

    # West
    assert np.all(mech_bc[0, west] > 0)
    # East
    assert np.all(mech_bc[0, east] < 0)
    # North
    assert np.all(mech_bc[1, north] < 0)
    # South
    assert np.all(mech_bc[1, south] > 0)
    # Top
    assert np.all(mech_bc[2, top] < 0)
    # Bottom
    assert np.all(mech_bc[2, btm] > 0)

    # 2.b. Test that Dirichlet BCs are zero
    np.allclose(mech_bc[:, dir_faces], 0)

    # 3. Test linear increase of mechanical boundary conditions (lithostatic conditions)


@trace(logger)
def test_decomposition_of_stress():
    """ Test the solutions acquired when decomposing stress to
    purely compressive and purely rotational components.

    Setup:
    1. Get the stress tensor and split in diagonal and off-diagonal components.
    2. Acquire solutions for each split tensor and the full tensor separately (3 cases).
    3. Compare solutions quantitatively (linearity of solution) and qualitatively (Paraview)

    TODO: Alternatively, separate the stress tensor to hydrostatic and deviatoric.
    """

    # Import stress tensor
    stress = gts.isc_modelling.stress_tensor()

    # Get normal and shear stresses
    normal_stress = np.diag(np.diag(stress).copy())
    shear_stress = stress - normal_stress

    # 1. Full stress tensor
    _folder_root = "test_decomposition_of_stress"
    params = {'shearzone_names': None, 'stress': normal_stress}
    setup = _prepare_tests(f"{_folder_root}", params=params, prepare_simulation=False)
    setup.create_grid()

    # -- Safely copy the grid generated above --
    # Ensures the implicit in-place variable changes across grid_buckets (gb.copy() is insufficient)
    # Get the grid_bucket .msh path
    path_to_gb_msh = f"{setup.viz_folder_name}/gmsh_frac_file.msh"
    # Re-create the grid buckets
    gb_n = pp.fracture_importer.dfm_from_gmsh(path_to_gb_msh, dim=3, network=setup._network)
    gb_s = pp.fracture_importer.dfm_from_gmsh(path_to_gb_msh, dim=3, network=setup._network)

    # 1. Pure normal stress
    params = {'shearzone_names': None, 'stress': normal_stress}
    setup_n = _prepare_tests(f"{_folder_root}/normal_stress", params=params, prepare_simulation=False)
    setup_n.set_grid(gb_n)

    # 2. Pure shear stress
    params = {'shearzone_names': None, 'stress': shear_stress}
    setup_s = _prepare_tests(f"{_folder_root}/shear_stress", params=params, prepare_simulation=False)
    setup_s.set_grid(gb_s)

    # --- Run simulations ---
    params = None  # Use default Newton solver parameters

    # 1. Full stress tensor
    pp.run_stationary_model(setup, params)

    # 2. Pure normal stress
    pp.run_stationary_model(setup_n, params)

    # 3. Pure shear stress
    pp.run_stationary_model(setup_s, params)

    # --- Compare results ---




def _prepare_tests(path_head, params=None, prepare_simulation=True):
    """ Helper method to create grids, etc. for test methods
    """
    if params is None:
        params = {}

    _this_file = Path(os.path.abspath(__file__)).parent
    _results_path = _this_file / f"results/test_mechanics_class_methods/{path_head}"
    _results_path.mkdir(parents=True, exist_ok=True)  # Create path if not exists
    __setup_logging(_results_path)
    logger.info(f"Path to results: {_results_path}")

    # --- DOMAIN ARGUMENTS ---
    default_params = {
        'mesh_args':
            {'mesh_size_frac': 10, 'mesh_size_min': .1 * 10, 'mesh_size_bound': 6 * 10},
        'bounding_box':
            {'xmin': -20, 'xmax': 80, 'ymin': 50, 'ymax': 150, 'zmin': -25, 'zmax': 75},
        'shearzone_names':
            ["S1_1", "S1_2", "S1_3", "S3_1", "S3_2"],
        'folder_name':
            _results_path,
        'solver':
            'direct',
        'stress':
            gts.isc_modelling.stress_tensor(),
        'length_scale':
            1,
        'scalar_scale':
            1,
    }
    default_params.update(params)

    setup = gts.ContactMechanicsISC(default_params)
    if prepare_simulation:
        setup.prepare_simulation()
    return setup


