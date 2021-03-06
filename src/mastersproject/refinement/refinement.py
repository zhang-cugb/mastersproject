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


@trace(logger)
def refine_mesh(
        in_file: str, out_file: str, dim: int,
        network: Union[pp.FractureNetwork3d, pp.FractureNetwork2d],
        num_refinements: int = 1,
) -> List[pp.GridBucket]:
    """ Refine a mesh by splitting, using gmsh

    Parameters
    ----------
    in_file : str
        path to .geo file to read
    out_file : str
        path to new .msh file to store mesh in, excluding the ending '.msh'.
    dim : int {2, 3}
        Dimension of domain to mesh
    network : Union[pp.FractureNetwork2d, pp.FractureNetwork3d]
        PorePy class defining the fracture network that is described by the .geo in_file
    num_refinements : int : Optional. Default = 1
        Number of refinements
    """

    try:
        import gmsh
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "To run gmsh python api on your system, "
            "download the relevant gmsh*-sdk.* from http://gmsh.info/bin/. "
            "Then, Add the 'lib' directory from the SDK to PYTHONPATH: \n"
            "export PYTHONPATH=${PYTHONPATH}:path/to/gmsh*-sdk.*/lib"
        )

    from porepy.fracs.simplex import tetrahedral_grid_from_gmsh
    from porepy.fracs.meshing import grid_list_to_grid_bucket

    assert os.path.isfile(in_file)

    # Run gmsh
    gmsh.initialize()
    gmsh.open(in_file)
    gmsh.model.mesh.generate(dim=dim)
    if out_file[-4:] == ".msh":
        out_file = out_file[:-4]

    # Save coarsest grid
    fname = f"{out_file}_0.msh"
    gmsh.write(fname)
    grid_list_ref = tetrahedral_grid_from_gmsh(network=network, file_name=fname)
    gb_list = [grid_list_to_grid_bucket(grid_list_ref)]

    for i in range(num_refinements):
        gmsh.model.mesh.refine()  # Refined grid

        fname = f"{out_file}_{i+1}.msh"
        gmsh.write(fname)

        # Create grid bucket from refined grid output
        grid_list_ref = tetrahedral_grid_from_gmsh(network=network, file_name=fname)
        gb_ref = grid_list_to_grid_bucket(grid_list_ref)
        gb_list.append(gb_ref.copy())

    gmsh.finalize()
    return gb_list


@trace(logger)
def gb_coarse_fine_cell_mapping(
        gb: pp.GridBucket, gb_ref: pp.GridBucket, tol=1e-8
):
    """ Wrapper for coarse_fine_cell_mapping to construct mapping for grids in GridBucket.

    Parameters
    ----------
    gb : pp.GridBucket
        Coarse grid bucket
    gb_ref : pp.GridBucket
        Refined grid bucket
    tol : float, Optional
        Tolerance for point_in_poly* -methods

    Returns
    -------
    mapping : list of tuples with entries (pp.GridBucket, pp.GridBucket, sps.csc_matrix)
        The first entry is the coarse grid.
        The second entry is the refined grid.
        The third entry is the mapping from coarse to fine cells
    """

    grids = gb.get_grids()
    grids_ref = gb_ref.get_grids()

    assert len(grids) == len(grids_ref), "Weakly check that GridBuckets refer to same domains"
    assert np.array_equal(np.append(*gb.bounding_box()), np.append(*gb_ref.bounding_box())), \
        "Weakly check that GridBuckets refer to same domains"

    # This method assumes a consistent node ordering between grids. At least assign one.
    gb.assign_node_ordering(overwrite_existing=False)
    gb_ref.assign_node_ordering(overwrite_existing=False)

    n_grids = len(grids)
    # mappings = [None]*n_grids
    mappings = {'gb': gb, 'gb_ref': gb_ref}

    for i in np.arange(n_grids):
        g, g_ref = grids[i], grids_ref[i]
        node_num, node_num_ref = gb._nodes[g]['node_number'], gb_ref._nodes[g_ref]['node_number']

        assert node_num == node_num_ref, "Weakly check that grids refer to same domain."

        mapping = coarse_fine_cell_mapping(g, g_ref, tol=tol)

        mappings[(g, g_ref)] = {'node_number': node_num,
                                'data': gb.node_props(g),
                                'data_ref': gb_ref.node_props(g_ref)}

    return mappings


@trace(logger)
def coarse_fine_cell_mapping(g: pp.Grid, g_ref: pp.Grid, tol=1e-8):
    """ Construct a mapping between cells of a grid and its refined version

    Assuming a regular and a refined mesh, where the refinement is executed by splitting.
    I.e. a cell in the refined grid is completely contained within a cell in the
    coarse grid.

    Parameters
    ----------
    g : pp.Grid
        Coarse grid
    g_ref : pp.Grid
        Refined grid
    tol : float, Optional
        Tolerance for pp.geometry_property_checks.point_in_polyhedron()

    Returns
    -------
    coarse_fine : sps.csc_matrix
        Column major sparse matrix mapping from coarse to fine cells.
    """

    assert g.num_cells < g_ref.num_cells, "Wrong order of input grids"
    assert g.dim == g_ref.dim, "Grids must be of same dimension"

    cell_nodes = g.cell_nodes()
    slices = zip(cell_nodes.indptr[:-1], cell_nodes.indptr[1:])  # start/end row pointers for each column

    # Create sps.csc_matrix mapping coarse cells to fine cell centers
    indptr = np.array([0])
    indices = np.empty(0)

    cells_ref = g_ref.cell_centers.copy()  # Cell centers in fine grid
    test_cells_ptr = np.arange(g_ref.num_cells)  # Pointer to cell centers
    nodes = g.nodes.copy()

    if g.dim == 1:
        nodes = nodes.copy()
        tangent = pp.map_geometry.compute_tangent(nodes)
        reference = [1, 0, 0]
        R = pp.map_geometry.project_line_matrix(nodes, tangent, tol=tol, reference=reference)
        nodes = R.dot(nodes)[0, :]
        cells_ref = R.dot(cells_ref)[0, :]

    elif g.dim == 2:  # Pre-processing for efficiency
        nodes = nodes.copy()
        R = pp.map_geometry.project_plane_matrix(nodes, check_planar=False)
        nodes = np.dot(R, nodes)[:2, :]
        cells_ref = np.dot(R, cells_ref)[:2, :]

    # Loop through every coarse cell
    for st, nd in slices:

        nodes_idx = cell_nodes.indices[st:nd]
        num_nodes = nodes_idx.size

        if g.dim == 1:
            assert (num_nodes == 2)
            line = np.sort(nodes[nodes_idx])
            test_points = cells_ref[test_cells_ptr]
            in_poly = np.searchsorted(line, test_points, side='left') == 1

        elif g.dim == 2:
            assert (num_nodes == 3), "We assume simplexes in 2D (i.e. 3 nodes)"
            polygon = nodes[:, nodes_idx]
            test_points = cells_ref[:, test_cells_ptr]
            in_poly = pp.geometry_property_checks.point_in_polygon(
                polygon, test_points, tol=tol)

        elif g.dim == 3:
            # Make polyhedron from node coordinates
            # Polyhedron defined as a list of nodes defining its (convex) faces.
            # Assumes simplexes: Every node except one defines every face.
            assert (num_nodes == 4), "We assume simplexes in 3D (i.e. 4 nodes)"
            node_coords = nodes[:, nodes_idx]

            ids = np.arange(num_nodes)
            polyhedron = [node_coords[:, ids != i] for i in np.arange(num_nodes)]
            test_points = cells_ref[:, test_cells_ptr]  # Test only points not inside another polyhedron.
            in_poly = pp.geometry_property_checks.point_in_polyhedron(
                polyhedron=polyhedron, test_points=test_points, tol=tol
            )

        else:
            logger.warning(f"A grid of dimension {g.dim} encountered. Skip!")
            continue

        # Update pointer to which cell centers to use as test points
        in_poly_ids = test_cells_ptr[in_poly]  # id of cells inside this polyhedron
        test_cells_ptr = test_cells_ptr[~in_poly]  # Keep only cells not inside this polyhedron

        # Update mapping
        indices = np.append(indices, in_poly_ids)
        indptr = np.append(indptr, indptr[-1] + in_poly_ids.size)

    data = np.ones(indices.size)

    coarse_fine = sps.csc_matrix((data, indices, indptr))

    assert (indices.size == g_ref.num_cells), "Every fine cell should be inside exactly one coarse cell"
    return coarse_fine

