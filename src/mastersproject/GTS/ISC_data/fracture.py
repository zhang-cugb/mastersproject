import numpy as np
import pandas as pd

import GTS as gts
import porepy as pp


def convex_plane(shearzone, coords='gts', path=None):
    """ Calculate shear-zone planes from shear-zones.

    Given a reference to a shear-zone, calculate the best-fit plane.

    Parameters:
    shearzone (str or list): Names of shearzones to construct planes of (S1_1, S1_2, ...)
    coords (str, optional): Coordinate system to use. Default: 'gts'.
    path (pathlib.Path): Path/to/01BasicInputData/

    Returns:
    dict: the pointcloud projected to a convex, fitted plane.

    """
    isc = gts.ISCData(path=path)

    if isinstance(shearzone, str):
        shearzone = [shearzone]
    elif shearzone is None:
        shearzone = ['S1_1', 'S1_2', 'S1_3', 'S3_1', 'S3_2']

    shearzones = {sz: {} for sz in shearzone}

    for sz in shearzones:
        sz_cloud = isc.get_shearzone(sz=sz, coords=coords)
        fp = gts.FitPlane(sz_cloud)

        shearzones[sz]['proj'] = fp.proj
        shearzones[sz]['n'] = fp.n
        shearzones[sz]['vertices'] = gts.convex_hull(fp.proj)

    return shearzones


def fracture_network(shearzone, export: bool = False, path=None, **network_kwargs):
    """ Make a fracture network from a selection of shear-zones.

    Parameters:
        shearzones (str or list): Shearzones to make fracture network of.
            if 'None': use all known shearzones.
        export (bool): Export network to vtk.
        network_kwargs: export kwargs
            domain (dict): domain boundaries, dict of 'xmin', 'xmax', ...
        path (pathlib.Path): Path/to/01BasicInputData/
    """
    if isinstance(shearzone, str):
        shearzone = [shearzone]
    elif shearzone is None:
        shearzone = ['S1_1', 'S1_2', 'S1_3', 'S3_1', 'S3_2']

    convex = convex_plane(shearzone, coords='gts', path=path)

    fractures = [pp.Fracture(convex[sz]['vertices']) for sz in shearzone]
    network = pp.FractureNetwork3d(fractures)

    domain = network_kwargs.get('domain', None)
    if domain is not None:
        network.impose_external_boundary(domain=domain)

    name = network_kwargs.get('name', None)
    if export:
        if name is None:
            name = 'NEW_fracture_network.vtu'
        if name[-4:] != '.vtu':
            name = name + '.vtu'
        network.to_vtk(name)

    return network