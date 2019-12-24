import numpy as np
import pandas as pd
from pathlib import Path

import GTS as gts
import porepy as pp


# def convex_plane(shearzone, coords='gts', path=None):
#     """ Calculate shear-zone planes from shear-zones.
#
#     Given a reference to a shear-zone, calculate the best-fit plane,
#     and return convex vertices of the planes.
#
#     Parameters:
#     shearzone (str or list): Names of shearzones to construct planes of (S1_1, S1_2, ...)
#     coords (str, optional): Coordinate system to use. Default: 'gts'.
#     path (pathlib.Path): Path/to/01BasicInputData/
#
#     Returns:
#     dict: the pointcloud projected to a convex, fitted plane.
#
#     """
#     isc = gts.ISCData(path=path)
#
#     if isinstance(shearzone, str):
#         shearzone = [shearzone]
#     elif shearzone is None:
#         shearzone = ['S1_1', 'S1_2', 'S1_3', 'S3_1', 'S3_2']
#
#     shearzones = {sz: {} for sz in shearzone}
#
#     for sz in shearzones:
#         sz_cloud = isc.get_shearzone(sz=sz, coords=coords)
#         fp = gts.FitPlane(sz_cloud)
#
#         shearzones[sz]['proj'] = fp.proj
#         shearzones[sz]['n'] = fp.n
#         shearzones[sz]['vertices'] = gts.convex_hull(fp.proj)
#
#     return shearzones
#
#
# def fracture_network(shearzone, export: bool = False, path=None, **network_kwargs):
#     """ Make a fracture network from a selection of shear-zones.
#
#     Parameters:
#         shearzones (str or list): Shearzones to make fracture network of.
#             if 'None': use all known shearzones.
#         export (bool): Export network to vtk.
#         network_kwargs: export kwargs
#             domain (dict): domain boundaries, dict of 'xmin', 'xmax', ...
#         path (pathlib.Path): Path/to/01BasicInputData/
#     """
#     if isinstance(shearzone, str):
#         shearzone = [shearzone]
#     elif shearzone is None:
#         shearzone = ['S1_1', 'S1_2', 'S1_3', 'S3_1', 'S3_2']
#
#     convex = convex_plane(shearzone, coords='gts', path=path)
#
#     fractures = [pp.Fracture(convex[sz]['vertices']) for sz in shearzone]
#     network = pp.FractureNetwork3d(fractures)
#
#     domain = network_kwargs.get('domain', None)
#     if domain is not None:
#         network.impose_external_boundary(domain=domain)
#
#     name = network_kwargs.get('name', None)
#     if export:
#         if name is None:
#             name = 'NEW_fracture_network.vtu'
#         if name[-4:] != '.vtu':
#             name = name + '.vtu'
#         network.to_vtk(name)
#
#     return network


def convex_plane(shearzone_names, coord_system='gts', path=None):
    """ Compute vertices for the convex polygon of the projected point cloud
    to the plane of best fit for each shear-zone is shearzone_names.

    Data imported from gts.ISCData()

    Parameters:
    shearzone_names : str or list
        Names of shear-zones to construct convex planes of
        Input values, e.g.: 'S1_1', 'S1_2', 'S1_2', or a list of them.
    coord_system : str, Default: 'gts'
        Name of coordinate system to use
        Input values: Either 'gts' or 'swiss'
    path : pathlib.Path
        Path/to/01BasicInputData/

    Returns:
    convex_shearzones : pd.DataFrame
        Convex polygon of projected points to best fit plane

    """
    isc = gts.ISCData(path=path)

    if isinstance(shearzone_names, str):
        shearzone_names = [shearzone_names]
    elif shearzone_names is None:
        shearzone_names = ['S1_1', 'S1_2', 'S1_3', 'S3_1', 'S3_2']
    assert (isinstance(shearzone_names, list))

    results = []
    for sz in shearzone_names:
        point_cloud = isc.get_shearzone(sz=sz, coords=coord_system)
        proj = gts.plane_from_points(point_cloud)  # projection

        convex_vertices = gts.convex_hull(proj)

        frame = pd.DataFrame(data=convex_vertices.T, columns=('x_proj', 'y_proj', 'z_proj'))
        frame['shearzone'] = sz
        results.append(frame)
    df = pd.concat(results, ignore_index=True)
    return df


def fracture_network(shearzone_names, export: bool = False, path=None, **network_kwargs):
    """ Make a fracture network from a selection of shear-zones.

    Parameters:
        shearzone_names : str or list
            Shearzones to make fracture network of.
            if 'None': use all known shearzones.
        export : bool
            Export network to vtk.
        path : pathlib.Path or str
            Path/to/01BasicInputData/
            Alternatively, provide 'windows' or 'linux' default paths for each system.
        network_kwargs : kwargs
            domain : dict
                keys 'xmin', 'xmax', etc. of domain boundaries.

    """
    if isinstance(path, str):
        if path == 'linux':
            _root = Path.cwd()  # should be path/to/mastersproject/src/mastersproject
        elif path == 'windows':
            _root = Path('C:/Users/Haakon/OneDrive/Dokumenter/FORSKNING/mastersproject/src/mastersproject')
        else:
            raise ValueError("Invalid input path type")
        path = _root / 'GTS/01BasicInputData'

    if isinstance(shearzone_names, str):
        shearzone_names = [shearzone_names]
    elif shearzone_names is None:
        shearzone_names = ['S1_1', 'S1_2', 'S1_3', 'S3_1', 'S3_2']
    assert (isinstance(shearzone_names, list))

    convex = convex_plane(shearzone_names, coord_system='gts', path=path)
    fractures = [pp.Fracture(
        convex.loc[convex.shearzone == sz,
                   ('x_proj', 'y_proj', 'z_proj')].to_numpy().T
    ) for sz in shearzone_names]

    network = pp.FractureNetwork3d(fractures)

    # domain_default = {'xmin': -6, 'xmax': 80, 'ymin': 55, 'ymax': 150, 'zmin': 0, 'zmax': 50}
    domain = network_kwargs.get('domain', None)
    if domain is not None:
        network.impose_external_boundary(domain=domain)

    name = network_kwargs.get('name', None)
    if export:
        if name is None:
            name = 'RIGHT_fracture_network.vtu'
        if name[-4:] != '.vtu':
            name = name + '.vtu'
        network.to_vtk(name)

    return network
