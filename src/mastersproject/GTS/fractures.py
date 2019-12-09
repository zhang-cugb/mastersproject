import porepy as pp

from GTS.fit_plane import get_shearzone_planes, convex_hull

from GTS.temp_data_extraction import ShearzoneInterception

def get_fractures():
    """ Method to get fractures from the provided ISC data"""

    shearzones = get_shearzone_planes()
    fracs = {}

    for s in list(shearzones.keys()):
        fracs[s] = convex_hull(shearzones[s]['proj'])

    return fracs


def get_fractures_manual():
    """ Method to get fractures from manually loaded ISC shearzone intercept data"""
    planes = ShearzoneInterception().interpolate_shearzones()
    fracs = {}

    for s in list(planes.keys()):
        fracs[s] = convex_hull(planes[s]['proj'])

    return fracs


def export_network(manual=False):
    """ Export fracture network for visualization

    Parameters:
        manual (bool): True if manual fractures are fetched.
    """
    if manual:
        fracs = get_fractures_manual()
        name = 'fracture_network_manual.vtu'
    else:
        fracs = get_fractures()
        name = 'fracture_network.vtu'

    fractures = [pp.Fracture(fracs[key]) for key in list(fracs.keys())]
    network = pp.FractureNetwork3d(fractures)
    network.to_vtk(name)
    return network



