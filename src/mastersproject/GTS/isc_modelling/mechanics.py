import porepy as pp
import numpy as np
from porepy.models.contact_mechanics_model import ContactMechanics
from porepy.models.abstract_model import AbstractModel


import GTS as gts
from GTS.prototype_1.mechanics.isotropic_setup import IsotropicSetup


class ContactMechanicsISC(ContactMechanics):
    """ Implementation of ContactMechanics for ISC"""

    def __init__(self, **kwargs):
        """ Initialize the mechanics class for an ISC geometry."""
        super().__init__()

        # Mesh size arguments
        default_mesh_args = {'mesh_size_frac': 10, 'mesh_size_min': 10}
        self.mesh_args = kwargs.get('mesh_args', default_mesh_args)

        # Bounding box of the domain
        default_box = {'xmin': -6, 'xmax': 80, 'ymin': 55, 'ymax': 150, 'zmin': 0, 'zmax': 50}
        self.box = kwargs.get('box', default_box)

    def create_grid(self):
        """ Create a GridBucket of a 3D domain with fractures
        defined by the ISC dataset.

        The method requires the following attribute:
            mesh_args (dict): Containing the mesh sizes.

        The method assigns the following attributes to self:
            gb (pp.GridBucket): The produced grid bucket.
            box (dict): The bounding box of the domain, defined through minimum and
                maximum values in each dimension.
            Nd (int): The dimension of the matrix, i.e., the highest dimension in the
                grid bucket.

        """
        network = gts.fracture_network(shearzone_names=None, export=True, path='linux', domain=self.box)
        self.gb = network.mesh(mesh_args=self.mesh_args)
        pp.contact_conditions.set_projections(self.gb)
        self.Nd = self.gb.dim_max()


class ContactMechanicsIsotropicISC(IsotropicSetup):

    def __init__(self, **kwargs):
        """ Initialize the mechanics class for an ISC geometry."""
        super().__init__()

        # Mesh size arguments
        default_mesh_args = {'mesh_size_frac': 10, 'mesh_size_min': 10}
        self.mesh_args = kwargs.get('mesh_args', default_mesh_args)

        # Bounding box of the domain
        default_box = {'xmin': -6, 'xmax': 80, 'ymin': 55, 'ymax': 150, 'zmin': 0, 'zmax': 50}
        self.box = kwargs.get('box', default_box)

    def create_grid(self):
        """ Create a GridBucket of a 3D domain with fractures
        defined by the ISC dataset.

        The method requires the following attribute:
            mesh_args (dict): Containing the mesh sizes.

        The method assigns the following attributes to self:
            gb (pp.GridBucket): The produced grid bucket.
            box (dict): The bounding box of the domain, defined through minimum and
                maximum values in each dimension.
            Nd (int): The dimension of the matrix, i.e., the highest dimension in the
                grid bucket.

        """
        network = gts.fracture_network(shearzone_names=None, export=True, path='linux', domain=self.box)
        self.gb = network.mesh(mesh_args=self.mesh_args)
        pp.contact_conditions.set_projections(self.gb)
        self.Nd = self.gb.dim_max()


def run_model(setup: AbstractModel):
    """
    Set up and run a model.

    Parameters
    setup : pp.AbstractModel
        A model that (ultimately) inherits from AbstractModel.

    """
    params = {'folder_name': 'GTS/isc_modelling/results/isotropic_setup_viz'}
    model = setup(params=params)
    model.prepare_simulation()
    model.init_viz()

    tol = 1e-10

    # Get fracture grid(s):
    # Set zero values there to facilitate export.
    for i in [1, 2]:
        g_list = model.gb.grids_of_dimension(i)
        for g in g_list:
            data = model.gb.node_props(g)
            data[pp.STATE]['u_'] = np.zeros((3, g.num_cells))

    # Get the 3D data.
    g3 = model.gb.grids_of_dimension(3)[0]
    d3 = model.gb.node_props(g3)
    # Solve problem
    x = model.assemble_and_solve_linear_system(tol)
    model.after_newton_convergence(x, None, None)  # Distribute solution

    # Get the state, transform it, and save to another state variable
    sol3 = d3[pp.STATE][model.displacement_variable]
    trsol3 = np.reshape(np.copy(sol3), newshape=(g3.dim, g3.num_cells), order='F')
    d3[pp.STATE][model.displacement_variable + '_'] = trsol3

    # Export solution
    model.export_step()

    return model