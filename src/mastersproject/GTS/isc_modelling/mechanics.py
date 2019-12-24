import porepy as pp
import numpy as np
from porepy.models.contact_mechanics_model import ContactMechanics

import GTS as gts


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
        network = gts.fracture_network(shearzone_names=None, export=True, path='linux', domain=self.mesh_args)
        self.gb = network.mesh(mesh_args=self.mesh_args)
        self.Nd = self.gb.dim_max()
