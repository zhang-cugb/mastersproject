import numpy as np
from .geological_model_visualization import GeologicalModel

def get_shearzone_planes():
    """ Get the intersection points for all shearzones in
        the geological model.
    """
    gm = GeologicalModel()
    pts = gm.export_intersections()

class FitPlane:

    def __init__(self, p: np.ndarray):
        """
        p (np.ndarray): Array of points (3 x n).
        """
        assert(p.shape[0] == 3)
        self.p = p
        self.Np = p.shape[1]

        self.plane, self.n = self.fit_plane

    def fit_plane(self):
        """ Fit point cloud to plane"""
        p = self.p
        A = np.vstack((p[0:2], np.ones(self.Np))).T
        b = p[3]

        x = np.inv(A.T * A) * (A.T * b)
        n = np.hstack((x[0:2], -1))
        return x, n

    def project_p(self):
        """ Project points to plane"""
