import numpy as np
import scipy.linalg
import porepy as pp
from scipy.spatial import ConvexHull

from GTS.geological_model_visualization import GeologicalModel


def get_shearzone_planes():
    """ Get the intersection points for all shearzones in
        the geological model.
    """
    gm = GeologicalModel()
    pts = gm.export_intersections()
    planes = list(pts.keys())

    shearzones = {k: {} for k in planes}

    for plane in planes:
        fp = FitPlane(pts[plane])
        shearzones[plane]['proj'] = fp.proj

        shearzones[plane]['n'] = fp.n

    return shearzones


class FitPlane:
    """ Fit a 2D plane to a cloud of 3D points.

    Computed quantities:
        plane (np.array 1x3): Equation of plane: a*X + b*Y + c = Z
        n (np.array 1x3): Unit normal vector of plane
        proj (np.array 3xn): Point cloud projected to the plane.
    """

    def __init__(self, p: np.ndarray):
        """
        p (np.ndarray): Array of points (3 x n).
        """
        assert(p.shape[0] == 3)
        self.p = p
        self.Np = p.shape[1]

        self.plane = self.fit_plane()
        n = np.hstack((self.plane[0:2], -1))
        self.n = n / np.linalg.norm(n)  # Normalized normal vector

        # Projections
        self.proj = self.project_p()

    def fit_plane(self):
        """ Fit point cloud to plane

        Plane given by a*X + b*Y + c = Z
        """
        p = self.p  # - np.mean(self.p, axis=1)  # Subtract mean for numerical accuracy
        A = np.c_[(p[0], p[1], np.ones(self.Np))]
        C, *_ = scipy.linalg.lstsq(A, p[2])  # coefficients
        # C = np.linalg.inv(A.T.dot(A)).dot(A.T).dot(p[2]) # <-- Same as above

        return C

    def project_p(self):
        """ Project points to plane

        x_p = x_0 - (n . x_0) n
        """
        p = self.p.T
        l = np.atleast_2d(p.dot(self.n)).T
        return (p - l * self.n).T


def convex_hull(p: np.ndarray) -> np.ndarray:
    """ Construct the convex hull of a planar set of 2D points

    Given a point cloud in a plane, compute its convex hull, and
    return the vertices in ccw order of this hull.

    Parameters:
        p (np.ndarray, 3 x n): Planar point cloud
    Returns:
        np.ndarray, 3 x m: Vertices of convex hull, sorted ccw.

    The procedure is as follows:
    * For a given point set, assert that they lie in a plane (or map them to one).
    (1) Map the points to the plane z = 0
    (2) Construct the convex hull (algorithm guarantees ccw ordering.
    (3) Rotate points back to original plane.

    The following methods from PorePy will be utilized directly:
    project_plane_matrix(...)
    """

    assert(pp.geometry_property_checks.points_are_planar(p))

    # Move to mean
    mean = np.mean(p, axis=1)
    pm = np.subtract(p.T, mean).T  # Use these for calculations, then transform back

    # Rotation matrix
    R = pp.map_geometry.project_plane_matrix(pm, check_planar=False)

    pm_2d = np.dot(R, pm)[:2, :]

    # Compute convex hull
    hull = ConvexHull(pm_2d.T)
    vertices_pl = pm_2d[:, hull.vertices]  # Guaranteed to be ccw
    vertices_pl = np.vstack((vertices_pl, np.zeros(vertices_pl.shape[1])))

    # Rotate back
    vertices = np.dot(R.T, vertices_pl)

    # Translate
    vertices = np.add(vertices.T, mean).T

    return vertices


""" 
The two methods below
    - _fit_normal_to_points
    - plane_from_points
constitute the new (23.11.2019) way of computing plane of best fit.
Public usage through plane_from_points.

"""


def _fit_normal_to_points(points: np.ndarray) -> np.array:
    """ Compute a normal from a collection of points.

    Source: http://www.ilikebigbits.com/2015_03_04_plane_from_points.html

    Parameters:
    points : np.ndarray (3,n)
        Array of points

    Returns
    normal : np.array (3,)
        Normalized normal vector to plane
    """
    assert ((points.shape[0] == 3) & (points.shape[1] >= 3))
    pts = points.copy()
    N = pts.shape[1]

    # Calculate centroid and shift points
    centroid = np.atleast_2d(np.sum(pts, axis=1) / N).T
    pts = pts - centroid
    x, y, z = pts[0], pts[1], pts[2]

    # Compute the dot products, x*x, x*y, x*z, y*y, y*z, z*z
    xx = np.dot(x, x)
    xy = np.dot(x, y)
    xz = np.dot(x, z)
    yy = np.dot(y, y)
    yz = np.dot(y, z)
    zz = np.dot(z, z)

    # Compute the determinants
    det_x = yy * zz - yz * yz
    det_y = xx * zz - xz * xz
    det_z = xx * yy - xy * xy

    det_max = max(det_x, det_y, det_z)
    if det_max <= 0.0:
        return None  # The points do not span a plane

    # Pick path with best conditioning
    if det_x == det_max:
        normal = np.array([det_x, xz * yz - xy * zz, xy * yz - xz * yy])
    elif det_y == det_max:
        normal = np.array([xz * yz - xy * zz, det_y, xy * xz - yz * xx])
    else:  # det_z == det_max
        normal = np.array([xy * yz - xz * yy, xy * xz - yz * xx, det_z])

    return normal / np.sqrt((normal ** 2).sum())  # Normalize


def plane_from_points(points: np.ndarray) -> np.ndarray:
    """ Compute a plane for a given point cloud.

    Compute the fitted plane to a point cloud, returning
    the points projected to the plane

    Parameters:
    points : nd.ndarray (3, n)
        Array of n 3D points

    Returns:
    proj : nd.ndarray (3, n)
        Array of n 3D points, projected to the best fit plane

    """
    assert ((points.shape[0] == 3) & (points.shape[1] >= 3)), "Wrong input shape."

    N = points.shape[1]
    centroid = np.sum(points, axis=1) / N
    pointsR = (points.T - centroid).T
    normal = _fit_normal_to_points(pointsR)
    assert (np.isclose(np.sqrt((normal ** 2).sum()), 1))  # Assert normalized

    proj = centroid + pointsR.T - np.outer(np.dot(pointsR.T, normal), normal)

    # Relative error of projected points
    error = np.sqrt(((proj - points.T) ** 2).sum(axis=1))
    length = np.sqrt(((points.T) ** 2).sum(axis=1))
    rel_error = error / length
    print(f"Sum of pointwise relative errors: {rel_error.sum():.4f}")

    return proj.T
