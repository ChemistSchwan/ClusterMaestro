import numpy as np
import numpy.linalg as npl
from .constants import bond_dist_dict
from scipy.linalg import expm
from math import *
import math


def rot_mat(unit_vec, alpha):
    """Function containing the general rotation matrix in 3D.

    Parameters
    ----------

    unit_vec: 1x3 array
        Unit vector specifying the axis of rotation.
    alpha: float
        Angle of rotation.

    Returns
    -------

    3x3 array
        Rotation matrix.

    """
    alpha = alpha * np.pi / 180.0
    nx = unit_vec[0]
    ny = unit_vec[1]
    nz = unit_vec[2]

    rot_mat = np.array(
        [
            nx * nx * (1 - np.cos(alpha)) + np.cos(alpha),
            nx * ny * (1 - np.cos(alpha)) - nz * np.sin(alpha),
            nx * nz * (1 - np.cos(alpha)) + ny * np.sin(alpha),
            ny * nx * (1 - np.cos(alpha)) + nz * np.sin(alpha),
            ny * ny * (1 - np.cos(alpha)) + np.cos(alpha),
            ny * nz * (1 - np.cos(alpha)) - nx * np.sin(alpha),
            nz * nx * (1 - np.cos(alpha)) - ny * np.sin(alpha),
            nz * ny * (1 - np.cos(alpha)) + nx * np.sin(alpha),
            nz * nz * (1 - np.cos(alpha)) + np.cos(alpha),
        ]
    )
    rot_mat = np.reshape(rot_mat, (3, 3))
    return rot_mat


def parallel_rot_mat(vector1, vector2):
    """Computes the cotation matrix that rotates one vector onto another vector.

    Parameters
    ----------

    vector1: 1x3 array
        Starting vector.
    vector2: 1x3 array
        Target vector.

    Returns
    -------

    3x3 array
        Rotation matrix.
    """
    vector1 = vector1 / npl.norm(vector1)
    vector2 = vector2 / npl.norm(vector2)

    v = np.cross(vector1, vector2)
    c = np.dot(vector1, vector2)
    v_x = np.zeros((3, 3))

    v_x[0][1] = -v[2]
    v_x[0][2] = v[1]
    v_x[1][0] = v[2]
    v_x[1][2] = -v[0]
    v_x[2][0] = -v[1]
    v_x[2][1] = v[0]

    rot = np.identity(3) + v_x + np.dot(v_x, v_x) * 1.0 / (1.0 + c)
    return rot


def fitPlaneLTSQ(XYZ):
    """Computes the fit plane for a set of points in 3D. returns the normal vector of that plane as well as the
    last paramter needed to reconstruct other plane representations.

    Parameters
    ----------

    XYZ: Nx3 array
        Points in 3D space.

    Returns
    -------

    1x3 array
        Normal vector.

    float
        Third parameter for parametric representation of the plane.
    """

    rows, _ = XYZ.shape
    G = np.ones((rows, 3))

    G[:, 0] = XYZ[:, 0]  # X
    G[:, 1] = XYZ[:, 1]  # Y

    Z = XYZ[:, 2]
    # resid, rank, s
    (a, b, c), *_ = npl.lstsq(G, Z, rcond=None)

    normal = np.array([a, b, -1])
    nn = npl.norm(normal)
    normal = normal / nn
    return normal, c


def xyz_to_polar(vector):
    """Computes the polar representation of a vector.

    Parameters
    ----------

    vector: 1x3 array
        Cartesian representation of a vector.

    Returns
    -------

    float
        Radial component.
    float
        Polar component.
    float
        Azimuthal component.
    """
    radial = npl.norm(vector)
    polar = np.arccos(vector[2] / radial)
    azimuthal = np.arctan(vector[1] / vector[0])
    return radial, polar, azimuthal


def polar_to_xyz(radial: float, polar: float, azimuthal: float):
    """Computes the cartesian representation of a vector.

    Parameters
    ----------

    radial: float
        Radial component.
    polar: float
        Polar component.
    azimuthal: float
        Azimuthal component.

    Returns
    -------

    1x3 array
        Cartesian representation of a vector.
    """
    x = radial * np.sin(polar) * np.cos(azimuthal)
    y = radial * np.sin(polar) * np.sin(azimuthal)
    z = radial * np.cos(polar)
    vector = np.array([x, y, z])
    return vector


def linefit_3D(XYZ):
    """Computes a line fit in 3D from the main eigenvector of a single value decomposition.

    Parameters
    ----------
    XYZ: Nx3 array
        Points in 3D space.
    """
    # Returns
    # -------
    #
    #    1x3 array
    #        Direction vector of the fit line.
    #    """
    uu, dd, vv = npl.svd(XYZ)
    return vv[0]


def angle_between_vectors(vector1, vector2, limit=True):
    """Computes the angle between two vectors.

    Parameters
    ----------

    vector1: 1x3 array
        First vector.
    vector2: 1x3 array
        Second vector.
    limit: Bool
        Limits the angle to 180 degrees (computes the inner angle).

    Returns
    -------

    float
        Angle in degrees.
    """
    projection = np.dot(vector1, vector2) / (npl.norm(vector1) * npl.norm(vector2))
    alpha = np.arccos(projection)
    if limit and alpha >= (0.5 * np.pi):
        alpha = np.pi - alpha
    alpha = alpha * 180 / np.pi
    return alpha


def dihedralAngle(vector1, vector2, vector3, vector4):
    """Computes the dihedral angle of the four given vectors.

    Parameters
    ----------

    vector1: 1x3 array
        First vector
    vector2: 1x3 array
        Second Vector
    vector3: 1x3 array
        Third vector
    vector4: 1x3 array
        Fourth vector
    """
    p0 = vector1
    p1 = vector2
    p2 = vector3
    p3 = vector4
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def dist(atom1, atom2):
    return npl.norm(atom1.position - atom2.position)


def dist_vec(atom1, atom2):
    return npl.norm(atom1 - atom2)


def NN(atom1, atom2):
    distance = npl.norm(atom1.position - atom2.position)
    bondDist = max(bond_dist_dict[atom1.symbol], bond_dist_dict[atom2.symbol])
    return bondDist > distance


def rotate(vec, axis, angle):
    """Rodrigues rotation formula.

    Rotates vector <vec> <angle> degrees around axis <axis>

    Parameters
    ----------

    vec: 1x3 array
        Startvector
    axis: 1x3 array
        Rotationaxis
    angle: float
        Angle for rotation (in degrees)
    """
    axis = axis / np.linalg.norm(axis)
    term1 = vec * np.cos(angle)
    term2 = (np.cross(axis, vec)) * np.sin(angle)
    term3 = axis * ((1 - np.cos(angle)) * axis.dot(vec))
    return term1 + term2 + term3


def rotateVector(vector1, vector2, vector3):
    """
    Rotates vector 1 to match the direction of vector2. Vector3 is the starting orientation of the fragment.
    This function is used often in the initilization of systems to rotate substituents or molecules to match a certain direction.

    Parameters
    ----------

    vector1: 1x3 array

    vector2: 1x3 array

    vector3: 1x3 array
    """
    c = np.dot(vector1, vector2) / np.linalg.norm(vector1) / np.linalg.norm(vector2)
    angle = np.arccos(np.clip(c, -1, 1))
    cross = np.cross(vector1, vector2)
    M0 = expm(np.cross(np.eye(3), cross / np.linalg.norm(cross) * angle))
    return np.dot(M0, vector3)


def planeAngle(p0, p1):
    """
    Calculates the angle between two planes.

    Parameters
    ----------

    p0: 1x4 array
        First plane.
    p1: 1x4 array
        Second plane.
    """
    a1, b1, c1, d1 = p0
    a2, b2, c2, d2 = p1
    d = a1 * a2 + b1 * b2 + c1 * c2
    e1 = math.sqrt(a1 * a1 + b1 * b1 + c1 * c1)
    e2 = math.sqrt(a2 * a2 + b2 * b2 + c2 * c2)
    d = d / (e1 * e2)
    A = math.degrees(math.acos(d))
    return A


def planeIntersect(a, b):
    """
    a, b   4-tuples/lists
           Ax + By +Cz + D = 0
           A,B,C,D in order

    output: 2 points on line of intersection, np.arrays, shape (3,)
    """
    a_vec, b_vec = np.array(a[:3]), np.array(b[:3])

    aXb_vec = np.cross(a_vec, b_vec)

    A = np.array([a_vec, b_vec, aXb_vec])
    d = np.array([-a[3], -b[3], 0.0]).reshape(3, 1)

    p_inter = np.linalg.solve(A, d).T

    return p_inter[0] - (p_inter + aXb_vec)[0]


def planeFit(array):
    """
    Returns list. First entry is a 1x4 array, which contains the equation form of a plane:

            Ax + By +Cz + D = 0
            A,B,C,D in order

    Second entry is a fit-value.

    Parameters
    ----------

    array: array
        array with points
    """
    from scipy.optimize import leastsq

    def f_min(X, p):
        plane_xyz = p[0:3]
        distance = (plane_xyz * X).sum(axis=1) + p[3]
        return distance / np.linalg.norm(plane_xyz)

    def residuals(params, signal, X):
        return f_min(X, params)

    p0 = [0, 0, 1, 0]
    sol = leastsq(residuals, p0, args=(None, array))[0]
    return [sol, (f_min(array, sol) ** 2).sum()]


def planeFit_alt(array):
    """
    Returns a 1x4 array with the array. The only difference to planeFit is the representation of the outcome.

    Parameters
    ----------

    array: array
        array with points
    """
    p0 = [0, 0, 1, 0]

    def f_min(X, p):
        plane_xyz = p[0:3]
        distance = (plane_xyz * X).sum(axis=1) + p[3]
        return distance / np.linalg.norm(plane_xyz)

    def residuals(params, signal, X):
        return f_min(X, params)

    from scipy.optimize import leastsq

    sol = leastsq(residuals, p0, args=(None, array))[0]
    return sol


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Parameters
    ----------

    y: 1d array
        Array with possible NaNs.
    """
    return np.isnan(y), lambda z: z.nonzero()[0]


def latticeConvert(lattice):
    """
    Converts a lattice-array (3x3) into two 1x3 arrays, one eith the lengths and one with the angles.

    Parameters
    ----------

    lattice: 3x3 array
    """
    lvec = np.array([np.linalg.norm(lat) for lat in lattice])
    angle = np.array(
        [
            np.degrees(math.acos(np.dot(lat1, lat2) / (lvec[i] * lvec[j])))
            for j, lat2 in enumerate(lattice)
            for i, lat1 in enumerate(lattice)
            if i > j
        ]
    )  # gamma,beta,alpha

    returnAngle = np.array([angle[2], angle[1], angle[0]])

    return lvec, returnAngle


def pairCorrelationFunction_3D(x, y, z, S, rMax, dr):
    """
    Source: *https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation*

    Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S.  This simple
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects.  If no such particles exist, an error is
    returned.  Try a smaller rMax...or write some code to handle edge effects! ;)

    Parameters
    ----------
    x: np.array
        An array of x positions of centers of particles
    y: np.array
        An array of y positions of centers of particles
    z: np.array
        An array of z positions of centers of particles
    S: float
        Length of each side of the cube in space
    rMax: float
        Outer diameter of largest spherical shell
    dr: float
        Increment for increasing radius of spherical shell

    Returns
    -------
    g(r): np.array
        A numpy array containing the correlation function g(r)
    radii: np.array
        A numpy array containing the radii of the spherical shells used to compute g(r)
    reference_indices: np.array
        Indices of reference particles
    """
    from numpy import zeros, sqrt, where, pi, mean, arange, histogram
    from joblib import Parallel, delayed
    import multiprocessing

    coreNumber = int(multiprocessing.cpu_count() / 2)

    bools1 = x > rMax
    bools2 = x < (S - rMax)
    bools3 = y > rMax
    bools4 = y < (S - rMax)
    bools5 = z > rMax
    bools6 = z < (S - rMax)
    (interior_indices,) = where(bools1 * bools2 * bools3 * bools4 * bools5 * bools6)
    num_interior_particles = len(interior_indices)
    if num_interior_particles < 1:
        raise RuntimeError(
            "No particles found for which a sphere of radius rMax\
                will lie entirely within a cube of side length S.  Decrease rMax\
                or increase the size of the cube."
        )
    edges = arange(0.0, rMax + 1.1 * dr, dr)
    num_increments = len(edges) - 1
    #    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = len(x) / S**3

    def firstLoopContent(p):
        index = interior_indices[p]
        d = sqrt((x[index] - x) ** 2 + (y[index] - y) ** 2 + (z[index] - z) ** 2)
        d[index] = 2 * rMax
        (result, bins) = histogram(d, bins=edges, normed=False)
        #        g[p,:] = result / numberDensity
        return result / numberDensity

    g = np.array(
        Parallel(n_jobs=coreNumber)(
            delayed(firstLoopContent)(p) for p in range(num_interior_particles)
        )
    )
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i + 1]) / 2.0
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[:, i]) / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3))
    return (g_average, radii, interior_indices)


def pairCorrelationFunction_NP(x, y, z, S, rMax, dr):
    """
    Source: *https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation*

    Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S.  This simple
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects.  If no such particles exist, an error is
    returned.  Try a smaller rMax...or write some code to handle edge effects! ;)

    Parameters
    ----------
    x: np.array
        An array of x positions of centers of particles
    y: np.array
        An array of y positions of centers of particles
    z: np.array
        An array of z positions of centers of particles
    S: float
        Length of each side of the cube in space
    rMax: float
        Outer diameter of largest spherical shell
    dr: float
        Increment for increasing radius of spherical shell

    Returns
    -------
    g(r): np.array
        A numpy array containing the correlation function g(r)
    radii: np.array
        A numpy array containing the radii of the spherical shells used to compute g(r)
    reference_indices: np.array
        Indices of reference particles
    """
    from numpy import zeros, sqrt, where, pi, mean, arange, histogram
    from joblib import Parallel, delayed
    import multiprocessing

    coreNumber = int(multiprocessing.cpu_count() / 2)

    edges = arange(0.0, rMax + 1.1 * dr, dr)
    num_increments = len(edges) - 1
    radii = zeros(num_increments)

    def firstLoopContent(p):
        index = p
        d = sqrt((x[index] - x) ** 2 + (y[index] - y) ** 2 + (z[index] - z) ** 2)
        d[index] = 2 * rMax
        (result, bins) = histogram(d, bins=edges, normed=False)
        return result

    g = []
    for p in range(len(x)):
        g.append(firstLoopContent(p))
    g = np.array(g)
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i + 1]) / 2.0
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[:, i]) / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3))
    return (g_average, radii)


class SphericalPoint(object):
    def __init__(self, r, theta, phi):
        # radial coordinate, zenith angle, azimuth angle
        self.r = r
        self.theta = theta
        self.phi = phi

    def degrees(self, prnt=False):
        if prnt:
            print(
                "SphericalPoint(%.4f, %.4f deg, %.4f deg)"
                % (self.r, degrees(self.theta) % 360, degrees(self.phi) % 360)
            )
        return [self.r, degrees(self.theta) % 360, degrees(self.phi) % 360]

    def __str__(self):
        return "(%0.4f, %0.4f, %0.4f)" % (self.r, self.theta, self.phi)

    def __repr__(self):
        return "SphericalPoint(%f, %f, %f)" % (self.r, self.theta, self.phi)

    def to_cartesian(self):
        x = self.r * cos(self.phi) * sin(self.theta)
        y = self.r * sin(self.phi) * sin(self.theta)
        z = self.r * cos(self.theta)
        return Point(x, y, z)


class Point(object):
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "(%0.4f, %0.4f, %0.4f)" % (self.x, self.y, self.z)

    def __repr__(self):
        return "Point(%f, %f, %f)" % (self.x, self.y, self.z)

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, f):
        return Point(self.x * f, self.y * f, self.z * f)

    def dist(self, other):
        p = self - other
        return (p.x**2 + p.y**2 + p.z**2) ** 0.5

    def toSpherical(self):
        r = self.dist(Point(0, 0, 0))
        theta = atan2(hypot(self.x, self.y), self.z)
        phi = atan2(self.y, self.x)
        return SphericalPoint(r, theta, phi)
