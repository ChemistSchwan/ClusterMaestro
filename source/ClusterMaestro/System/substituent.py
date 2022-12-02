"""The Substituent-class handles the structures attached to the Core-class. Together they form a
whole Molecule.
"""

import numpy as np
import numpy.linalg as npl
import math
from .structure import Structure
from ClusterMaestro.util.linalg_lib import fitPlaneLTSQ, NN, dihedralAngle
from numpy import cross, eye, dot
from scipy.linalg import expm, norm
from copy import deepcopy
from ase import neighborlist


class Substituent(Structure):
    """Substituent object. Molecular structure, that is attached to a core-object.

    Inherits from Structure but has some additional parameters.

    Parameters
    ----------

    connection:\* 1x3 array
        Vector specifying the connection from the Core to the Substituent.
        Currentliy just the function is there, not the object-property.
    dihedral: float
        Relative, normalized dihedral angle of the substituent relative to the core.
    """

    def find_normal(self):
        """
        Computes the normal vector of the plane spanned by the Substituent.
        For nonplanar Substituents this may still be used to define its orientation.
        """

        data = np.array([atom.position for atom in self])
        normal, _ = fitPlaneLTSQ(data)
        return normal

    def get_dihedral(self, core):
        """
        Computes the relative dihedral angle of a Substituent. Internally used
        for e.g. substituent-rotation for dimer-creation.

        Parameters
        ----------

        core: Core
            The core of the system.
        """
        positionList = np.zeros((4, 3), dtype="float64")
        positionList[2] = self[0].position
        positionList[3] = self[1].position
        for i in range(len(core)):
            if NN(core[i], self[0]):
                positionList[1] = core[i].position
                positionList[0] = core[i - 1].position
        angle = dihedralAngle(
            positionList[0], positionList[1], positionList[2], positionList[3]
        )
        self.dihedral = angle

    def get_connection(self, core):
        """
        Find the rotation of the substituent. Connection is then a normalized vector from
        core-center-of-mass to the substituent.

        Parameters
        ----------

        core: Core
            The core of the system.
        """
        connection = self[0].position - core.center_of_mass()
        connection /= npl.norm(connection)
        self.connection = connection

    def get_connectedAtom(self, core):
        """
        Find the actual vector of the bonding Atoms.

        Parameters
        ----------

        core: Core
            The core of the system.
        """
        vector = np.array([10, 10, 10])
        for i in self:
            for j in core:
                if np.linalg.norm(i.position - j.position) < np.linalg.norm(vector):
                    vector = np.array(i.position - j.position)
                    coreAtom = j
                    subAtom = i
        return vector

    def torsion(self, alpha):
        """
        Rotates the Substituent by a specified angle around the Core-Substituent bond

        Parameters
        ----------

        alpha: float
            angle of rotation
        """

        def rotate_as(vec, axis, angle):
            """
            Rotates the Atoms in the Structure by applying a rotation matrix.

            Parameters
            ----------

            rot_mat: 3x3 array
                rotation matrix
            """
            axis = axis / np.linalg.norm(axis)
            term1 = vec * np.cos(angle)
            term2 = (np.cross(axis, vec)) * np.sin(angle)
            term3 = axis * ((1 - np.cos(angle)) * axis.dot(vec))
            return term1 + term2 + term3

        def rotate_new(vec, axis, angle):
            """
            Rotates vector, but better.
            """
            # print(axis)
            def M(axis, theta):
                return expm(cross(eye(3), axis / norm(axis) * theta))

            M0 = M(axis, angle)
            return dot(M0, vec)

        shift = deepcopy(self[0].position)
        self.positions -= shift

        connection = self.connection
        connection = connection / npl.norm(connection)
        for atom in self:
            atom.position = rotate_new(atom.position, connection, math.radians(alpha))
            atom.position += shift

    def getBondLength(self, atoms=["C", "C"]):
        """
        Returns a list of NN-bondlengths of the substituent-structure as array.

        Parameters
        ----------

        atoms: list
            List of atoms, whose bonds are considered in the analysis.
        """
        cutOff = neighborlist.natural_cutoffs(self)
        neighborList = neighborlist.NeighborList(
            cutOff, self_interaction=False, bothways=True
        )
        neighborList.update(self)
        matrix = neighborList.get_connectivity_matrix()
        bondList = []
        for i, atom1 in enumerate(self):
            for j, atom2 in enumerate(self):
                if matrix[i, j] == 1:
                    at1_sym = atom1.symbol
                    at2_sym = atom2.symbol
                    at1_0_at2_1 = at1_sym == atoms[0] and at2_s == atoms[1]
                    at1_1_at2_0 = at1_sym == atoms[1] and at2_s == atoms[0]
                    if at1_0_at2_1 or at1_0_at2_1:
                        bondList.append(np.linalg.norm(atom1.position - atom2.position))
        return np.array(bondList)

    def getBondAngle(self, atoms=["C", "C", "C"]):
        """
        Returns a list of NN-bondangles in degree of the substitunet structure.

        Parameters
        ----------

        atoms: list
            List of atoms, whose bonds are considered in the analysis. The bondangle
            will be calculated for atoms in that specific order.
        """
        import math

        def unit_vector(vector):
            return vector / np.linalg.norm(vector)

        def angle_between(v1, v2):
            v1_u = unit_vector(v1)
            v2_u = unit_vector(v2)
            return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

        cutOff = neighborlist.natural_cutoffs(self)
        neighborList = neighborlist.NeighborList(
            cutOff, self_interaction=False, bothways=True
        )
        neighborList.update(self)
        matrix = neighborList.get_connectivity_matrix()
        angleList = []
        for i, atom1 in enumerate(self):
            for j, atom2 in enumerate(self):
                for k, atom3 in enumerate(self):
                    if j != i and i != k and i > k:
                        if matrix[i, j] == 1 and matrix[j, k] == 1:
                            at1_s = atom1.symbol
                            at2_s = atom2.symbol
                            at3_s = atom3.symbol
                            bool1 = (
                                at1_s == atoms[0]
                                and at2_s == atoms[1]
                                and at3_s == atoms[2]
                            )
                            bool2 = (
                                at1_s == atoms[2]
                                and at2_s == atoms[1]
                                and at3_s == atoms[0]
                            )
                            if bool1 or bool2:
                                vector1 = atom2.position - atom1.position
                                vector2 = atom2.position - atom3.position
                                angle = angle_between(vector1, vector2)
                                angleList.append(math.degrees(angle))
        return np.array(angleList)
