"""Core-Class that handles the core structure of a cluster molecule. Usually, this is usually an
adamantane-type structure and most functions are purpose-written for that.
"""

import numpy as np
import numpy.linalg as npl
from scipy.spatial import ConvexHull
import ase
from ase import neighborlist
from .structure import Structure


class Core(Structure):
    """Core object.

    Inherits from Structure object. This represents the base structure of a molecule, where
    the substituents attach to.

    Parameters
    ----------

    See Structure object.

    """

    def volume(self):
        """Computes the volume of the Core via a convex hull method."""

        positions = np.array([atom.position for atom in self])
        hull = ConvexHull(positions)
        return hull.volume

    def getCoreBondlength(self, atoms=["C", "C"]):
        """
        Returns a list of NN-bondlengths of the core-structure as array, prints the average.

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
                    sym1 = atom1.symbol
                    sym2 = atom2.symbol
                    if (sym1 == atoms[0] and sym2 == atoms[1]) or (
                        sym1 == atoms[1] and sym2 == atoms[0]
                    ):
                        bondList.append(np.linalg.norm(atom1.position - atom2.position))

        return np.array(bondList)

    def getBondAngle(self, atoms=["C", "C", "C"]):
        """
        Returns a list of NN-bondangles in degree of the substitunet structure.

        Parameters
        ----------

        atoms: list
            List of atoms, whose bonds are considered in the analysis. The bondangle will be
            calculated for atoms in that specific order.
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
                            if (
                                atom1.symbol == atoms[0]
                                and atom2.symbol == atoms[1]
                                and atom3.symbol == atoms[2]
                            ) or (
                                atom1.symbol == atoms[2]
                                and atom2.symbol == atoms[1]
                                and atom3.symbol == atoms[0]
                            ):
                                vector1 = atom2.position - atom1.position
                                vector2 = atom2.position - atom3.position
                                angle = angle_between(vector1, vector2)
                                angleList.append(math.degrees(angle))
        return np.array(angleList)
