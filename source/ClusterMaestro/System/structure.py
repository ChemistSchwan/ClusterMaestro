"""The Structure-class represents the base for all other classes. Based on the ASE.Atoms class,
it extends some functionalities that can be used across many length scales. Everything more
specific is handled in the subclasses.
"""

from copy import deepcopy
import os, math, ase
import numpy as np
import numpy.linalg as npl

import matplotlib.pyplot as plt

import ClusterMaestro as cm
from ClusterMaestro.util.constants import bond_dist_dict, plot_dict
from ClusterMaestro.util.linalg_lib import planeFit


class Structure(ase.Atoms):
    """Structure Object.

    Inherits from the ASE-Atoms class.
    See "https://wiki.fysik.dtu.dk/ase/ase/atoms.html" for more information.
    """

    def removeAtom(system, tag):
        """
        Removes atom according to the tag **(not index!)** of the atoms in the system.
        """
        for i, atom in enumerate(system):
            if atom.tag == tag:
                delete = deepcopy(i)
        system.pop(delete)
        return system

    def __add__(self, other: "iterable containing Atoms"):
        molecule = self.copy()
        try:
            for atom in other:
                molecule += atom
        except ValueError:
            raise Exception("Addition is supported only for iterables.")
        return molecule

    def XYZ(self):
        """Atomic positions of the Structure.

        Returns
        -------

        Nx3 array
            atom positions
        """

        return np.reshape(np.array([atom.position for atom in self]), (-1, 3))

    def connectivity_matrix(self):
        """Computes the connectivity matrix of the Structure. Returns an NxN matrix
        (N: number of Atoms) with elements i,j = 1 if atom i and atom j are in bonding
        distance and 0 if not. Does not count bonding of Atom i with itself.

        This connectivity matrix is based on adjusted bondlengths and has proven to be more
        reliable for the tested systems.
        """
        N = len(self)
        connectivity = np.zeros((N, N), dtype=int)
        dr = np.zeros((N, N, 3))
        for ii, atom_i in enumerate(self):
            for jj, atom_j in enumerate(self):
                dr[ii, jj, :] = atom_i.position - atom_j.position

        distance = npl.norm(dr, axis=2)

        for ii, atom_i in enumerate(self):
            for jj, atom_j in enumerate(self):
                rad_i = bond_dist_dict[atom_i.symbol]
                rad_j = bond_dist_dict[atom_j.symbol]
                # epsilon = rad_i if rad_i >= rad_j else rad_j
                epsilon = (rad_i + rad_j) * 0.55
                if atom_i.index != atom_j.index and distance[ii, jj] <= epsilon:
                    connectivity[ii, jj] = 1
        return connectivity

    def connectivityMatrix(self):
        """
        Computes the connectivity matrix using the ase-routine.
        """
        cutOff = ase.neighborlist.natural_cutoffs(self)
        neighborList = ase.neighborlist.NeighborList(
            cutOff, self_interaction=False, bothways=True
        )
        neighborList.update(self)
        matrix = neighborList.get_connectivity_matrix(sparse=False)
        return matrix

    def center_of_mass(self):
        """Computes the center of mass of the Structure

        Returns
        -------

        1x3 array
            center of mass coordinates
        """
        CoM = np.array([0.0, 0.0, 0.0])
        non_H_count = 0
        for atom in self:
            if atom.symbol != "H":
                CoM += atom.position
                non_H_count += 1
        CoM /= non_H_count
        return CoM

    def torsion(self, substituent, alpha):
        """Rotates the Substituent by a specified angle around the Core-Substituent bond

        Parameters
        ----------

        alpha: float
            angle of rotation
        """

        def rotate_as(vec, axis, angle):
            """Rotates the Atoms in the Structure by applying a rotation matrix.

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

        connection = substituent.connection
        for atom in substituent:
            atom.position = rotate_as(atom.position, connection, math.radians(alpha))

        for i in range(len(substituent)):
            self[substituent[i].tag].position = substituent[i].position

    def plot(self):
        """Simple plot-function using a matplotlib 3D-plot.
        Not suitable for publication-ready figures, but sufficient for a quick structure-check.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        for i in self:
            ax.scatter(
                i.position[0],
                i.position[1],
                i.position[2],
                s=plot_dict[i.symbol][0],
                color=plot_dict[i.symbol][1],
            )
        conMat = self.connectivity_matrix()
        for i in range(len(conMat)):
            for j in range(len(conMat[i])):
                if conMat[i, j] == 1:
                    ax.plot(
                        [self[i].position[0], self[j].position[0]],
                        [self[i].position[1], self[j].position[1]],
                        [self[i].position[2], self[j].position[2]],
                        color=plot_dict[self[i].symbol][1],
                    )
        plt.show()

    def _renew_positions(self, structure):
        """This Function coorects the positions of the given structure and fits them to
        the more recent structures. Used internally e.g. for the dimer-creation after a rotation.

        Parameters
        ----------

        structure: List
            Parts of the System, that will be renewed.
            Options: substituents, core, all, molecules
        """
        for structPart in structure:
            if structPart == "substituents":
                for i in self:
                    for j in self.substituents:
                        for k in j:
                            if k.tag == i.tag:
                                k.position = i.position
            elif structPart == "core":
                for i in self:
                    for j in self.core:
                        if i.tag == j.tag:
                            j.position = i.position
            elif structPart == "molecules":
                for i in self:
                    for j in self.molecules:
                        for k in j:
                            if i.tag == j.tag:
                                i.position = j.position
            elif structPart == "all":
                for i in self:
                    for j in self.substituents:
                        for k in j:
                            if k.tag == i.tag:
                                i.position = k.position
                    for j in self.core:
                        if i.tag == j.tag:
                            i.position = j.position

    def appendToXYZ(self, xyzFile):
        """
        Appends the current system to an existing xyz-File to create a visualizable trajectory.
        If the xyz-File does not yet exist, it will be created.

        Parameters
        ----------

        xyzFile: str
            xyz-File, where the coordinated will be added.
        """
        with open(xyzFile, "a") as target:
            target.write("%d\n\n" % len(self))
            for i in self:
                target.write(
                    "%2s%17.8f%17.8f%17.8f%9d\n"
                    % (i.symbol, i.position[0], i.position[1], i.position[2], i.tag)
                )

    def aromatic_carbons(self) -> "list of Atoms":
        """
        Identifies all the aromatic carbon Atoms in the system based on the amount of
        bonding partners each carbon Atom has. Additionally checks, if these atoms are
        in a plane in order to exclude unsaturated carbon atoms.
        """
        connectivity = self.connectivity_matrix()
        aromCs = []
        for i, col in enumerate(connectivity):
            NNList = [self[i].position]
            for j, val in enumerate(col):
                if val != 0 and self[i].symbol in ["C", "Si", "Ge", "Sn"]:
                    NNList.append(self[j].position)
            if len(NNList) in [2, 3, 4]:
                try:
                    plane, certain = planeFit(np.array(NNList))
                    if certain < 0.01:
                        aromCs.append(self[i])
                except:
                    aromCs.append(self[i])
        return aromCs

    def findUnsatAtoms(self):
        """
        Finds all unsaturated carbon atoms with hydrogen. Also pays attention on the
        hybridization of the carbon atom.
        """
        connectivity = self.connectivity_matrix()
        unsatAtoms = []
        aromCs = self.aromatic_carbons()
        aromTags = []
        for i in aromCs:
            aromTags.append(i.index)
        for i, col in enumerate(connectivity):
            NNList = [self[i].position]
            for j, val in enumerate(col):
                if val != 0 and self[i].symbol in ["C", "Si", "Ge", "Sn"]:
                    NNList.append(self[j].position)
            if self[i].index in aromTags and len(NNList) < 4:
                if self[i].symbol in ["C", "Si", "Ge", "Sn"]:
                    unsatAtoms.append(self[i])
            elif self[i].index not in aromTags and len(NNList) < 5:
                if self[i].symbol in ["C", "Si", "Ge", "Sn"]:
                    unsatAtoms.append(self[i])
        return unsatAtoms

    def saturateCH(self, write=False):
        """
        Looks for unsaturated carbon atoms using *findUnsatAtom* and adds hydrogen atoms in
        the approximately right location.

        Parameters
        ----------

        write: bool
            Whether or not an output is written.
        """
        HBond = {"C": 1.07797, "Si": 1.43615, "Ge": 1.52382, "Sn": 1.72072}
        unsatAtoms = self.findUnsatAtoms()
        if write:
            for i in unsatAtoms:
                print(i)
        indexList = [i.index for i in unsatAtoms]
        connectivity = self.connectivityMatrix()
        for i, col in enumerate(connectivity):
            if i in indexList:
                bondingAtomVecs = []
                for j, val in enumerate(col):
                    if val == 1:
                        bondingAtomVecs.append(self[i].position - self[j].position)
                vecH = sum(bondingAtomVecs)
                addH = ase.Atoms("H")
                addH.positions = self[i].position
                vecH = vecH / np.linalg.norm(vecH) * HBond[self[i].symbol]
                addH.positions += vecH
                if write:
                    print(addH)
                self += addH

    def optimize(
        self,
        acc="normal",
        threads=1,
        GFN=2,
        output=False,
        showString=False,
        folder=None,
        inputFile=None,
    ):
        """
        Optimizes the structure using the xTB-package by Grimme et. al. The exact code
        can be found in *ClusterMaestro.util.optimize.opt_xtb*.

        This function **RETURNS** the optimized structure as a new object!

        Parameters
        ----------

        acc: str
            Level of accuracy. Default: normal
            Options: crude, sloppy, loose, lax, normal, tight, vtight, extreme
        threads: int
            Number of threads the optimization is running on. Default: 1
        GFN: int
            GFN-algorithm. Default: 2
        output: Bool
            Whether or not the file "xtbopt.out" is being written and stored.
            Defalt: False
        showString: Bool
            Whether or not the string for the xtbinput is beeing printed.
        folder: str or None
            If not None: keep folder, where the optimization is performed.
            str is the name of that folder
        inputFile: str or None
            If not None: File with detailed input (see xTB Documentation for
                    information on how to set up a proper inputFile).
        """
        newMol = cm.util.optimize.opt_xtb(
            self,
            acc=acc,
            threads=threads,
            GFN=GFN,
            output=output,
            showString=showString,
            folder=folder,
            inputFile=inputFile,
        )
        return newMol

    def writeEnergy(self, filename):
        """
        Writes the structure as in the ase-Function, but adds the energy to the
        xyz-File, if there is one.

        Parameters
        ----------

        filename: str
            Filename of the resulting file.
        """
        self.write("%s_" % filename, format="xyz")
        with open("%s_" % filename, "r") as source:
            with open(filename, "w") as target:
                for number, line in enumerate(source):
                    if number == 1:
                        target.write("%20.10f\n" % self.energy)
                    else:
                        target.write(line)
        os.remove("{}_".format(filename))

    def rotateSystem(self, axis, angle):
        """
        Rotates the whole structure around the given axis by the given angle.

        Parameters
        ----------

        axis: str
            Axis, around which the rotation will happen.
        angle: float
            Angle of the rotation.
        """
        if axis == "x":
            axis = np.array([1, 0, 0])
        elif axis == "y":
            axis = np.array([0, 1, 0])
        elif axis == "z":
            axis = np.array([0, 0, 1])
        angle = math.radians(angle)
        for atom in self:
            atom.position = cm.util.linalg_lib.rotate(atom.position, axis, angle)

    def getBondTypes(self):
        """
        Creates a connectivity matrix and extractes the bond types.

        For a phenyl ring, this would be: [["C", "C"], ["C", "H"]]

        and for [PhSi]4S6, this is: [["Si", "S"], ["Si", "C"], ["C", "C"], ["C", "H"]]
        """
        connectivity = self.connectivity_matrix()
        allConnections = []
        for i, atom1 in enumerate(self):
            for j, atom2 in enumerate(self):
                if connectivity[i, j] == 1:
                    sym1, sym2 = atom1.symbol, atom2.symbol
                    if [sym1, sym2] not in allConnections and [
                        sym2,
                        sym1,
                    ] not in allConnections:
                        allConnections.append([sym1, sym2])
        return allConnections

    def getAngleTypes(self):
        """
        Creates a connectivity matrix and extractes the angle types.

        For a phenyl ring, this would be: [["C", "C", "C"], ["C", "C", "H"]]

        and for [PhSi]4S6, this is: [['S', 'Si', 'S'], ['Si', 'S', 'Si'], ['C', 'Si', 'S'],
                ['C', 'C', 'C'], ['C', 'C', 'Si'], ['H', 'C', 'C']]
        """
        connectivity = self.connectivity_matrix()
        allAngles = []
        for i, atom1 in enumerate(self):
            for j, atom2 in enumerate(self):
                for k, atom3 in enumerate(self):
                    if j != i and i != k and i > k:
                        if connectivity[i, j] == 1 and connectivity[j, k] == 1:
                            s1, s2, s3 = atom1.symbol, atom2.symbol, atom3.symbol
                            if [s1, s2, s3] not in allAngles and [
                                s3,
                                s2,
                                s1,
                            ] not in allAngles:
                                allAngles.append([s1, s2, s3])
        return allAngles

    def elements(self):
        """
        Finds all unique elements in a given structure object and exports them as list.
        """
        elements = []
        for i in self:
            if i.symbol not in elements:
                elements.append(i.symbol)
        return elements
