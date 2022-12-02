"""The Solid-class can handle large structures with multiple molecules and with periodic
boundary conditions. Many functions are similar to the Dimer class, but written for periodic
solids.
"""

from copy import deepcopy
import os, shutil, ase
import numpy as np
from ase import neighborlist
from ase.io import read

import matplotlib.pyplot as plt
from matplotlib import style
import pylab as pl

import ClusterMaestro as cm
import ClusterMaestro.System.molecule
from ClusterMaestro.System.structure import Structure
from ClusterMaestro.util.constants import atomMasses


class Solid(Structure):
    """
    Solid state object.

    Inherits from the Structure object.
    Adds functions, which come in handy for the analysis of larger-scale solid state systems.

    Currently ONLY supports CUBIC CELLS properly. You may try with other cell shapes, but
    proper functionalitiy is not guaranteed.

    Parameters
    ----------

    molecules: List
        List of complete molecules within the solid state. Molecules, which are only partly
        in the unit-cell are beeing completed.
    tags: Array
        Array of actual tags of the atoms.
        Stays consistent during modification of the molecules.
    """

    def __init__(self, molecule, **kwargs):
        super(self.__class__, self).__init__(**kwargs)
        self.arrays = molecule.arrays
        self.set_tags(range(len(self)))
        self.tags = self.arrays["tags"]
        self.molecules = []
        self.cell = molecule.cell
        self._pbc = []

    def __repr__(self):
        return "CMSolid: Atoms: {}, cell: {}".format(len(self), self.cell)

    def _renew_positions_mols(self):
        moleculelength = len(self.molecules[0])
        for i, mol in enumerate(self.molecules):
            for j, atom in enumerate(mol):
                self.positions[i * moleculelength + j] = atom.position

    def getCell(self, fileName):
        """
        Finds the cell of the structure from a turbomole-type file
        (e.g. from a periodic xtb output) and writes it to the cell
        parameter of the Solid-object.

        Parameters
        ----------

        fileName: str
            Name of the file, where the cell parameters are.
        """
        readLattice, latticeLine = False, 0
        lattice = []
        with open(fileName, "r") as source:
            for i, line in enumerate(source):
                if readLattice and i < latticeLine + 4:
                    lattice.append(np.array(line.split(), dtype=np.float64))
                if "$lattice bohr" in line:
                    readLattice = True
                    latticeLine = i
        print(np.array(lattice, dtype=np.float64))
        self.cell = np.array(lattice, dtype=np.float64)

    def expandSystem(self):
        extendSystem = deepcopy(self)
        extendSystem.write("ex_1.xyz", format="xyz")
        counter = 1

        def exAllMols(system, coordinates, numMols):
            addMols = []
            for i in system.molecules[:numMols]:
                i.positions += coordinates
                addMols.append(i)
            for i in addMols:
                system.molecules.append(i)

        if len(self.molecules) > 0:
            mols = True
        else:
            mols = False
        numMols = len(self.molecules)
        for coordinate in range(2):
            addSystem = deepcopy(self)
            addSystem.positions += self.cell[coordinate]
            extendSystem += addSystem
            exAllMols(extendSystem, self.cell[coordinate], numMols)
            addSystem.positions -= 2 * self.cell[coordinate]
            extendSystem += addSystem
            exAllMols(extendSystem, -self.cell[coordinate], numMols)
        addSystem = deepcopy(self)
        addSystem.positions += self.cell[0] + self.cell[1]
        extendSystem += addSystem
        exAllMols(extendSystem, self.cell[0] + self.cell[1], numMols)
        addSystem = deepcopy(self)
        addSystem.positions += self.cell[0] - self.cell[1]
        extendSystem += addSystem
        exAllMols(extendSystem, self.cell[0] - self.cell[1], numMols)
        addSystem = deepcopy(self)
        addSystem.positions += self.cell[1] - self.cell[0]
        extendSystem += addSystem
        exAllMols(extendSystem, self.cell[1] - self.cell[0], numMols)
        addSystem = deepcopy(self)
        addSystem.positions -= self.cell[0] + self.cell[1]
        extendSystem += addSystem
        exAllMols(extendSystem, -1 * self.cell[1] - self.cell[0], numMols)
        addSystem = deepcopy(extendSystem)
        addSystem.positions += self.cell[2]
        extendSystem += addSystem
        exAllMols(extendSystem, self.cell[2], numMols * 9)
        addSystem.positions -= 2 * self.cell[2]
        extendSystem += addSystem
        exAllMols(extendSystem, -2 * self.cell[2], numMols * 9)
        extendSystem.cell *= 3
        return extendSystem

    def analyze(
        self,
        atomsPerMolecule=None,
        writeStructures=False,
        core=None,
        substituent=None,
        poscar=False,
    ):
        """
        Analyses the structure by running the *seperateMolecules*-function.
        Consequently, this function analyses the molecules
        of the given solid state. A short summary of the properties of the
        given system will be printed.

        Parameters
        ----------

        atomsPerMolecule: int or None
            Number of atoms in each molecule. If the atoms are not in the right
            order (first molecule 1, then molecule 2 and so on),
            do NOT set this parameter. The function will default into a slower,
            but more versitile algorithm.
        writeStructure: Bool
            Whether or not the individual molecules are written in files.
        core, substituents:
            see completeMolecules
        poscar: bool
            Does the structure come from a poscar or not? This is important for the
            algorithm, which completes the molecules. *poscar=False* is much faster,
            but can return wrong structures, if the atoms are not in the right order.
        """
        if atomsPerMolecule == None:
            print("Other algorithm")
        else:
            self.seperateMolecules(atomsPerMolecule)
            density = self.calcDensity()
            print("The Density of the system is %.5f g/cm^3\n" % density)
            atomsPerMolList = []
            warningBool = False
            string1, string2 = "%18s" % "Molecule:", "%18s" % "Number of Atoms:"
            for i, molecule in enumerate(self.molecules):
                string1 += "%5d" % i
                string2 += "%5d" % len(molecule)
                if len(molecule) != atomsPerMolecule:
                    warningBool = True
            print(string1)
            print(string2)
            if warningBool:
                print(
                    "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                )
                print(
                    "++                            WARNING                            ++"
                )
                print(
                    "++        Not all molecules have the right number of atoms       ++"
                )
                print(
                    "++   Check for harsh structure distortions/ unwanted reactions   ++"
                )
                print(
                    "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
                )
            if writeStructures:
                if not poscar:
                    completeStructure = self.completeMolecules(
                        core=core, substituent=substituent
                    )
                else:
                    completeStructure = self.completeMolecules(
                        core=core, substituent=substituent, brute=True
                    )
                completeStructure.write("completeStructure.xyz", format="xyz")
                cm.io.write.poscar(
                    "completeStructure.vasp",
                    completeStructure,
                    box=np.array([self.cell[0, 0], self.cell[1, 1], self.cell[2, 2]]),
                )
                for i, molecule in enumerate(self.molecules):
                    molecule.write("mol_%.3d.xyz" % i, format="xyz")

    def analyzeFAST(self, atomsPerMolecule, writeStructures=False):
        """
        Analyses the structure by running the *seperateMolecules*-function.
        Consequently, this function analyses the molecules of the given solid state.
        A short summary of the properties of the given system will be printed.

        Parameters
        ----------

        atomsPerMolecule: int
            Number of atoms in each molecule.
        """
        self = cm.lib.solidFuncs.analyze(self, atomsPerMolecule, writeStructures)

    def seperateMolecules_OLD(self, atomsPerMolecule):
        """
        From a structure-object (or one of its inherit), which was read-in as a
        solid state (PBC), the molecules, which are partly outside the box are
        completed and seperated from each other. They are returned in a list of molecules.

        Only works, if the atoms are in a defined order (e.g. created with the included
        init-function and treated with a code, that does not mess up the order).

        Parameters
        ----------

        system: cm.system
            The read-in structure of the system.
        atomsPerMolecule: int
            How many atoms are there per molecule?
        """
        expandedSystem = self.expandSystem()
        expAtomList = []
        for i, atom in enumerate(expandedSystem):
            at_pos = atom.position
            cell = self.cell
            if (
                at_pos[0] < cell[0, 0] + 10
                and at_pos[1] < cell[1, 1] + 10
                and at_pos[2] < cell[2, 2] + 10
            ):
                expAtomList.append(atom)
        expandedSystem = cm.System.solid.Solid(ase.Atoms(expAtomList))
        cutOff = neighborlist.natural_cutoffs(expandedSystem)
        neighborList = neighborlist.NeighborList(
            cutOff, self_interaction=False, bothways=True
        )
        neighborList.update(expandedSystem)
        matrix = neighborList.get_connectivity_matrix()
        newMolecules = []

        def findClosest(extendSystem, CenterPoint, dontLook):
            distance = 100
            dontLook = np.array(dontLook)
            dontLook.flatten()
            for i in extendSystem:
                if i.tag not in dontLook:
                    dist = np.linalg.norm(centerPoint - i.position)
                    if dist < distance:
                        closest = i
                        distance = dist
            return closest

        centerPoint = np.array(
            [sum(self.cell[0] / 2), sum(self.cell[1] / 2), sum(self.cell[2] / 2)]
        )
        dummySystem = deepcopy(self)
        for i, atom in enumerate(expandedSystem):
            atom.tag = atom.tag % len(self)
        for i in range(0, len(self), atomsPerMolecule):
            closest = cm.util.general.findClosest(
                dummySystem, centerPoint, newMolecules
            )
            newMolecules.append([closest.tag])
            for k in newMolecules[-1]:
                result = np.argwhere(matrix[k] == 1)
                for number in result[:, -1]:
                    if number not in newMolecules[-1]:
                        newMolecules[-1].append(number)
        retMolecules = []
        for i in range(len(newMolecules)):
            retMolecules.append([])
            for j in newMolecules[i]:
                retMolecules[i].append(expandedSystem[j])
        for i in range(len(retMolecules)):
            retMolecules[i] = cm.System.molecule.Molecule(ase.Atoms(retMolecules[i]))
        self.molecules = retMolecules

    def seperateMolecules(self, atomsPerMolecule):
        """
        From a structure-object (or one of its inherit), which was read-in as a solid state
        (PBC), the molecules, which are partly outside the box are completed and seperated
        from each other. They are returned in a list of molecules.

        Only works, if the atoms are in a defined order (e.g. created with the included
        init-function and treated with a code, that does not mess up the order).

        Parameters
        ----------

        system: cm.system
            The read-in structure of the system.
        atomsPerMolecule: int
            How many atoms are there per molecule?
        """
        for i in range(0, len(self), atomsPerMolecule):
            newMol = []
            for j in range(atomsPerMolecule):
                newMol.append(deepcopy(self[i + j]))
            self.molecules.append(cm.System.molecule.Molecule(ase.Atoms(newMol)))

    def completeMolecules_OLD(self, atomsPerMolecule):
        """
        From a structure-object (or one of its inherited classes), which was read-in as a
        solid state (PBC), the molecules, which are partly outside the box are completed.
        This makes the visual evaluation of the system easyer and functions as a
        starting point for further evaluation.

        Only works, if the atoms are in a defined order (e.g. created with the included
        init-function and treated with a code, that does not mess up the order).

        Parameters
        ----------

        system: cm.system
            The read-in structure of the system.
        atomsPerMolecule: int
            How many atoms are there per molecule?
        """
        if len(self.molecules) == 0:
            self.seperateMolecules(atomsPerMolecule)
        system = self.molecules
        solidState = cm.System.solid.Solid(ase.Atoms())
        for i in system:
            solidState += i
        return solidState

    def matchPattern(self, core="Ad", substituent="Ph", threads=4):
        """
        Matches the pattern (what is the core/ what are the substituents of the molecule?)
        according to a given structure. This comes in handy, if the *find_substituents*
        or *find_core*-function fails.

        Parameters
        ----------

        core: str
            Core structure. See :ref:`Structure parts`
        substituents: str
            Substituent structure. See :ref:`Structure parts`
        """
        from joblib import Parallel, delayed

        matchMol = cm.initilize.init.dimer(core=core, substituent=substituent)
        moleculeSize = len(matchMol.molecules[0])
        molecules = deepcopy(self)
        molecules = cm.util.general.chunkAtoms(molecules, moleculeSize)
        for i in range(len(molecules)):
            molecules[i] = ClusterMaestro.System.molecule.Molecule(
                ase.Atoms(molecules[i])
            )
        molecules[1].tags += len(molecules[0])
        self.molecules = molecules
        Parallel(n_jobs=threads)(
            delayed(i.matchPattern)(core=core, substituent=substituent)
            for i in self.molecules
        )
        for i in self.molecules:
            i.matchPattern(core=core, substituent=substituent)

    def calcDensity(self):
        """
        Calculates the density of the system in g /cm^3. Onyl supports cells with only 90 deg angles!
        """
        weight = 0
        for i in self:
            weight += atomMasses[i.symbol]
        volume = np.linalg.det(self.cell)
        weight *= 1.66053906660
        return weight / volume

    def makeDense(self, atomsPerMolecule=None, core=None, substituent=None):
        """
        Lammps often returns the structure in a very un-dense manner. In order to
        extract e.g. Dimers from this structure, the molecules need to be shifted
        into the box, which is done with this function.

        Parameters
        ----------

        atomsPerMolecule: int
            see Solid.analyze
        core: Core
            see Solid.analyze
        substituent: Substituent
            see Solid.analyze
        """

        def pointInBox(point, box):
            returnVal = []
            for i, val in enumerate(point):
                if val < 0 or val > box[i, i]:
                    returnVal.append(i)
            if len(returnVal) == 0:
                return [True, None]
            else:
                return [False, returnVal]

        def shiftMolVal(point, box, vals):
            newPoint = np.zeros((3))
            for i in vals:
                if point[i] < box[i, i]:
                    newPoint[i] = +box[i, i]
                elif point[i] > box[i, i]:
                    newPoint[i] = -box[i, i]
            return newPoint

        def shiftMolInBox(molecule, shiftVal):
            for i in molecule:
                i.position += shiftVal
            return molecule

        self.analyze(
            atomsPerMolecule=atomsPerMolecule, core=core, substituent=substituent
        )
        for j, molecule in enumerate(self.molecules):
            inBox = pointInBox(molecule.center_of_mass(), self.cell)
            if not inBox[0]:
                shiftVal = shiftMolVal(molecule.center_of_mass(), self.cell, inBox[1])
                self.molecules[j] = shiftMolInBox(molecule, shiftVal)
        self._renew_positions_mols()

    def getBox(self, dataFile):
        """
        Extracts the box from a lammps-input file.

        Parameters
        ----------

        dataFile: str
            Lammps-inputfile (endung: '.dat')
        """
        with open(dataFile, "r") as source:
            for line in source:
                if "lo" in line:
                    print(line)

    def expandMolecules(self):
        """
        Expands the system according to the molecules in it.
        """
        molSystem = cm.System.solid.Solid(ase.Atoms())
        for molecule in self.molecules:
            molSystem += molecule
            molSystem.molecules.append(molecule)
        molSystem.cell = self.cell
        molSystem = molSystem.expandSystem()
        return molSystem

    def extractDimers(
        self, optimize=False, acc="normal", threads=1, GFN=2, maxCoreDist=12
    ):
        """
        Extracts the directly connected dimers from a solid state. This solid state should
        already be analysed!

        Parameters
        ----------

        optimize: bool
            Whether or not the extracted dimer structure will be optimized (using xTB and their
            default parameters).
        acc: str
            See *ClusterMaestro.util.optimize*.
        threads: int
            See *ClusterMaestro.util.optimize*.
        GFN: int, "ff" or list
            See *ClusterMaestro.util.optimize*. If *GFN* is a List, a sequence will be made, where
            the *GFN* parameters are executed one by one.
        maxCoreDist: float
            Maximum accepted core-core distance, which is still acceptible as a dimer.
        """
        if os.path.exists("dimers"):
            shutil.rmtree("dimers")
        os.mkdir("dimers")
        dimerCounter = 1
        for i, mol1 in enumerate(self.molecules):
            for j, mol2 in enumerate(self.molecules):
                distance = np.linalg.norm(
                    mol1.core.center_of_mass() - mol2.core.center_of_mass()
                )
                if distance < maxCoreDist and i > j:
                    print(f"{dimerCounter} Dimers found and extracted")
                    dimer = deepcopy(mol1)
                    dimer += mol2
                    if optimize:
                        if type(GFN) == list:
                            for metric in GFN:
                                dimer = dimer.optimize(
                                    acc=acc, threads=threads, GFN=metric
                                )
                        else:
                            dimer = dimer.optimize(acc=acc, threads=threads, GFN=GFN)
                    if optimize:
                        dimer.writeEnergy(
                            "dimers/{:05d}.xyz".format(dimerCounter), format="xyz"
                        )
                    else:
                        dimer.write(
                            "dimers/{:05d}.xyz".format(dimerCounter), format="xyz"
                        )
                    dimerCounter += 1
                    if dimerCounter == 27:
                        break

    def pairDistFuncMolecules(
        self, rMax=10, dr=0.1, plotName="pairDistFuncMolecule.pdf", figSize=(4, 3)
    ):
        """
        This function plots the pair distribution function of the given System with respect to
        whole molecules and NOT atoms. It is based on
        *https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation*.

        Returns also the return arrays from the *pairCorrelationFunction_3D*.

        The solid needs to be analysed already!

        Parameters
        ----------

        rMax: float
            Maximum radius to which the PDF will be created.
        plotName: str
            Name of the resulting plot.
        figSize: Tuple
            Dimentions of the resulting figure.
        """
        style.use("seaborn")

        positions = np.empty((len(self.molecules), 3))
        for i, molecule in enumerate(self.molecules):
            positions[i] = molecule.center_of_mass()
            print(molecule)
        dummy = ase.Atoms()
        for i in positions:
            dummy += ase.Atom(symbol="Xe", position=i)
        newSolid = cm.System.solid.Solid(dummy)
        newSolid.cell = self.cell

        while 4 * rMax > newSolid.cell[0, 0]:
            print("exapnding cell for pairDistFuncMolecules ...")
            newSolid = newSolid.expandSystem()

        g_r, r, reference_indices = cm.util.linalg_lib.pairCorrelationFunction_3D(
            newSolid.positions[:, 0],
            newSolid.positions[:, 1],
            newSolid.positions[:, 2],
            newSolid.cell[0, 0],
            rMax,
            dr,
        )

        # Plot defining region:
        plt.figure(figsize=figSize)
        plt.plot(r, g_r)  # , color='black')
        plt.xlabel("Radius / ${\AA}$")
        plt.ylabel("G(r)")
        plt.xlim((0, rMax))
        plt.ylim((0, 1.05 * g_r.max()))
        plt.tight_layout()
        plt.savefig(plotName)

        return g_r, r

    def pairDistFuncSubs(
        self, rMax=10, dr=0.1, plotName="pairDistFuncSubstituents.pdf", figSize=(4, 3)
    ):
        """
        This function plots the pair distribution function of the given System with
        respect to the substituents only and NOT atoms or molecules.
        It is based on
        *https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation*.

        Returns also the return arrays from the *pairCorrelationFunction_3D*.

        The solid needs to be analysed already!

        Parameters
        ----------

        rMax: float
            Maximum radius to which the PDF will be created.
        plotName: str
            Name of the resulting plot.
        figSize: Tuple
            Dimentions of the resulting figure.
        """
        style.use("seaborn")
        positions = []
        for i, molecule in enumerate(self.molecules):
            for j in molecule.substituents:
                positions.append(j.center_of_mass())
        positions = np.array(positions)
        dummy = ase.Atoms()
        for i in positions:
            dummy += ase.Atom(symbol="Xe", position=i)
        newSolid = cm.System.solid.Solid(dummy)
        newSolid.cell = self.cell
        while 4 * rMax > newSolid.cell[0, 0]:
            print("exapnding cell for pairDistFuncMolecules ...")
            newSolid = newSolid.expandSystem()

        g_r, r, reference_indices = cm.util.linalg_lib.pairCorrelationFunction_3D(
            newSolid.positions[:, 0],
            newSolid.positions[:, 1],
            newSolid.positions[:, 2],
            newSolid.cell[0, 0],
            rMax,
            dr,
        )

        # Plot defining region:
        plt.figure(figsize=figSize)
        plt.plot(r, g_r)  # , color='black')
        plt.xlabel("Radius / ${\AA}$")
        plt.ylabel("G(r)")
        plt.xlim((0, rMax))
        plt.ylim((0, 1.05 * g_r.max()))
        plt.tight_layout()
        plt.savefig(plotName)
        return g_r, r

    def pairDistFuncAtoms(
        self,
        rMax=10,
        dr=0.1,
        excludeAtom=None,
        plotName="pairDistFuncAtoms.pdf",
        orgT=False,
        orgE=False,
        figSize=(4, 3),
    ):
        """
        This function plots the pair distribution function of the given System with respect
        to every single atom in the System and NOT the molecules.
        It is based on
        *https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation*.

        Returns also the return arrays from the *pairCorrelationFunction_3D*.

        The solid needs to be analysed already!

        Parameters
        ----------

        rMax: float
            Maximum radius to which the PDF will be created.
        excludeAtom: None or list
            Exclude one or more atom types by giving a list of the atom symbols
            (e.g. ["H"] or ["H", "C"])
        plotName: str
            Name of the resulting plot.
        orgT: Bool
            Calculates the tetrel equlivalent within the organic structure.
        orgE: Bool
            Calculates the chalcogen equivalent within the organic structure.
        figSize: Tuple
            Dimentions of the resulting figure.
        """
        style.use("seaborn")

        newSolid = deepcopy(self)

        while 4 * rMax > newSolid.cell[0, 0]:
            print("expanding cell for pairDistFuncAtoms ...")
            newSolid = newSolid.expandSystem()

        positionArray = []
        if not orgT and not orgE:
            if excludeAtom is not None:
                for atom in newSolid:
                    if atom.symbol not in excludeAtom:
                        positionArray.append(atom.position)
                positionArray = np.array(positionArray)
            else:
                positionArray = newSolid.positions
        elif orgT:
            for molecule in newSolid.molecules:
                neighborList = molecule.core.connectivityMatrix()
                for i, atom in enumerate(molecule.core):
                    if atom.symbol != "H":
                        if np.sum(neighborList[i]) == 3:
                            positionArray.append(deepcopy(atom.position))
        elif orgE:
            for molecule in newSolid.molecules:
                neighborList = molecule.core.connectivityMatrix()
                for i, atom in enumerate(molecule.core):
                    if atom.symbol != "H":
                        if np.sum(neighborList[i]) == 4:
                            positionArray.append(deepcopy(atom.position))

            positionArray = np.array(positionArray)

        g_r, r, reference_indices = cm.util.linalg_lib.pairCorrelationFunction_3D(
            positionArray[:, 0],
            positionArray[:, 1],
            positionArray[:, 2],
            newSolid.cell[0, 0],
            rMax,
            dr,
        )

        # Plot defining region:
        plt.figure(figsize=figSize)
        plt.plot(r, g_r)  # , color='black')
        plt.xlabel("Radius / ${\AA}$")
        plt.ylabel("G(r)")
        plt.xlim((0, rMax))
        plt.ylim((0, 1.05 * g_r.max()))
        plt.tight_layout()
        plt.savefig(plotName)

        return g_r, r

    def getNextNeighbors(self, rMax=12.0, plot=False, plotName="NN_Hist.pdf"):
        """
        Analyses a solid in regards to the next neighbors of the solid. The solid must be analysed.

        Parameters
        ----------

        rMax: float
            Every molecule in a radius of r <= rMax is considered a neighbor. Please choose
            this function wisely by using the *pairDistFuncMolecules* function.
        plot: bool
            If True, a histogram plot will be created.
        plotName: str
            Name of resulting histogram plot.
        """
        NNarray = np.zeros(len(self.molecules))
        for i, mol in enumerate(self.molecules):
            for j, otherMol in enumerate(self.molecules):
                if i != j:
                    if (
                        np.linalg.norm(mol.center_of_mass() - otherMol.center_of_mass())
                        < rMax
                    ):
                        NNarray[i] += 1
        if plot:
            n_bins = 15
            fig, ax = plt.subplots()
            ax.hist(NNarray, bins=n_bins, align="left", rwidth=1, range=(0, 15))
            plt.xlim(0, 15)
            plt.xticks(range(0, 21, 2), range(0, 21, 2))
            plt.xlabel("Number of next neighbors")
            plt.ylabel("Appearances of this number of NN")
            plt.savefig(plotName)
        return NNarray

    def viewNeighborPositions(
        self, reference, core="Ad", substituent="Ph", maxCoreDist=10
    ):
        """
        Uses *rotateDimerOnRef* to fit all two-molecule systems of the *Dimer* and created
        dummy-helium molecules at the positions of their center of mass. If the folder
        *dimers* with extracted dimer structures already exists, this function will use the
        structures from that folder.

        Parameters
        ----------

        reference: Dimer
            Reference-molecule, which will be used for the rotation.
        maxCoreDist: float
            Maximum accepted core-core distance, which is still acceptible as a dimer.
        """
        if not os.path.exists("dimers"):
            self.extractDimers(threads=4, maxCoreDist=maxCoreDist)
        centers = ase.Atoms()
        counter = 0
        for dimerNum in os.listdir("dimers"):
            print(counter)
            dimer = cm.System.Dimer(read(f"dimers/{dimerNum}"))
            dimer.matchPattern(core=core, substituent=substituent)
            rotDimer = dimer.rotateDimerOnRef(reference)
            center = ase.Atom(
                symbol="H", position=rotDimer.molecules[1].center_of_mass()
            )
            centers += center
            counter += 1
        referenceMol = reference.molecules[0]

        return centers, referenceMol

    def getDistributionPlot(
        self, reference, maxCoreDist=10, core="Ad", substituent="Ph"
    ):
        """
        Creates a 2D plot for the distribution of other molecules around a reference molecule.
        If the function *viewNeighborPositions* was executed before and the return objects
        were saved as "centers.xyz" and "referenceMol.xyz", this function will not perform
        this again, but read the two files.

        Parameters
        ----------

        maxCoreDist: float
            Maximum accepted core-core distance, which is still acceptible as a dimer.
        reference: Dimer
            Reference-molecule, which will be used for the rotation.
        """
        plt.clf()
        if os.path.exists("centers.xyz") and os.path.exists("referenceMol.xyz"):
            surroundings = read("centers.xyz")
            ref = cm.System.Molecule(read("referenceMol.xyz"))
        else:
            surroundings, ref = self.viewNeighborPositions(
                reference, core=core, substituent=substituent, maxCoreDist=maxCoreDist
            )
            surroundings.write("centers.xyz")
            ref.write("referenceMol.xyz")
        vals = [[], []]
        distances = []
        for i in surroundings:
            distances.append(np.linalg.norm(i.position))
        distances = np.array(distances)
        for i in surroundings:
            test = cm.util.linalg_lib.Point(i.position[0], i.position[1], i.position[2])

            sph = test.toSpherical()

            sphCoords = sph.degrees()
            vals[0].append(sphCoords[1])
            vals[1].append(sphCoords[2])
        plt.scatter(
            vals[0], vals[1], s=2, c=distances, cmap=plt.get_cmap("cividis"), alpha=0.7
        )
        plt.colorbar(label="core-core distance / $\AA$")

        plt.xlabel(r"polar angle $\theta$ / °")
        plt.ylabel(r"azimut angle $\varphi$ / °")
        plt.xlim(0, 180)
        plt.ylim(0, 360)

        plt.savefig("angleDistrib.pdf")
        plt.savefig("angleDistrib.png", dpi=500)

    def structureFactor(
        self,
        plotName="SF.pdf",
        fast=True,
        arrayFile="sf.npy",
        pbc=True,
        yRange=[-5, 3],
        figSize=(4, 3),
    ):
        """
        Calculates the structure factor. This function may take a while. Saves the
        resulting array in a file named "sf.npy".

        Roughly based on "https://pabloalcain.github.io/physics/structure-factor/".

        Parameters
        ----------

        plotName: str
            Name of the resulting figure.
        fast: bool
            Takes the faster, complied version for the SF calculation.
        pbc: bool
            whether or not periodic boundary conditions are used.
        figSize: Tuple
            Dimentions of the resultung figure.
        """

        def ssf(x, size, q, pbc=False):
            """
            From a series of positions x in a cubic box of length size we get
            the structure factor for momentum q
            """
            #      print("Another Iteration done")
            natoms = np.shape(x)[0]
            sf = 0.0
            for i in range(natoms):
                for j in range(i + 1, natoms):
                    dx = x[j] - x[i]
                    if pbc:
                        for k in range(3):
                            if dx[k] > size / 2:
                                dx[k] -= size
                            if dx[k] < -size / 2:
                                dx[k] += size
                    r = np.linalg.norm(dx)
                    sf += 2 * np.sin(q * r) / (q * r)
            sf /= natoms
            sf += 1
            return sf

        q = np.linspace(0, 11.0, 110)[1:]
        if not os.path.isfile("sf.npy"):
            if fast:
                import ctypes as ct

                _DIRNAME = os.path.dirname(__file__)
                libssf = ct.CDLL(
                    os.path.join(_DIRNAME, "../lib/strukturFaktor/libssf.so")
                )
                ssf_c = libssf.ssf
                ssf_powder_c = libssf.ssf_powder
                ssf_powder_c_PBC = libssf.ssf_powder_periodic

                def structureFactor(x, size, k, rep=2, lebedev=194):
                    natoms = np.shape(x)[0]
                    npoints = len(k)
                    tmp = (ct.c_double * (npoints * 2))()
                    x_p = x.ctypes.data_as(ct.c_void_p)
                    k_p = k.ctypes.data_as(ct.c_void_p)
                    ssf_c.argtypes = [
                        ct.c_void_p,
                        ct.c_int,
                        ct.c_double,
                        ct.c_int,
                        ct.c_int,
                        ct.c_int,
                        ct.c_void_p,
                        ct.c_void_p,
                    ]
                    ssf_c(x_p, natoms, size, npoints, lebedev, rep, k_p, tmp)
                    ssf = np.frombuffer(tmp, dtype=np.double, count=npoints * 2)
                    return ssf.reshape((npoints, 2))

                def structureFactorPowder(x, size, k, periodic):
                    natoms = np.shape(x)[0]
                    npoints = len(k)
                    tmp = (ct.c_double * npoints)()
                    x_p = x.ctypes.data_as(ct.c_void_p)
                    k_p = k.ctypes.data_as(ct.c_void_p)
                    if not periodic:
                        ssf_powder_c.argtypes = [
                            ct.c_void_p,
                            ct.c_int,
                            ct.c_double,
                            ct.c_int,
                            ct.c_void_p,
                            ct.c_void_p,
                        ]
                        ssf_powder_c(x_p, natoms, size, npoints, k_p, tmp)
                        ssf = np.frombuffer(tmp, dtype=np.double, count=npoints)
                    else:
                        ssf_powder_c_PBC.argtypes = [
                            ct.c_void_p,
                            ct.c_int,
                            ct.c_double,
                            ct.c_int,
                            ct.c_void_p,
                            ct.c_void_p,
                        ]
                        ssf_powder_c_PBC(x_p, natoms, size, npoints, k_p, tmp)
                        ssf = np.frombuffer(tmp, dtype=np.double, count=npoints)
                    return ssf

                sf = structureFactorPowder(
                    self.positions, size=self.cell[0, 0], k=q, periodic=pbc
                )
            else:
                size = self.cell[0, 0]
                x = self.positions
                sf = [ssf(x, size, _, pbc) for _ in q]
            sf = np.array(sf)
            np.save(arrayFile, sf)
        style.use("seaborn")
        sf = np.load(arrayFile)
        fig, ax = pl.subplots(figsize=figSize)
        ax.plot(q, sf)
        ax.axhline(y=0, c="k", ls="--")
        ax.set_xlabel(r"Q / $\AA^{-1}$")
        ax.set_ylabel("Structure factor S(Q)")
        ax.legend()
        plt.ylim(yRange[0], yRange[1])
        fig.tight_layout()
        fig.savefig(plotName)
        pl.close()

    def getBondLengthDist(self, part="core", atoms=["H", "H"]):
        """
        Returns a distribution of the bondlengths within the core-structure of the system.
        The system needs to be analyzed.

        Parameters
        ----------

        part: str
            Currently: "core" or "substituent". This will only take one of the structure parts
            into account.
        atoms: list
            Only bonds between these two atoms types are considered.
        """
        averageList, maxList, minList, totalList = [], [], [], []
        for i in self.molecules:
            if part == "core":
                bondList = i.core.getCoreBondlength(atoms=atoms)
            elif part == "substituent":
                for j in i.substituents:
                    bondList = j.getBondLength(atoms=atoms)
            averageList.append(bondList.sum() / len(bondList))
            maxList.append(bondList.max())
            minList.append(bondList.min())
            for j in bondList:
                totalList.append(j)
        totalList = np.array(totalList)
        return totalList

    def getBondAngles(self, part="core", atoms=["C", "C", "C"]):
        """
        Returns a distribution of the bondangles within the core-structure of the system.
        The system needs to be analyzed.

        Parameters
        ----------

        part: str
            Currently: "core" or "substituent". This will only take one of the structure parts
            into account.
        atoms: list
            Only bonds between these three atom types in this order are considered.
        """
        averageList, maxList, minList, totalList = [], [], [], []
        for i in self.molecules:
            if part == "core":
                bondList = i.core.getBondAngle(atoms=atoms)
            elif part == "substituent":
                for j in i.substituents:
                    bondList = j.getBondAngle(atoms=atoms)
            averageList.append(bondList.sum() / len(bondList))
            for j in bondList:
                totalList.append(j)
        totalList = np.array(totalList)
        return totalList
