"""The Dimer-class handles all structures with more than one molecule without periodic
boundary conditions. Also structures with three or more molecules can be handled
"""

import os, shutil, ase
import numpy as np
import numpy.linalg as npl
from scipy import sparse
from scipy.spatial import ConvexHull
import rmsd
import matplotlib.pyplot as plt
from matplotlib import style

from copy import deepcopy
from .structure import Structure
from .molecule import Molecule

import ClusterMaestro as cm
import ClusterMaestro.util.linalg_lib as ll


class Oligomer(Structure):
    """
    Dimer object.

    Inherits from the Structure object.
    Adds functions to the ASE-Atom object, that are helpful to
    analyse and handle molecular cluster materials.

    Although this class is named "Dimer", it is able to deal with any number
    of molecules in an agglomerate and is frequently used by the developer
    to deal with large particles.

    Parameters
    ----------

    molecules: List
        List of molecules within the dimer-object.
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

    def _mixed_analysis(
        self, algorithm, core1, substituent1, core2, substituent2, core_cutoff, elements
    ):
        """
        This function handles the case, that a dimer consisting of different monomers is
        being analysed.

        Be careful with the elements list.
        """
        if algorithm == "old":
            except_string = "The old algorithm does not work with different molecules"
            except_string += "in the system. \nPlease choose a different algorithm"
            raise ValueError(except_string)
        ref_mol1 = cm.initilize.init.monomer(core=core1, substituent=substituent1)
        ref_mol2 = cm.initilize.init.monomer(core=core2, substituent=substituent2)

        for molecule in self.molecules:
            if len(molecule) == len(ref_mol1):
                if algorithm == "try":
                    molecule.analyze(
                        core=core1,
                        substituent=substituent1,
                        core_cutoff=core_cutoff,
                        elements=elements,
                    )
                elif algorithm == "new":
                    molecule.analyze(core=core1, substituent=substituent1)
            if len(molecule) == len(ref_mol2):
                if algorithm == "try":
                    molecule.analyze(
                        core=core2,
                        substituent=substituent2,
                        core_cutoff=core_cutoff,
                        elements=elements,
                    )
                elif algorithm == "new":
                    molecule.analyze(core=core2, substituent=substituent2)

    def analyze(
        self,
        algorithm="try",
        core="Ad",
        substituent="Ph",
        core2=None,
        substituent2=None,
        core_cutoff=25,
        elements=[],
        justMols=False,
    ):
        """
        Identifies the individual molecules and their respective Substituents and Core.

        Parameters
        ----------

        algorithm: str
            Algorithm to find the core- and substituent structures.
            *old*: using old method, as implemented by MS, based on *find_core* (see Molecule).
            *new*: using new method, as implemented by SS, based on *analyze* (see Molecule)
                (usually more reliable).
            *try*: using first the new one and then, if this fails uses the old one.
        core: str
            Core-structure. Will be used by the new algorithm and testing for accuracy.
            See :ref:`Structure parts`
        substituent: str
            Substituent-structure.  Will be used for testing of accuracy.
            See :ref:`Structure parts`
        core2: str or None
            Core-structure. Will be used by the new algorithm and testing for accuracy.
            Only define this variable, if a mixed dimer is present.
            See :ref:`Structure parts`
        substituent2: str or None
            Substituent-structure.  Will be used for testing of accuracy.
            Only define this variable, if a mixed dimer is present.
            See :ref:`Structure parts`
        core_cutoff: int
            The N densest atoms to keep after the initial KDE. See find_cores for more info.
        elements: list
            List with the elements, the core is made out of. Only works, if the substituents
            contain NO element from that list. (e.g. organic substituents, inorganic core)
        justMols: bool
            Whether or not the analyze function will only identify the molecules and not the parts.
            This might be enough depending on the analysis
        """
        connectivity = self.connectivity_matrix()
        n_components, component_list = sparse.csgraph.connected_components(connectivity)
        for molIdx in range(n_components):
            molecule = Molecule(
                ase.Atoms(
                    [
                        self[i]
                        for i in range(len(component_list))
                        if component_list[i] == molIdx
                    ]
                )
            )
            self.molecules.append(molecule)
        if justMols:
            return None
        if core2 != None and substituent2 != None:
            self._mixed_analysis(
                algorithm, core, substituent, core2, substituent2, core_cutoff, elements
            )
        else:
            for molecule in self.molecules:
                if algorithm == "try":
                    molecule.analyze(
                        core=core,
                        substituent=substituent,
                        core_cutoff=core_cutoff,
                        elements=elements,
                    )
                elif algorithm == "old":
                    molecule.analyzeOld(core_cutoff, elements)
                elif algorithm == "new":
                    molecule.analyze(core=core, substituent=substituent)

    def rm_outer_ph(self):
        """
        Computes the distances of the N Substituents from the center of mass of the Dimer.
        Returns a list of the N-2 closest Substituents as the outer two are assumed to
        not participate in the intermolecular bond.
        """
        CoM = self.center_of_mass()
        groups = []
        for molecule in self.molecules:
            groups = groups + [
                substituent.copy() for substituent in molecule.substituents
            ]

        N_groups = len(groups)
        ring_CoM_list = np.zeros((N_groups, 3))
        idxs = np.zeros(N_groups)

        for ii, ring in enumerate(groups):
            idxs[ii] = ii
            ring = cm.System.Structure(ring)
            ring_CoM_list[ii] = ring.center_of_mass()

        CoM_dists = np.array([npl.norm(ring_CoM - CoM) for ring_CoM in ring_CoM_list])

        ring_to_idx = dict(zip(CoM_dists, idxs))
        CoM_dists = np.sort(CoM_dists)
        CoM_dists = CoM_dists[:-2]

        inner_idxs = [ring_to_idx[i] for i in CoM_dists]
        return [groups[int(i)] for i in inner_idxs]

    def get_core_distance(self, limit=None):
        """
        Computes the core-core distance between the CoMs of the Molecules in the Dimer.

        Can also compute all the core-core distances in a system with more than two molecules.

        Parameters
        ----------

        limit: none or float
            If (and only if) the system contains more than 2 molecules, a limit can be given,
            which acts as an upper bound of the core-core distances.

        Returns
        -------

        float or list
            If a dimer is given, a float is returnt, while a list of floats is returned, if
            the system contains more than two molecules.
        """
        if len(self.molecules) < 2:
            text = "Could not get core-distance for this system. Less than 2 molecules were found."
            raise Exception(text)
        elif len(self.molecules) == 2:
            CoM1 = self.molecules[0].core.center_of_mass()
            CoM2 = self.molecules[1].core.center_of_mass()
            distance = npl.norm(CoM1 - CoM2)
            return distance
        else:
            distances = []
            CoMs = []
            for i, mol1 in enumerate(self.molecules):
                CoMs.append(mol1.core.center_of_mass())
            for i in range(len(self.molecules)):
                for j in range(len(self.molecules)):
                    if not i >= j:
                        distances.append(npl.norm(CoMs[i] - CoMs[j]))
            return distances

    def set_core_distance(self, distance):
        """
        Sets the core-distance to a given value, starting from a given structure.

        Parameters
        ----------

        distance: float
            Target-core-core distance
        """
        CoM1 = self.molecules[0].core.center_of_mass()
        CoM2 = self.molecules[1].core.center_of_mass()
        currentDist = npl.norm(CoM1 - CoM2)
        vector = (CoM1 - CoM2) / currentDist
        vector *= distance - currentDist
        self.molecules[1].positions -= vector
        for i in range(len(self)):
            if i < len(self.molecules[0]):
                self.molecules[0][i].tag = i
            else:
                self.molecules[1][i - len(self.molecules[0])].tag = i
        for i in range(len(self)):
            for j in range(len(self.molecules[1])):
                if self[i].tag == self.molecules[1][j].tag:
                    self[i].position = self.molecules[1][j].position
        self._renew_positions("all")

    def mean_pairwise_angle(self):
        """
        Computes the mean pariwise angle of the Dimer.

        The outer most substituents are thrown out. Then, for each Substituent the
        closest Substituent in the other Molecule is found. The angles between their
        normal vectors is computed and averaged.

        Returns
        -------

        float
            Mean pairwise angle.
        """
        groups = self.rm_outer_ph()

        N_groups = len(groups)
        CoMs = np.zeros((N_groups, 3))
        normal_vecs = np.zeros((N_groups, 3))

        for ii, ring in enumerate(groups):
            ring = cm.System.Substituent(ring)
            CoMs[ii] = ring.center_of_mass()
            normal_vecs[ii] = ring.find_normal()

        dists = np.zeros((N_groups, N_groups))
        for ii, com_i in enumerate(CoMs):
            for jj, com_j in enumerate(CoMs):
                dists[ii, jj] = npl.norm(com_i - com_j)
                if ii == jj:
                    dists[ii, jj] = 100.0

        closest_rings = np.zeros(N_groups)
        for ii, dist in enumerate(dists):
            closest_rings[ii] = np.argmin(dist)

        pw_angles = np.zeros(N_groups)
        for ii, ev in enumerate(closest_rings):
            groups[ii] = cm.System.Substituent(groups[ii])
            groups[int(ev)] = cm.System.Substituent(groups[int(ev)])
            normal1 = groups[ii].find_normal()
            normal2 = groups[int(ev)].find_normal()
            alpha = ll.angle_between_vectors(normal1, normal2)
            pw_angles[ii] = alpha
        avg_angle = np.sum(pw_angles) / N_groups
        return avg_angle

    def bonding_substituents(self, epsilon=5.5):
        """
        Old function, not in active support.
        """

        mol1, mol2 = self.molecules

        bonding_phs = 0
        for Ph1 in mol1.substituents:
            Ph1 = cm.System.Substituent(Ph1)
            for Ph2 in mol2.substituents:
                Ph2 = cm.System.Substituent(Ph2)
                if npl.norm(Ph1.center_of_mass() - Ph2.center_of_mass()) <= 5.5:
                    bonding_phs += 1
        return bonding_phs

    def get_core_tags(self):
        tags = []
        for i in self.molecules:
            for j in i.core:
                if j.symbol != "H":
                    tags.append(j.tag)
        return tags

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
            molecules[i] = Molecule(ase.Atoms(molecules[i]))
        molecules[1].tags += len(molecules[0])
        self.molecules = molecules
        Parallel(n_jobs=threads)(
            delayed(i.matchPattern)(core=core, substituent=substituent)
            for i in self.molecules
        )
        for i in self.molecules:
            i.matchPattern(core=core, substituent=substituent)

    def extractDimers(
        self, optimize=False, acc="normal", threads=1, GFN=2, maxCoreDist=12
    ):
        """
        Extracts the directly connected dimers from a Dimer. This Dimer should
        already be analysed!

        Parameters
        ----------

        optimize: bool
            Whether or not the extracted dimer structure will be optimized
            (using xTB and their default parameters).
        acc: str
            See *ClusterMaestro.util.optimize*.
        threads: int
            See *ClusterMaestro.util.optimize*.
        GFN: int, "ff" or list
            See *ClusterMaestro.util.optimize*. If *GFN* is a list, a sequence will be made,
            where the *GFN* parameters are executed one by one.
        maxCoreDist: float
            Maximum accepted core-core distance, which is still acceptible as a dimer.
        """
        print("This is a dimer, not a solid")
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
                    print("{} Dimers found and extracted".format(dimerCounter))
                    dimer = deepcopy(mol1)
                    dimer += mol2
                    if optimize:
                        if type(GFN) == list:
                            for metric in GFN:
                                print("Optimizing")
                                dimer = dimer.optimize(
                                    acc=acc, threads=threads, GFN=metric
                                )
                        else:
                            dimer = dimer.optimize(acc=acc, threads=threads, GFN=GFN)
                    if optimize:
                        dimer.writeEnergy(
                            "dimers/{:05d}.xyz".format(dimerCounter)
                        )  # , format="xyz")
                    else:
                        dimer.write(
                            "dimers/{:05d}.xyz".format(dimerCounter), format="xyz"
                        )
                    dimerCounter += 1
                    if dimerCounter == 27:
                        break

    def analyzeParticle(self, returnSurface=False):
        """
        This function creates a sphere inside the given partice, excluding outer molecules.
        From this sphere, the relevant molecules can be extracted. The very outer molecules
        have no physical meaning for an amorphous material, as they are surface molecules
        and not bulk molecules.

        The system should be analyzed already, and should be made from a large number of
        molecules (> 100), otherwise most molecules are considered surface molecules.

        Parameters
        ----------

        returnSuface: bool
            Returns both, the core system as well as the removed suface moleculesi, if *True*.

        Returns
        -------
        : ClusterMaestro.System.dimer.Dimer
            New system with only relevant atoms and its density as an object parameter (.density).
        """
        center = self.center_of_mass()
        centerList, distanceList = [], []
        dummySystem = ase.Atoms()
        for mol in self.molecules:
            dummySystem += ase.Atom(symbol="Xe", position=mol.center_of_mass())
        hull = ConvexHull(dummySystem.positions)

        outerDistList = []
        for atom in dummySystem:
            if atom.index in hull.vertices:
                atom.symbol = "Au"
                outerDistList.append(np.linalg.norm(atom.position - center))

        radius = min(outerDistList)
        volume = (4 / 3) * np.pi * radius**3
        massSum = 0

        for atom in self:
            if np.linalg.norm(atom.position - center) < radius:
                massSum += cm.util.constants.atomMasses[atom.symbol]
        massSum *= 1.66053906660
        density = massSum / volume
        print(
            "The approximate density of the System is: {:.3f} g / cm^3".format(density)
        )

        coreSystem, outerSystem = ase.Atoms(), ase.Atoms()
        molecules, cutMols = [], []
        for i, mol in enumerate(self.molecules):
            if dummySystem[i].symbol == "Xe":
                coreSystem += mol
                molecules.append(mol)
            else:
                outerSystem += mol
                cutMols.append(mol)

        coreSystem = cm.System.Oligomer(coreSystem)
        coreSystem.molecules = molecules
        coreSystem.density = density

        if returnSurface:
            outerSystem = cm.System.Oligomer(outerSystem)
            outerSystem.molecules = molecules
            outerSystem.density = density
            return coreSystem, outerSystem
        else:
            return coreSystem

    def rotateDimerOnRef(self, reference):
        """
        Rotates a dimer system in a way, that the core of the first molecule matches the core
        of a reference system. It is recommended, that this reference system is created
        on-the-fly using the *init* package of the ClusterMaestro module.

        Parameters
        ----------
        reference: Oligomer
            Reference-molecule, which will be used for the rotation.
        dimer: Oligomer
            Dimer-molecule, which will be rotated.
        """
        molecule_1 = deepcopy(reference.positions)
        rotDimer = deepcopy(self)
        molecule_2 = deepcopy(self.positions)
        weights = np.zeros(len(molecule_1))
        for atom in self.molecules[0].core:
            weights[atom.tag] = 1
        rotDimer.positions = rmsd.kabsch_fit(molecule_2, molecule_1, W=weights)
        for i in rotDimer:
            for j in rotDimer.molecules:
                for k in j:
                    if i.tag == k.tag:
                        k.position = i.position
        return rotDimer

    def structureFactor(
        self,
        plotName="SF.pdf",
        fast=True,
        arrayFile="sf.npy",
        yRange=[-5, 3],
        width=7,
        height=5,
    ):
        """
        Calculates the structure factor. This function may take a while.
        Saves the resulting array in a file named "sf.npy".

        Roughly based on "https://pabloalcain.github.io/physics/structure-factor/".

        Parameters
        ----------

        plotName: str
            Name of the resulting figure.
        fast: bool
            Uses the faster, complied version for the SF calculation.
        """
        import pylab as pl
        import matplotlib.pyplot as plt

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
        if not os.path.isfile(arrayFile):
            if fast:
                import time
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

                def structureFactorPowder(x, size, k):
                    natoms = np.shape(x)[0]
                    npoints = len(k)
                    tmp = (ct.c_double * npoints)()
                    x_p = x.ctypes.data_as(ct.c_void_p)
                    k_p = k.ctypes.data_as(ct.c_void_p)
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
                    return ssf

                sf = structureFactorPowder(self.positions, size=self.cell[0, 0], k=q)
            else:
                size = self.cell[0, 0]
                x = self.positions
                sf = [ssf(x, size, _) for _ in q]
            sf = np.array(sf)
            np.save(arrayFile, sf)
        from matplotlib import style

        style.use("seaborn")
        sf = np.load(arrayFile)
        fig, ax = pl.subplots(figsize=(width, height))
        ax.plot(q, sf)
        ax.axhline(y=0, c="k", ls="--")
        ax.set_xlabel(r"Q / $\AA^{-1}$")
        ax.set_ylabel("Structure factor S(Q)")
        #        ax.legend()
        plt.ylim(yRange[0], yRange[1])
        fig.tight_layout()
        fig.savefig(plotName)
        pl.close()

    def pairDistFuncAtoms(
        self,
        rMax=10,
        dr=0.1,
        excludeAtom=None,
        plotName="pairDistFuncAtoms.pdf",
        orgT=False,
        orgE=False,
    ):
        """
        This function plots the pair distribution function of the given System with
        respect to every single atom in the System and NOT the molecules. It is based
        on *https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation*.

        Returns also the return arrays from the *pairCorrelationFunction_3D*.

        The solid needs to be analysed already!

        Parameters
        ----------

        rMax: float
            Maximum radius to which the PDF will be created.
        excludeAtom: None or list
            Exclude one or more atom types by giving a list of the atom
            symbols (e.g. ["H"] or ["H", "C"])
        plotName: str
            Name of the resulting plot.
        orgT: Bool
            Calculates the tetrel equlivalent within the organic structure.
        orgE: Bool
            Calculates the chalcogen equivalent within the organic structure.
        """
        import matplotlib.pyplot as plt
        from matplotlib import style

        style.use("seaborn")

        newSolid = deepcopy(self)

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
            positionArray = np.array(positionArray)
        elif orgE:
            for molecule in newSolid.molecules:
                neighborList = molecule.core.connectivityMatrix()
                for i, atom in enumerate(molecule.core):
                    if atom.symbol != "H":
                        if np.sum(neighborList[i]) == 4:
                            positionArray.append(deepcopy(atom.position))
            positionArray = np.array(positionArray)
        g_r, r = cm.util.linalg_lib.pairCorrelationFunction_NP(
            positionArray[:, 0], positionArray[:, 1], positionArray[:, 2], 1, rMax, dr
        )

        # Plot defining region:
        plt.figure(figsize=(4, 3))
        plt.plot(r, g_r)  # , color='black')
        plt.xlabel("Radius / ${\AA}$")
        plt.ylabel("G(r)")
        plt.xlim((0, rMax))
        plt.ylim((0, 1.05 * g_r.max()))
        plt.tight_layout()
        plt.savefig(plotName)
        return g_r, r

    def pairDistFuncMolecules(
        self, rMax=10, dr=0.1, plotName="pairDistFuncMolecule.pdf"
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

        """
        style.use("seaborn")
        positions = np.empty((len(self.molecules), 3))
        for i, molecule in enumerate(self.molecules):
            positions[i] = molecule.center_of_mass()
        dummy = ase.Atoms()

        for i in positions:
            dummy += ase.Atom(symbol="Xe", position=i)
        newSolid = cm.System.Dimer(dummy)

        g_r, r = cm.util.linalg_lib.pairCorrelationFunction_NP(
            newSolid.positions[:, 0],
            newSolid.positions[:, 1],
            newSolid.positions[:, 2],
            1,
            rMax,
            dr,
        )

        # Plot defining region:
        plt.figure(figsize=(4, 3))
        plt.plot(r, g_r)  # , color='black')
        plt.xlabel("Radius / ${\AA}$")
        plt.ylabel("G(r)")
        plt.xlim((0, rMax))
        plt.ylim((0, 1.05 * g_r.max()))
        plt.tight_layout()
        plt.savefig(plotName)
        return g_r, r

    def pairDistFuncSubs(
        self, rMax=10, dr=0.1, plotName="pairDistFuncSubstituents.pdf"
    ):
        """
        This function plots the pair distribution function of the given System with
        respect to the substituents only and NOT atoms or molecules. It is based on
        *https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation*.

        Returns also the return arrays from the *pairCorrelationFunction_3D*.

        The solid needs to be analysed already!

        Parameters
        ----------

        rMax: float
            Maximum radius to which the PDF will be created.
        plotName: str
            Name of the resulting plot.
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
        newSolid = cm.System.Dimer(dummy)

        g_r, r = cm.util.linalg_lib.pairCorrelationFunction_NP(
            newSolid.positions[:, 0],
            newSolid.positions[:, 1],
            newSolid.positions[:, 2],
            1,
            rMax,
            dr,
        )

        # Plot defining region:
        plt.figure(figsize=(4, 3))
        plt.plot(r, g_r)  # , color='black')
        plt.xlabel("Radius / ${\AA}$")
        plt.ylabel("G(r)")
        plt.xlim((0, rMax))
        plt.ylim((0, 1.05 * g_r.max()))
        plt.tight_layout()
        plt.savefig(plotName)
        return g_r, r

    def getBondLengthDist(self, part="core", atoms=["H", "H"]):
        """
        Returns a distribution of the bondlengths within the core-structure of the system.
        The system needs to be analyzed.

        Parameters
        ----------

        part: str
            Currently: "core" or "substituent". This will only take one of the structure
            parts into account.
        atoms: list
            Only bonds between these two atoms types are considered.
        """
        averageList, totalList = [], []
        for i in self.molecules:
            if part == "core":
                bondList = i.core.getCoreBondlength(atoms=atoms)
            elif part == "substituent":
                for j in i.substituents:
                    bondList = j.getBondLength(atoms=atoms)
            averageList.append(bondList.sum() / len(bondList))
            for k in bondList:
                totalList.append(k)
        totalList = np.array(totalList)
        return totalList

    def getBondAngles(self, part="core", atoms=["C", "C", "C"]):
        """
        Returns a distribution of the bondangles within the core-structure of the system.
        The system needs to be analyzed.

        Parameters
        ----------

        part: str
            Currently: "core" or "substituent". This will only take one of the structure
            parts into account.
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
            for k in bondList:
                totalList.append(k)
        totalList = np.array(totalList)
        return totalList

    def getBondLengthCoreSub(self):
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
        totalList = []
        for i in self.molecules:
            for j in i.substituents:
                totalList.append(np.linalg.norm(j.get_connectedAtom(i.core)))
        totalList = np.array(totalList)
        return totalList
