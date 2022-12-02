"""The Molecule-Class handles complete structures with core- and substituents.
"""

import ase, rmsd
from ase.io import read
import numpy as np
import numpy.linalg as npl
from scipy import stats, sparse
from copy import deepcopy

import ClusterMaestro as cm
from ClusterMaestro.util.constants import get_scaling, inv_ptable
import ClusterMaestro.util.linalg_lib as ll
from .core import Core
from .structure import Structure
from .substituent import Substituent


class Molecule(Structure):
    """
    Molecule object.

    Inherits from the Structure object.
    Adds functions to the Structure object, that are helpful to analyse
    molecular cluster materials.

    Parameters
    ----------

    tags: Array
        Array of actual tags of the atoms.
        Stays consistent during modification of the molecules.
    core: Structure
        Isolated core-structure of the molecule without substituents.
    substituents: List of Substituents
        Isolated substituent-structures of the molecule without core-atoms.
        Substitions are organized in a list of Molecule-objects.
    """

    def __init__(self, molecule, **kwargs):
        super(self.__class__, self).__init__(**kwargs)
        self.arrays = molecule.arrays
        self.set_tags(range(len(self)))
        self.tags = self.arrays["tags"]
        self.substituents = []
        self.core = None

    def analyzeOld(self, core_cutoff=25, elements=[]):
        """
        Analyses the structure of a monomer and finds core and substituents
        according to some gaussian-functions.

        Parameters
        ----------
        core_cutoff: int
            N densest Atoms to keep for the iterative cleaning process
            (see _core_atoms for details)
        elements: list
            List with the elements, the core is made out of. Only works,
            if the substituents contain NO element from that list.
            (e.g. organic substituents, inorganic core)
        """
        self.find_core(core_cutoff, elements)
        self.find_substituents()

    def analyzeCluster(self, core="Ad", substituent="Ph"):
        """
        Analyses the structure of a monomer and finds core and substituents according to their
        connectivity matrix as analyzed by the *AgglomerativeClustering* algorithm as implemented
        in sklearn. Works pretty well for large substituents.

        Parameters
        ----------

        core: str
            Core structure. See :ref:`Structure parts`
        substituent: str
            Substituent structure. See :ref:`Structure parts`
        """
        from sklearn.cluster import AgglomerativeClustering

        connectivity = self.connectivity_matrix()

        refMono = cm.initilize.init.monomer(core=core, substituent=substituent)
        n_clusters = 1 + len(refMono.substituents)
        cluster = AgglomerativeClustering(
            n_clusters=n_clusters, connectivity=connectivity
        )
        clustering = cluster.fit(connectivity)
        labels = clustering.labels_
        parts, subs = [], []
        for i in range(n_clusters):
            parts.append([])
        for i in range(len(labels)):
            parts[labels[i]].append(self[i])

        for i in parts:
            if len(i) == len(refMono.core):
                self.core = cm.System.Core(i)
            else:
                subs.append(cm.System.Substituent(i))
        self.substituents = subs

    def analyze(self, core, substituent, core_cutoff=25, elements=[]):
        """
        Analyses the structure of a monomer and finds core and substituents by
        using the *analyzeCluster* and the *analyzeOld* method and checking their validity.

        Core and substituent structures MUST be given for comparison. Currently only
        works for molecules with uniform substituents.

        Parameters
        ----------

        core: str
            Core structure. See :ref:`Structure parts`
        substituent: str
            Substituent structure. See :ref:`Structure parts`
        core_cutoff: int
            N densest Atoms to keep for the iterative cleaning process
            (see _core_atoms for details)
        elements: list
            List with the elements, the core is made out of. Only works, if
            the substituents contain NO element from that list. (e.g. organic
            substituents, inorganic core)
        """
        if core is None or substituent is None:
            raise Exception(
                "Core and/or substituent structure are not known! Please define structure parts for a proper analysis"
            )
        match = False
        reference = cm.initilize.init.monomer(core=core, substituent=substituent)
        if len(elements) == 0:
            coreElems, subElems = (
                reference.core.elements(),
                reference.substituents[0].elements(),
            )
            if len(set(coreElems).intersection(set(subElems))) == 0:
                elements = coreElems
        if len(elements) == 0:
            self.analyzeOld(core_cutoff, elements)
            if len(self.core) == len(reference.core):
                for i in range(len(self.substituents)):
                    if len(self.substituents[i]) != len(reference.substituents[i]):
                        break
                    if i == len(reference.substituents) - 1:
                        match = True
        else:
            newCore = ase.Atoms()
            conMolecule = ase.Atoms()
            for i in self:
                if i.symbol in elements:
                    newCore.append(deepcopy(i))
                else:
                    conMolecule.append(deepcopy(i))
            self.core = cm.System.Core(newCore)
            conMolecule = cm.System.Structure(conMolecule)
            connectivity = conMolecule.connectivity_matrix()
            n_components, component_list = sparse.csgraph.connected_components(
                connectivity
            )
            for i in range(n_components):
                self.substituents.append([])
            for i, atom in enumerate(conMolecule):
                self.substituents[component_list[i]].append(deepcopy(atom))
            for i in range(len(self.substituents)):
                self.substituents[i] = Substituent(self.substituents[i])
            if len(self.core) == len(reference.core):
                for i in range(len(self.substituents)):
                    if len(self.substituents[i]) != len(reference.substituents[i]):
                        break
                    if i == len(reference.substituents) - 1:
                        match = True
        if not match:
            self.analyzeCluster(core=core, substituent=substituent)
            if len(self.core) == len(reference.core):
                for i in range(len(self.substituents)):
                    if len(self.substituents[i]) != len(reference.substituents[i]):
                        break
                    if i == len(reference.substituents) - 1:
                        match = True
        if not match:
            print("THIS DID NOT WORK!!! THE MOLECULE HAS NOT BEEN PROPERLY ANALYZED!!!")
        return match

    def _core_atoms(self, core_cutoff):
        """
        Finds all the Atoms of the Core.

        Starts by separating the hydrogen Atoms from the rest. Then a kernel density
        estimation is perfomed and all but the N densest Atoms are discarded. Partial
        Substituents are then cleaned up in an iterative process
        (see _clean_core for details).

        Parameters
        ----------

        core_cutoff: int
            N densest Atoms to keep for the iterative cleaning process

        Returns
        -------

        Core
            The Core of the molecule without hydrogen atoms
        list
            The hydrogen Atoms of the Molecule
        """
        stripped = Core([atom for atom in self if atom.symbol != "H"])
        hydrogens = [atom for atom in self if atom.symbol == "H"]

        if core_cutoff < len(stripped):
            xyz = stripped.XYZ()
            values = xyz.T

            weights = [inv_ptable[atom.symbol] for atom in stripped]
            kde = stats.gaussian_kde(values, weights=weights)
            density = kde(values)

            idxs = np.argpartition(density, -core_cutoff)[-core_cutoff:]
            stripped = Core([stripped[idx] for idx in idxs])

        core_atoms = self._clean_core(stripped)

        if len(core_atoms) == 0:
            raise Exception(
                "ERROR: No core atoms left! Increase the core_cutoff parameter"
            )
            # This could be automatized to automatically increase the cutoff until it's non zero
        return core_atoms, hydrogens

    def _cut_dangling_atoms(self, core):
        """
        Part of the Core-finding process.

        Removes all Atoms with <= 1 bonding partner as to shave of the substituents.

        Parameters
        ----------

        core: Core
            Core for cleanup
        """

        connectivity = core.connectivity_matrix()

        core_atoms = []
        connection_count = np.sum(connectivity, axis=1)

        for ii, con in enumerate(connection_count):
            if con > 1:
                core_atoms.append(core[ii])
        return Core(core_atoms)

    def _clean_core(self, core):
        """
        Part of the Core-finding process.
        Iterative process where each iteration all the Atoms with <= 1 bonding partners
        are removed from the Core.

        Parameters
        ----------

        core: Core
            Core for cleanup

        Returns
        -------

        Core
            The Core of the molecule"""
        l0 = len(core)
        core_atoms = self._cut_dangling_atoms(core)
        l1 = len(core_atoms)

        max_iterations = 10
        ii = 0
        while l0 != l1:
            if ii >= max_iterations:
                raise Exception(
                    "Cleaning up cores did not converge. Resulting cores will have unwanted extra atoms"
                )
            l0 = l1
            core_atoms = self._cut_dangling_atoms(core_atoms)
            l1 = len(core_atoms)
            ii = ii + 1
        return core_atoms

    def find_core(self, core_cutoff=25, elements=[]):
        """
        Finds all the Core Atoms (see _core_atoms for details) and reattaches the hydrogen
        Atoms to adamantane Cores.

        Parameters
        ----------

        core_cutoff: int
            N densest Atoms to keep for the iterative cleaning process
            (see _core_atoms for details)
        elements: list
            List with the elements, the core is made out of. Only works, if the substituents
            contain NO element from that list. (e.g. organic substituents, inorganic core)
        """
        if len(elements) != 0:
            core = ase.Atoms()
            for i in self:
                if i.symbol in elements:
                    core += i
            core = cm.System.core.Core(core)
        else:
            core, hydrogens = self._core_atoms(core_cutoff)
            # Reattach hydrogen atoms
            if core[0].symbol == "C":
                for atom in core:
                    for hydrogen in hydrogens:
                        dist = npl.norm(atom.position - hydrogen.position)
                        if dist <= 1.3:
                            core = core + [hydrogen]
        self.core = core

    def aromatic_carbons(self) -> "list of Atoms":
        """
        Identifies all the aromatic carbon Atoms in the Dimer based on the amount
        of bonding partners each carbon Atom has.
        Should be changed to use connectivity matrix.
        """

        connectivity = self.connectivity_matrix()
        con_count = np.sum(connectivity, axis=1)

        arom_cs = []
        for ii, con in enumerate(con_count):
            if (con == 3 or con == 2) and self[ii].symbol == "C":
                arom_cs.append(self[ii])
        return arom_cs

    def find_substituents(self) -> "list of Substituents":
        """
        Identifies all the Substituents in the Molecule.
        Core must be known.
        """
        substituents = []
        if self.core is None:
            raise Exception("No Core-structure defined!")
        coreList, connectList = [], []
        for i in self.core:
            coreList.append(i.tag)
        for i, atom in enumerate(self):
            if atom.tag not in coreList:  # and atom.symbol == "C":
                for j in self.core:
                    rad_i = cm.util.constants.bond_dist_dict[atom.symbol]
                    rad_j = cm.util.constants.bond_dist_dict[j.symbol]
                    epsilon = (rad_i + rad_j) * 0.55
                    if np.linalg.norm(j.position - atom.position) < epsilon:
                        connectList.append(atom)
        usedList = []
        for i in coreList:
            usedList.append(i)
        for i, conAtom in enumerate(connectList):
            substituents.append([])
            substituents[i].append(conAtom)
            usedList.append(conAtom.tag)
            newAdded = True
            while newAdded:
                newAdded = False
                for j, substAtom in enumerate(substituents[i]):
                    for k, atom in enumerate(self):
                        if ll.NN(substAtom, atom) and atom.tag not in usedList:
                            substituents[i].append(atom)
                            usedList.append(atom.tag)
                            newAdded = True
            substituents[i] = Substituent(substituents[i])
            self.substituents.append(substituents[i])
        return substituents

    def XYZ(self):
        """
        Atomic positions of the Structure.

        Returns
        -------

        Nx3 array
            atom positions
        """
        return np.reshape(np.array([atom.position for atom in self]), (-1, 3))

    def replace_core(self, new_atom_spec, core_cutoff=25):
        """
        Assumes the Molecule's Core is adamantane and replaces it with a different one.
        Tertiary carbons are replaced with the metal and secondary ones with sulfur.
        First all the Substituents are translated away from the Molecule's center of mass.
        Then the core atoms are translated outward as well and replaced by the new Atom type.
        Lastly the Molecule's Atoms are updated.

        Should be generalized to support back and forth conversions between different Cores.

        .. Warning::
           This function is outdated. It is recommended to start a new structure from the
           *initialize* package. However, to some users this function might still be useful,
           so it will be left in.


        Parameters
        ----------

        new_atom_spec: str
            Name of the new main Core Atom. Currently supports Si, Ge and Sn
        """
        scale_rings, scale_S, scale_new_atom = get_scaling(new_atom_spec)

        ring_atoms = []
        for ring in self.substituents:
            ring = deepcopy(ring)
            core_center = self.core.center_of_mass()
            new_pos = (ring.center_of_mass() - core_center) * scale_rings + core_center

            trans_vec = new_pos - ring.center_of_mass()
            ring.translate(trans_vec)
            ring_atoms.append(ring)

        group_atoms = []
        for ring in ring_atoms:  # self.substituents:
            ring_atoms = [atom for atom in ring]
            group_atoms += ring_atoms

        new_core_list = []
        core = Core([atom for atom in self.core if atom.symbol != "H"])
        core_CoM = core.center_of_mass()

        for atom in core:
            CoM_dist = npl.norm(atom.position - core_CoM)
            new_atom = deepcopy(atom)
            if CoM_dist <= 1.71:
                new_atom.position = (
                    atom.position - core_CoM
                ) * scale_new_atom + core_CoM
                new_atom.symbol = new_atom_spec
            else:
                new_atom.position = (atom.position - core_CoM) * scale_S + core_CoM
                new_atom.symbol = "S"
            new_core_list.append(new_atom)

        new_core_list += group_atoms
        new_core_list = sorted(new_core_list, key=lambda a: a.tag)
        newMolecule = Molecule(ase.Atoms(new_core_list))
        newMolecule.analyzeOld(core_cutoff=core_cutoff)
        return newMolecule

    def matchPattern(self, core="Ad", substituent="Ph"):
        """
        Matches the pattern (what is the core/ what are the substituents of the molecule?)
        according to a given structure. This comes in handy, if the *find_substituents*
        or *find_core*-function fails and the molecule was created with ClusterMaestro.

        Parameters
        ----------

        core: str
            Core structure. See :ref:`Structure parts`
        substituents: str
            Substituent structure. See :ref:`Structure parts`
        """
        matchMol = cm.initilize.init.monomer(core=core, substituent=substituent)
        coresize, subsize = len(matchMol.core), len(matchMol.substituents[0])
        core, substituents = [], [[], [], [], []]
        for i, atom in enumerate(self):
            if i < coresize:
                core.append(atom)
            elif i < coresize + subsize:
                substituents[0].append(atom)
            elif i < coresize + 2 * subsize:
                substituents[1].append(atom)
            elif i < coresize + 3 * subsize:
                substituents[2].append(atom)
            elif i < coresize + 4 * subsize:
                substituents[3].append(atom)
        for i in substituents:
            i = Substituent(i)
        for i in range(len(substituents)):
            substituents[i] = cm.System.substituent.Substituent(
                ase.Atoms(substituents[i])
            )
        self.core = cm.System.core.Core(ase.Atoms(core))
        self.substituents = substituents

    def rotateMonomerOnRef(self, reference, mono):
        """
        Rotates a monomer in a way, that the core of the first molecule matches the core
        of a reference system. It is recommended, that this reference system is created
        on-the-fly using the *init* package of the ClusterMaestro module.

        Parameters
        ----------
        reference: Molecule
            Reference-molecule, which will be used for the rotation.
        dimer: Molecule
            Molecule, which will be rotated.
        """
        molecule_1 = deepcopy(reference.positions)
        rotMono = deepcopy(mono)
        molecule_2 = deepcopy(mono.positions)
        weights = np.zeros(len(molecule_1))
        for atom in mono.core:
            weights[atom.tag] = 1
        rotMono.positions = rmsd.kabsch_fit(molecule_2, molecule_1, W=weights)
        for i in rotMono:
            for j in rotMono.core:
                if i.tag == j.tag:
                    j.position = i.position
        for i in rotMono:
            for j in rotMono.substituents:
                for k in j:
                    if i.tag == k.tag:
                        k.position = i.position
        return rotMono
