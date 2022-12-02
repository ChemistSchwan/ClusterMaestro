import ase, inspect, os, sys, shutil
from ase.io import read
import ClusterMaestro, math
from ClusterMaestro.System.molecule import Molecule
from ClusterMaestro.System.substituent import Substituent
from ClusterMaestro.System.core import Core
from ClusterMaestro.System.molecule import Molecule
from ClusterMaestro.System.oligomer import Oligomer
import ClusterMaestro.util.linalg_lib as ll
from copy import deepcopy
import numpy as np


def monomer(molecule="AdPh4", core=None, substituent=None):
    """Initilize Function for Monomer-structures.

    Creates a new monomer of selected chemical formular. Default is AdPh\ :sub:`4`\ .
    To add more molecules, cores or substituents to the selection, use the *expand* submodule.

    All Core- and Substituent structures are pre-optimized by themselves(!) The rotation of
    the substituents and the bondlength between the core and the substituent is NOT optimized.

    Parameters
    ----------

    molecule: str
        Choose from fully pre-optimized molecules:
        AdPh4 (default)
    core: str
        Choose core structures. See :ref:`Structure parts`
    substituent: str or list
        Choose substituents. See :ref:`Structure parts`.

        If a list is given, multiple different substituents can be attached to the core-structure,
        which makes mixed molecules possible.
    """
    path = (
        inspect.getfile(ClusterMaestro).replace("__init__.py", "")
        + "initilize/structures"
    )
    if molecule == "AdPh4" and core == None and substituent == None:
        molecule = Molecule(read("%s/AdPh4.xyz" % path))
        core, substituents = [], [[], [], [], []]
        for i, atom in enumerate(molecule):
            if i < 22:
                core.append(atom)
            elif i < 33:
                substituents[0].append(atom)
            elif i < 44:
                substituents[1].append(atom)
            elif i < 55:
                substituents[2].append(atom)
            elif i < 66:
                substituents[3].append(atom)
        for i in range(len(substituents)):
            substituents[i] = Substituent(ase.Atoms(substituents[i]))
        molecule.core = Core(core)
        molecule.substituents = substituents
        return molecule
    else:
        core = read("%s/cores/%s.xyz" % (path, core), format="xyz")
        coresize = len(core)
        if type(substituent) == list:
            subList = []
            for i in substituent:
                subList.append(read("%s/substituents/%s.xyz" % (path, i), format="xyz"))
        else:
            subList = []
            subList.append(
                read("%s/substituents/%s.xyz" % (path, substituent), format="xyz")
            )
            subList = subList * 8
        subsize = [len(subList[0]), len(subList[1]), len(subList[2]), len(subList[3])]
        molecule = ase.Atoms()
        for atom in core:
            if atom.symbol != "F":
                molecule += atom
        subCounter = 0
        for atom in core:
            if atom.symbol == "F":
                distList = []
                coresize -= 1
                subCounter += 1
                for j in core:
                    if atom.index != j.index:
                        distList.append(np.linalg.norm(atom.position - j.position))
                NN = np.argmin(np.array(distList))
                rotVector = -atom.position + core[NN].position
                addSubstituent = deepcopy(subList[subCounter - 1])
                for i in addSubstituent:
                    i.position = ll.rotateVector(
                        np.array([0, 0, -1]), rotVector, i.position
                    )
                addSubstituent.positions += core[NN].position
                molecule += addSubstituent
        molecule = Molecule(molecule)
        core, substituents = [], []
        for i in range(subCounter):
            substituents.append([])
        sumAtoms = []
        for i, atom in enumerate(molecule):
            if i < coresize:
                core.append(atom)
            else:
                for j in range(len(substituents)):
                    if i < coresize + sum(subsize[: j + 1]) and i >= coresize + sum(
                        subsize[:j]
                    ):
                        substituents[j].append(atom)
        for i in substituents:
            i = Substituent(i)
        for i in range(len(substituents)):
            substituents[i] = Substituent(ase.Atoms(substituents[i]))
        molecule.core = Core(core)
        molecule.substituents = substituents
        return molecule


def dimer(
    molecule="AdPh4",
    core=None,
    substituent=None,
    distance=6.0,
    subAngle=None,
    rotAngle=None,
    xshift=0.0,
    yshift=0.0,
    inverse=True,
):
    """Initilize Function for Monomer-structures.

    Creates a new dimer of selected chemical formular. Default is AdPh4.
    To add more molecules, cores or substituents, go to the installation of the module
    (/ClusterMaestro/System/structures) and add them according to the python-script that is there.

    .. Warning::
       This function might not work as intended for ALL structures that may be included by the
       user. It works well for structures, which are already in the included database.

    All Core- and Substituent structures are pre-optimized by themselves(!).
    The rotation of the substituents is NOT (yet) optimized.

    Parameters
    ----------

    molecule: str
        Choose from fully pre-optimized molecules:
        AdPh4 (default)
    core: str
        Choose core structures from: See :ref:`Structure parts`
    substituent: str or list
        Choose substituents from: See :ref:`Structure parts`
    distance: float
        Distance of the two monomers in \AA
    subAngle: float or None
        Defines the rotation-angle of the substituents in the finished structure.
        Very important, if the dimer will be used for a potential-energy hyperplane-scan.
    rotAngle: float or None
        Defines the relative rotation angle between the two monomers. rotAngle=0 is an
        alternating structure.
    xshift: float
        One dimer gets shifted in x-direction by that amount
    yshift: float
        One dimer gets shifted in y-direction by that amount
    inverse: bool
        Whether or not the structure will be inversed.
    """
    molecule1 = monomer(molecule, core, substituent)
    coreSize, subSize = len(molecule1.core), len(molecule1.substituents[0])
    currentVec = molecule1.core.center_of_mass() - molecule1[1].position
    currentVec /= np.linalg.norm(currentVec)
    if len(molecule1.substituents) != 3:
        for i in range(len(molecule1)):
            molecule1[i].position = ll.rotateVector(
                currentVec, np.array([0, 0, -1]), molecule1[i].position
            )
    molecule1._renew_positions(structure=["substituents", "core"])
    if subAngle != None:
        for i in molecule1.substituents:
            i.get_connection(molecule1.core)
            i.get_dihedral(molecule1.core)
            i.torsion(subAngle - i.dihedral)
            i.get_connection(molecule1.core)
        molecule1._renew_positions(structure=["all"])
    molecule2 = deepcopy(molecule1)
    if subAngle != None:
        for i in molecule2.substituents:
            if len(molecule2.substituents) == 4:
                minusAng = 180 - (2 * (subAngle + 30))
                if not inverse:
                    minusAng -= minusAng
            else:
                minusAng = 180 - (2 * subAngle)
            i.get_connection(molecule2.core)
            i.get_dihedral(molecule2.core)
            i.torsion(subAngle - i.dihedral + minusAng)
            i.get_connection(molecule2.core)
        molecule2._renew_positions(structure=["all"])

    molecule2.tags += len(molecule1)
    for i in molecule2.substituents:
        for j in i:
            j.tag += len(molecule1)
    for i in molecule2.core:
        i.tag += len(molecule1)
    if inverse:
        molecule2.positions *= -1
    else:
        molecule2.positions[:, 2] *= -1
    if rotAngle != None:
        for k in range(len(molecule2.positions)):
            molecule2[k].position = ll.rotate(
                molecule2[k].position, np.array([0, 0, 1]), math.radians(rotAngle)
            )
    for i in range(len(molecule2)):
        molecule2[i].position[0] += xshift
    for i in range(len(molecule2)):
        molecule2[i].position[1] += yshift

    dimer = []
    for i in molecule1:
        i.position[2] += distance / 2
        dimer.append(i)
    for i in molecule2:
        i.position[2] -= distance / 2
        dimer.append(i)
    molecule1._renew_positions(structure=["substituents", "core"])
    molecule2._renew_positions(structure=["substituents", "core"])
    dimer = Oligomer(ase.Atoms(dimer))
    dimer.molecules = [molecule1, molecule2]
    counter = 0
    for i in dimer:
        i.tag = counter
        counter += 1
    return dimer


def mix_dimer(
    molecule1="AdPh4",
    core1=None,
    substituent1=None,
    molecule2="AdPh",
    core2=None,
    substituent2=None,
    distance=6.0,
    subAngle=None,
    rotAngle=None,
    xshift=0.0,
    yshift=0.0,
    inverse=True,
):
    """Initilize Function for Monomer-structures.

    Creates a dimer of selected chemical formular with two different molecules. This is an
    extention of the simple *dimer* function in this very same package.

    .. Warning::
       This function might not work as intended for ALL structures that may be included by
       the user. It works well for structures, which are already in the included database.

    All Core- and Substituent structures are pre-optimized by themselves(!). The rotation of the
    substituents is NOT (yet) optimized.

    Parameters
    ----------

    molecule1: str
        Choose from fully pre-optimized molecules for the first molecule of the dimer.
        AdPh4 (default)
    core1: str
        Choose core structures for the first molecule from: See :ref:`Structure parts`
    substituent1: str or list
        Choose substituents for the first molecule from: See :ref:`Structure parts`
    molecule1: str
        Choose from fully pre-optimized molecules for the second molecule of the dimer:
        AdPh4 (default)
    core1: str
        Choose core structures for the second molecule from: See :ref:`Structure parts`
    substituent1: str or list
        Choose substituents for the second molecule from: See :ref:`Structure parts`
    distance: float
        Distance of the two monomers in \AA
    subAngle: float or None
        Defines the rotation-angle of the substituents in the finished structure.
        Very important, if the dimer will be used for a potential-energy hyperplane-scan.
    rotAngle: float or None
        Defines the relative rotation angle between the two monomers. rotAngle=0 is an
        alternating structure.
    xshift: float
        One dimer gets shifted in x-direction by that amount
    yshift: float
        One dimer gets shifted in y-direction by that amount
    inverse: bool
        Whether or not the structure will be inversed.
    """
    molecule1 = monomer(molecule1, core1, substituent1)
    coreSize, subSize = len(molecule1.core), len(molecule1.substituents[0])
    currentVec = molecule1.core.center_of_mass() - molecule1[1].position
    currentVec /= np.linalg.norm(currentVec)
    if len(molecule1.substituents) != 3:
        for i in range(len(molecule1)):
            molecule1[i].position = ll.rotateVector(
                currentVec, np.array([0, 0, -1]), molecule1[i].position
            )
    molecule1._renew_positions(structure=["substituents", "core"])
    if subAngle != None:
        for i in molecule1.substituents:
            i.get_connection(molecule1.core)
            i.get_dihedral(molecule1.core)
            i.torsion(subAngle - i.dihedral)
            i.get_connection(molecule1.core)
        molecule1._renew_positions(structure=["all"])

    molecule2 = monomer(molecule2, core2, substituent2)
    coreSize, subSize = len(molecule2.core), len(molecule2.substituents[0])
    currentVec = molecule2.core.center_of_mass() - molecule2[1].position
    currentVec /= np.linalg.norm(currentVec)
    if len(molecule2.substituents) != 3:
        for i in range(len(molecule2)):
            molecule2[i].position = ll.rotateVector(
                currentVec, np.array([0, 0, -1]), molecule2[i].position
            )
    molecule2._renew_positions(structure=["substituents", "core"])

    if subAngle != None:
        for i in molecule2.substituents:
            if len(molecule2.substituents) == 4:
                minusAng = 180 - (2 * (subAngle + 30))
                if not inverse:
                    minusAng -= minusAng
            else:
                minusAng = 180 - (2 * subAngle)
            i.get_connection(molecule2.core)
            i.get_dihedral(molecule2.core)
            i.torsion(subAngle - i.dihedral + minusAng)
            i.get_connection(molecule2.core)
        molecule2._renew_positions(structure=["all"])

    molecule2.tags += len(molecule1)
    for i in molecule2.substituents:
        for j in i:
            j.tag += len(molecule1)
    for i in molecule2.core:
        i.tag += len(molecule1)
    if inverse:
        molecule2.positions *= -1
    else:
        molecule2.positions[:, 2] *= -1
    if rotAngle != None:
        for k in range(len(molecule2.positions)):
            molecule2[k].position = ll.rotate(
                molecule2[k].position, np.array([0, 0, 1]), math.radians(rotAngle)
            )
    for i in range(len(molecule2)):
        molecule2[i].position[0] += xshift
    for i in range(len(molecule2)):
        molecule2[i].position[1] += yshift

    dimer = []
    for i in molecule1:
        i.position[2] += distance / 2
        dimer.append(i)
    for i in molecule2:
        i.position[2] -= distance / 2
        dimer.append(i)
    molecule1._renew_positions(structure=["substituents", "core"])
    molecule2._renew_positions(structure=["substituents", "core"])
    dimer = Oligomer(ase.Atoms(dimer))
    dimer.molecules = [molecule1, molecule2]
    counter = 0
    for i in dimer:
        i.tag = counter
        counter += 1
    return dimer


def amorphSolid(
    molecule="AdPh4",
    core=None,
    substituent=None,
    number=1,
    box=[20, 20, 20],
    tolerance=2.0,
    density=None,
    optimize=False,
    useFile=None,
):
    """This function creates an initial stricture of an amorphous solid state of the
    defined molecule. For this task, the software-tool packmol is beeing used. Accordingly,
    this function is basically a wrapper for the packmol software package.

    One must keep in mind, that the actual box, in which the system is generated is
    1.5 Angstrom per edge smaller than the specified box due to the lack of implementation of
    PBC in packmol.

    Also keep in mind, that the packing is a complex problem, which can take a while.

    If you want to write the data into a file, use the ClusterMaestro.util.write package.
    It currently supports vasp-POSCAR and periodic Turbomole-input files. Also, all ASE-utilities
    can be used, as the returned object is an ase.Atoms-object.

    Parameters
    ----------

    molecule: str
        Choose from fully pre-optimized molecules:
        AdPh4 (default)
    core: str
        Choose core structures from: See :ref:`Structure parts`
    substituent: str or list
        Choose substituents from: See :ref:`Structure parts`
    number: int
        Number of molecules located in the box.
    box: 1x3 array
        Edge-length of the CUBIC box in angstrom.
    density: float or None
        Density of the final system in g/cm^3. If this parameter is given, the given box
        will be ignored.
    optimize: bool
        Whether or not the SINGLE molecule will be optimized before packing it in the box.
    useFile: str or None
        If a molecule should be read in and not be created in the process, give the name
        of this molecule here.
    """
    from ClusterMaestro.util.constants import atomMasses

    if os.path.exists(".tmp_packmol"):
        shutil.rmtree(".tmp_packmol")
    os.mkdir(".tmp_packmol")
    os.chdir(".tmp_packmol")
    molecule = monomer(molecule=molecule, core=core, substituent=substituent)
    if optimize:
        molecule = molecule.optimize()
    if useFile != None:
        molecule = Molecule(read("../{}".format(useFile)))
    molecule.write("tmp_molecule.xyz", format="xyz")
    if density != None:
        print(density)
        weight = 0
        for i in molecule:
            weight += atomMasses[i.symbol]
        weight *= 1.66053906660 * number
        volume = weight / density
        boxsite = volume ** (1 / 3)
        box = [boxsite, boxsite, boxsite]
        print("Length of cell: {} \AA".format(boxsite))
    weight = 0
    for i in molecule:
        weight += atomMasses[i.symbol]
    weight *= 1.66053906660 * number
    density = weight / (box[0] * box[1] * box[2])
    print("Density of the periodic system: %.5f g per cm^3" % density)
    print(
        "Box-dimentions of the periodic system: %.2f %.2f %.2f"
        % (box[0], box[1], box[2])
    )
    with open("input.inp", "w") as pmInp:
        pmInp.write("tolerance %.1f\nfiletype xyz\noutput " % tolerance)
        pmInp.write("result.xyz\n\nstructure tmp_molecule.xyz\n")
        pmInp.write(
            "  number %d\n  inside box 0. 0. 0. %.2f %.2f %.2f\nend structure"
            % (number, box[0] - 1.8, box[1] - 1.8, box[2] - 2)
        )
    os.system("packmol < input.inp > packmol.out")
    try:
        solidState = ClusterMaestro.System.solid.Solid(read("result.xyz", format="xyz"))
        solidState.cell = np.array([[box[0], 0, 0], [0, box[1], 0], [0, 0, box[2]]])
        if os.path.exists("result.xyz_FORCED"):
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
            print("++                     WARNING                     ++")
            print("++     Packmol ended without perfect packaging     ++")
            print("++            Doublecheck the structure!!          ++")
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
        os.chdir("..")
        shutil.rmtree(".tmp_packmol")
        return solidState
    except:
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("++                          ! FAIL !                            ++")
        print("++     Packmol was unable to put the molecules in the box!!     ++")
        print("++        Enlargen the box or reduce number of molecules!       ++")
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
