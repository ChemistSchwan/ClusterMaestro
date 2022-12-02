import ClusterMaestro as cm
import numpy as np
from copy import deepcopy
from ase import neighborlist
from ase.io import read
import ase, os, shutil, time
import matplotlib.pyplot as plt

# from MDAnalysis.transformations.wrap import unwrap
# import MDAnalysis as mda


def cleanTrajectory_OLD(trajFile, core, substituent, box):
    """
    *Currently way too slow*

    Takes a trajectory file and calculates the RMSD. also writes the completed molecules into a new trajectory file
    (for a better visualization).

    Parameters
    ----------

    trajFile: str
        File with many xyz-like coordinates of the evolving structure.
    core: str
        Core structure, See :ref:`Structure parts`.
    substituent: str
        Substituent structure, See :ref:`Structure parts`
    box: 1x3 array
        Dimentions of the cubic box.
    """
    if os.path.exists(".trajSeperate"):
        shutil.rmtree(".trajSeperate")
    os.mkdir(".trajSeperate")
    os.chdir(".trajSeperate")
    shutil.copyfile("../%s" % trajFile, trajFile)
    files = []
    firstLine = ""
    xyzString = ""
    with open(trajFile, "r") as source:
        for i, line in enumerate(source):
            if i == 0:
                firstLine = line
                xyzString = firstLine
            elif line == firstLine:
                files.append(xyzString)
                xyzString = firstLine
            else:
                xyzString += line
    files.append(xyzString)
    for i, string in enumerate(files):
        with open("%.4d.xyz" % i, "w") as target:
            target.write(string)
    counter = 0
    for i, string in enumerate(files):
        start = time.time()
        system = cm.System.solid.Solid(read("%.4d.xyz" % i, format="xyz"))
        system.cell = box
        completeSystem = system.completeMolecules(
            core=core, substituent=substituent, box=box
        )
        completeSystem.appendToXYZ("xyzTrajectory.xyz")
        counter += 1
        print("Structures done: %d" % counter)
        print("time needed: %.3f s" % (time.time() - start))
    shutil.copyfile("xyzTrajectory.xyz", "../xyzTrajectory.xyz")
    os.chdir("..")
    shutil.rmtree(".trajSeperate")


# def cleanTrajectory(fileName, trajFile, core, substituent, box):
#    """
#    *Massive speedup due to the usage of MDAnalysis (https://userguide.mdanalysis.org/index.html)*
#
#    Takes a trajectory file and calculates the RMSD. also writes the completed molecules into a new trajectory file
#    (for a better visualization).
#
#    Parameters
#    ----------
#
#    fileName: str
#        Name of the printed file.
#    trajFile: str
#        File with many xyz-like coordinates of the evolving structure.
#    core: str
#        Corestructure, as in init.
#    substituent: str
#        Substituentstructure, as in init
#    box: 1x3 array
#        Dimentions of the cubic box.
#    """
#    u = mda.Universe(trajFile)
#    bonds = cm.util.general.getBonds(len(u.atoms), core=core, substituent=substituent)
#    u.add_TopologyAttr('bonds', bonds)
#    u.atoms.dimensions = [box[0], box[1], box[2], 90, 90, 90]
#    for i, ts in enumerate(u.trajectory):
#        u.atoms.unwrap(compound='fragments')
#    ca = u.select_atoms('all')
#    with mda.Writer(fileName, ca.n_atoms) as w:
#        for ts in u.trajectory:
#            ca.atoms.unwrap(compound='fragments')
#            w.write(ca)


def getRMSD(trajFile):
    """
    Returns an 1D array containing the RMSD values for each timestep. It is highly recommended to use this function on the
    trajectory file, which was created by the *cleanTrajectory*-function.

    Parameters
    ----------

    trajFile: str
        Trajectoryfile, which will be analysed.
    """
    if os.path.exists(".trajRMSD"):
        shutil.rmtree(".trajRMSD")
    os.mkdir(".trajRMSD")
    os.chdir(".trajRMSD")
    shutil.copyfile("../%s" % trajFile, trajFile)
    files = []
    firstLine = ""
    xyzString = ""
    with open(trajFile, "r") as source:
        for i, line in enumerate(source):
            if i == 0:
                firstLine = line
                xyzString = firstLine
            elif line == firstLine:
                files.append(xyzString)
                xyzString = firstLine
            else:
                xyzString += line
    files.append(xyzString)
    for i, string in enumerate(files):
        with open("%.4d.xyz" % i, "w") as target:
            target.write(string)
    rmsdList = []
    for i, string in enumerate(files):
        if i == 0:
            firstSystem = cm.System.solid.Solid(read("%.4d.xyz" % i))
        system = cm.System.solid.Solid(read("%.4d.xyz" % i))
        rmsdList.append(cm.util.rmsd.simpleRMSD(firstSystem, system))
    os.chdir("..")
    shutil.rmtree(".trajRMSD")
    return np.array(rmsdList)


def plotRMSD(rmsdList, timestep=1, filename="RMSD.pdf"):
    """
    Simple plot-function for lists returned from the getRMSD-function.

    Parameters
    ----------

    rmsdList: array
        List with the rmsd-values.
    timestep: int
        timestep for each frame in the trajectory-file (!) in picoseconds(!).
    filename: str
        Name of the resulting imagefile.
    """
    import matplotlib.pyplot as plt

    xarray = np.arange(0, len(rmsdList), 1)
    xarray = xarray * timestep
    plt.plot(xarray, rmsdList)
    xTicks = plt.xticks()
    plt.ylabel("RMSD / $\AA$")
    plt.xlabel("time in ps")
    plt.savefig(filename)


def seperateMolecules(system, box, atomsPerMolecule):
    """
    From a structure-object (or one of its inherit), which was read-in as a solid state (PBC), the molecules, which are
    partly outside the box are completed and seperated from each other. They are returned in a list of molecules.

    Only works, if the atoms are in a defined order (e.g. created with the included init-function and treated with a code,
    that does not mess up the order).

    Parameters
    ----------

    system: cm.system
        The read-in structure of the system.
    box: 1x3 array
        Dimentions of the CUBIC box. should be the same as the initilized structure
    atomsPerMolecule: int
        How many atoms are there per molecule?
    """
    expandedSystem = deepcopy(system)
    for i in range(3):
        addSystem = deepcopy(system)
        addSystem.positions[:, i] += box[i]
        expandedSystem += deepcopy(addSystem)
        addSystem.positions[:, i] -= box[i] * 2
        expandedSystem += deepcopy(addSystem)
    expandedSystem.write("expSys.xyz", format="xyz")
    cutOff = neighborlist.natural_cutoffs(expandedSystem)
    neighborList = neighborlist.NeighborList(
        cutOff, self_interaction=False, bothways=True
    )
    neighborList.update(expandedSystem)
    matrix = neighborList.get_connectivity_matrix()
    newMolecules = []
    for i in range(0, len(system), atomsPerMolecule):
        newMolecules.append([i])
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
    return retMolecules


def completeMolecules(system, box, atomsPerMolecule):
    """
    From a structure-object (or one of its inherit), which was read-in as a solid state (PBC), the molecules, which are
    partly outside the box are completed. This makes the visual evaluation of the system easyer and functions as a
    starting point for further evaluation.

    Only works, if the atoms are in a defined order (e.g. created with the included init-function and treated with a code,
    that does not mess up the order).

    Parameters
    ----------

    system: cm.system
        The read-in structure of the system.
    box: 1x3 array
        Dimentions of the CUBIC box. should be the same as the initilized structure
    atomsPerMolecule: int
        How many atoms are there per molecule?
    """
    system = seperateMolecules(system, box, atomsPerMolecule)
    solidState = cm.System.structure.Structure()
    for i in system:
        solidState += i
    return solidState
