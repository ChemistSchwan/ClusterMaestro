import numpy as np
import sys, ase, os
from ase import neighborlist
import ClusterMaestro as cm
from copy import deepcopy


def getBonds(systemLength, core, substituent):
    """
    Finds the bonds for the molecules in the system. Molecules NEED to be created using the *ClusterMaestro.initilize.init*-function
    (for the order of the atoms, as this function does NOT create the bonds for a generic system, but for a defined system of molecules.)
    This is especially useful for solid state systems.

    Parameters
    ----------

    system: int
        Number of atoms in the ystem, which will be analysed.
    core: str
        Core-structure od the system.
    substituent: str
        Substituent of the system.
    """
    bonds = []
    testMol = cm.initilize.init.monomer(core=core, substituent=substituent)
    cutOff = neighborlist.natural_cutoffs(testMol)
    neighborList = neighborlist.NeighborList(
        cutOff, self_interaction=False, bothways=True
    )
    neighborList.update(testMol)
    matrix = neighborList.get_connectivity_matrix(sparse=False)
    for i, line in enumerate(matrix):
        for j, val in enumerate(line):
            if val == 1 and i < j:
                bonds.append((i, j))
    orgBonds = deepcopy(bonds)
    for i in range(len(testMol), systemLength, len(testMol)):
        for j in orgBonds:
            bonds.append((j[0] + i, j[1] + i))
    return bonds


def frange(start, stop=None, step=None):
    """A function very similar to the python-included range, but with support for floating point Numbers.

    Parameters
    ----------

    start: float
        start of the range
    stop: float
        end of the range
    step: float
        stepsize to go from start to stop
    """
    if stop == None:
        stop = start + 0.0
        start = 0.0
    if step == None:
        step = 1.0
    liste = []
    while True:
        if step > 0 and start >= stop:
            break
        elif step < 0 and start <= stop:
            break
        liste.append(start)
        #        yield float(start) # return float number
        start = start + step
    return liste


def split(arr, size):
    """Splits the given array in array of arrays, each with given size.

    Parameters
    ----------

    arr: array
        Starting array
    size: int
        Size of inner array in return.
    """
    arrs = []
    while len(arr) > size:
        pice = arr[:size]
        arrs.append(pice)
        arr = arr[size:]
    arrs.append(arr)
    return arrs


def chunks(l, n):
    """
    Splits a list *l* in a list of lists, each with size *n*.

    Parameters
    ----------

    l: list
        Starting list.
    n: int
        Size of each chunk
    """
    for i in range(0, len(l), n):
        yield l[i : i + n]


def chunkAtoms(atoms, n):
    """
    Splits an Atoms-Object *atoms* in a list of lists, each with size *n*.

    Parameters
    ----------

    l: list
        Starting list.
    n: int
        Size of each chunk
    """
    retList = []
    for i in range(0, len(atoms), n):
        # retList.append([atoms[i:i + n]])
        retList.append([])
        for j in range(n):
            retList[-1].append(atoms[i + j])
    return retList


def percent(val1, val2):
    """Returns the difference of the teo given values in percent.

    Parameters
    ----------

    val1: float
        First value
    val2: float
        Second value
    """
    percent = val1 - val2
    percent /= (val1 + val2) / 2
    return percent * 100


def distance(atom1, atom2):
    """
    Calculates the distance of the two given atoms.

    Parameters
    ----------

    atom1: Atom

    atom2: Atom
    """
    dist = np.linalg.norm(atom1.position - atom2.position)
    return dist


def fileSplitter(xyzFile):
    """
    This function splits a file with multiple xyz-Files in it into seperate xyz-files for further analysis.

    Parameters
    ----------

    xyzFile: str
        File, which will be splitted.
    """
    files = []

    firstLine = ""
    xyzString = ""

    with open(dataFile, "r") as source:
        for i, line in enumerate(source):
            if i == 0:
                firstLine = line
                xyzString = firstLine
            elif line == firstLine:
                files.append(xyzString)
                xyzString = firstLine
            else:
                xyzString += line

    for i, string in enumerate(files):
        with open("split_%d.xyz" % i, "w") as target:
            target.write(string)


def isFloat(something):
    try:
        float(something)
        return True
    except:
        return False


def isInt(something):
    try:
        int(something)
        return True
    except:
        return False


def findClosest(extendSystem, centerPoint, dontLook):
    """
    Find the closest Atom to *centerPoint* in *extendSystem* and ignore all atoms, whose tag are in *dontLook*.
    """
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


def deleteSecondLine(dataFile):
    """
    Deletes the second line of a given file. This might come in handy, as this line can cause problems for
    an xyz file, if the file is created and read with different codes.

    Parameters
    ----------

    dataFile: str
        File to be read in and modified.
    """
    with open(dataFile, "r") as source:
        with open("{}_tmp".format(dataFile), "w") as target:
            for i, line in enumerate(source):
                if i == 1:
                    target.write("\n")
                else:
                    target.write(line)
    os.remove(dataFile)
    os.rename("{}_tmp".format(dataFile), dataFile)


def deleteLastLine(dataFile):
    """
    Deletes the last line of a given file. This might come in handy, as this line can cause problems for
    an xyz file, if the file is created and read with different codes.

    Parameters
    ----------

    dataFile: str
        File to be read in and modified.
    """

    with open(dataFile, "r") as source:
        lines = source.readlines()
        lines = lines[:-1]
    with open("{}_tmp".format(dataFile), "w") as target:
        for line in lines:
            target.write(line)
    os.remove(dataFile)
    os.rename("{}_tmp".format(dataFile), dataFile)
