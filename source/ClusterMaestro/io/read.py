import ClusterMaestro as cm
import numpy as np
import ase


def readLmpTraj(filename, timestep=1):
    """
    Read in a lammps-dump of the format:

    dump d1 all     custom 10000 traj.lmp element xu yu zu

    Parameters
    ----------

    filename: str
        Name of the datafile.
    timestep = float
        Time between steps. Not in the trajectory file, needs to be added manually.
    """
    traj = cm.System.trajectory.Trajectory()
    init, bound, getAtoms, boundCounter, atomCounter, getNumber = (
        True,
        False,
        False,
        0,
        0,
        False,
    )
    boundaries = np.zeros((3))
    firstLine = ""
    with open(filename, "r") as source:
        for line in source:
            if line == firstLine:
                getAtoms = False
                atomSystem = ase.Atoms(
                    symbols=labelList, positions=atomCoords, cell=boundaries
                )
                atomSystem = cm.System.solid.Solid(atomSystem)
                atomNumber = 0
                traj.frames.append(atomSystem)
            if init:
                atomSystem = ase.Atoms()
                firstLine = line
                init = False
                getAtoms = False
                atomNumber = 0
            if getNumber:
                atomCoords = np.zeros((int(line), 3), dtype=np.float64)
                getNumber = False
            if "NUMBER OF ATOMS" in line:
                getNumber = True
            if getAtoms:
                splitLine = line.split()
                labelList.append(splitLine[0])
                atomCoords[atomNumber] = np.array(
                    [float(splitLine[1]), float(splitLine[2]), float(splitLine[3])]
                )
                atomNumber += 1
            if "ATOMS element" in line:
                getAtoms, labelList = True, []
            if bound:
                boundaries[boundCounter] = float(line.split()[1]) - float(
                    line.split()[0]
                )
                boundCounter += 1
                if boundCounter == 3:
                    bound, boundCounter = False, 0
            if "BOX BOUNDS" in line:
                bound = True
        atomSystem = ase.Atoms(symbols=labelList, positions=atomCoords, cell=boundaries)
        atomSystem = cm.System.solid.Solid(atomSystem)
        traj.frames.append(atomSystem)
    traj.timeStep = timestep
    return traj
