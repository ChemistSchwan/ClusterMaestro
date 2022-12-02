import ClusterMaestro as cm
from copy import deepcopy

def shiftMols(trajectory, atomsPerMol):
    """
    Cleans the remaining molecule-jumps from the *cleanTrajectory* Function.

    Parameters
    ----------

    atomsPerMol: int
        How many atoms do we have per dedecated molecule in the system?
    """
    for i in range(len(trajectory.frames)):
        if i != 0:
            if cm.util.rmsd.simpleRMSD(trajectory.frames[i], trajectory.frames[i-1]) > 2:
                print("\n")
                print("Shift of a molecule from timestep %d to %d" % (i-1, i))
                print("Running analyze-function...")
                sys1 = deepcopy(trajectory.frames[i])
                sys2 = deepcopy(trajectory.frames[i-1])
                sys1.analyze(atomsPerMol)
                sys2.analyze(atomsPerMol)
                for j in range(len(sys1.molecules)):
                    if cm.util.rmsd.simpleRMSD(sys1.molecules[j], sys2.molecules[j]) > 1:
                        centerOfMasses = sys1.molecules[j].center_of_mass() - sys2.molecules[j].center_of_mass()
                        for k, val in enumerate(centerOfMasses):
                            if abs(val) > 5:
                                if val > 0:
                                    for l in range(i, len(trajectory.frames)):
                                        trajectory.frames[l].positions[j*atomsPerMol:(j+1)*atomsPerMol,k] -= trajectory.box[k]
                                if val < 0:
                                    for l in range(i, len(trajectory.frames)):
                                        trajectory.frames[l].positions[j*atomsPerMol:(j+1)*atomsPerMol,k] += trajectory.box[k]
                                print("Shift of molecule %d was corrected." % j)
    return trajectory

