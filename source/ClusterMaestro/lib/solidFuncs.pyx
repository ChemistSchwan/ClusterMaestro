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
                sys1.analyzeFAST(atomsPerMol)
                sys2.analyzeFAST(atomsPerMol)
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


def analyze(solidState, atomsPerMolecule, writeStructures=False):
    """
    Analyses the structure by running the *seperateMolecules*-function. Consequently, this function analyses the molecules
    of the given solid state. A short summary of the properties of the given system will be printed.

    Parameters
    ----------

    atomsPerMolecule: int
        Number of atoms in each molecule.
    """
    solidState.seperateMolecules(atomsPerMolecule)
    density = solidState.calcDensity()
    print("The Density of the system is %.5f g/cm^3\n" % density)
    atomsPerMolList = []
    warningBool = False
    string1, string2 = "%18s" % "Molecule:", "%18s" % "Number of Atoms:"
    for i, molecule in enumerate(solidState.molecules):
        string1 += "%5d" % i
        string2 += "%5d" % len(molecule)
        if len(molecule) != atomsPerMolecule:
            warningBool = True
    print(string1)
    print(string2)
    if warningBool:
        print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("++                            WARNING                            ++")
        print("++        Not all molecules have the right number of atoms       ++")
        print("++   Check for harsh structure distortions/ unwanted reactions   ++")
        print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    if writeStructures:
        completeStructure = solidState.completeMolecules(atomsPerMolecule)
        completeStructure.write("completeStructure.xyz", format="xyz")
        cm.util.write.poscar("completeStructure.vasp", completeStructure,
                box=[solidState.cell[0,0], solidState.cell[1,1], solidState.cell[2,2]])
        for i, molecule in enumerate(solidState.molecules):
            molecule.write("mol_%.3d.xyz" % i, format="xyz")
    return solidState

def structureFactor(solid):
    """
    Calculates the structure factor. This function may take a while. Saves the resulting array in a file
    named "sf.npy".

    Parameters
    ----------

    plotName: str
        Name of the resulting figure.
    """
    import pylab as pl
    import matplotlib.pyplot as plt
    import numpy as np
    def ssf(x, size, q, pbc=False):
        """
        From a series of positions x in a cubic box of length size we get
        the structure factor for momentum q
        """
        print("Another Iteration done")
        natoms = np.shape(x)[0]
        sf = 0.0
        for i in range(natoms):
            x1 = x[i]
            for j in range(i+1, natoms):
                x2 = x[j]
                dx = x2 - x1
                if pbc:
                    for i in range(3):
                        if dx[i] >  size/2: dx[i] -= size
                        if dx[i] < -size/2: dx[i] += size
                r = np.linalg.norm(dx)
                sf += 2*np.sin(q*r)/(q*r)
        sf /= natoms
        sf += 1
        return sf
    q = np.linspace(0, 11.0, 56)[1:]
    size = solid.cell[0,0]
    x = solid.positions
    sf = [ssf(x, size, _, True) for _ in q]
    sf = np.array(sf)
    return sf

