import ClusterMaestro as cm
import numpy as np
from ase.io import read
import os, linecache


def subContanctTypes(molecule1, molecule2, distance=6.5):
    """
    Calculates the number of substituent contacts of type 1, 2 and 3 in a dimer. All substituents taken into account must have a distance
    (regarding their center of mass) of less than *distance* and additionally:

    Type 1: angle < 10°

    Type 2: angle > 80°

    Type 3: angle > 10° and < 80°

    Parameters
    ----------
    mnolecule1: Molecule
        Molecule structure as input. Must be analyzed already.
    mnolecule2: Molecule
        Molecule structure as input. Must be analyzed already.
    distance: float
        Maximum distance 2 substituents can be apart for beeing taken into account for this analysis.

    Return
    ------
    1x3 array with Type 1 contact in first, Type 2 ccontact in second and Type 3 contact in third place.
    """
    typeArray = np.zeros((3))
    for i, sub1 in enumerate(molecule1.substituents):
        for j, sub2 in enumerate(molecule2.substituents):
            distance = np.linalg.norm(sub1.center_of_mass() - sub2.center_of_mass())
            if distance < 6.5:
                plane1 = cm.util.linalg_lib.planeFit(sub1.positions)[0]
                plane2 = cm.util.linalg_lib.planeFit(sub2.positions)[0]
                angle = cm.util.linalg_lib.planeAngle(plane1, plane2)
                if angle < 10:
                    typeArray[0] += 1
                elif angle > 80 and angle < 90:
                    typeArray[1] += 1
                else:
                    typeArray[2] += 1
    return typeArray


def getFurthestAtoms(molecule1, molecule2):
    """
    Calculates the distance between the two atoms, which are the furthest away from each
    other minus the closest atom distance.
    This can give a view on the general orientation of the dimers.

    Parameters
    ----------
    mnolecule1: Molecule
        Molecule structure as input. Must be analyzed already.
    mnolecule2: Molecule
        Molecule structure as input. Must be analyzed already.
    """
    closest_dist, furthest_dist = 100, 0
    for atom1 in molecule1:
        for atom2 in molecule2:
            dist = np.linalg.norm(atom1.position - atom2.position)
            if dist > furthest_dist:
                furthest_dist = dist
            elif dist < closest_dist:
                closest_dist = dist
    return furthest_dist - closest_dist


def analyzeDimers(
    folderName,
    subDistance=6.5,
    featureset=["subContanctTypes", "furthestAtoms", "coreDistance"],
):
    """
    The folder should consist of all the dimer structures to be analyzed.

    This function is parallelized and will use half of the available threads, which usually matches the number of cpu cores.

    Parameters
    ----------
    folderName: str
        Name of the folder with the dimer structures in xyz format.
    subDistance: float
        See *distance* in *ClusterMaestro.analyze.featureCalc.subContanctTypes*
    featureset: list of strings
        Which features are to be calculated? Choose from: *subContanctTypes*, *furthestAtoms*, *coreDistance*
    """
    from joblib import Parallel, delayed
    import multiprocessing

    files = os.listdir(folderName)
    for i in range(len(files)):
        files[i] = "%s/%s" % (folderName, files[i])
    fileList = files
    fileList.sort()

    headerString = "%15s%15s" % ("filename", "energy")
    for feature in featureset:
        if "subContanctTypes" in featureset:
            headerString += "%15s%15s%15s" % ("first-Sub", "second-Sub", "third-Sub")
        if "coreDistance" in featureset:
            headerString += "core-core"
        if "furthestAtoms" in featureset:
            headerString += "totalDistDiff"
    headerString += "\n"

    with open("features.dat", "w") as target:
        target.write(headerString)

    def analyzeDimer(filename):
        with open("features.dat", "a") as target:
            system = cm.System.Oligomer(read(filename))
            system.analyze(25)
            system.energy = float(linecache.getline(filename, 2))
            target.write("%15s%15.5f" % (filename[-7:], system.energy))
            for feature in featureset:
                if "subContanctTypes" in featureset:
                    contactArray = subContanctTypes(
                        system.molecules[0], system.molecules[1], subDistance
                    )
                    target.write(
                        "%15d%15d%15d"
                        % (contactArray[0], contactArray[1], contactArray[2])
                    )

                if "coreDistance" in featureset:
                    coreDistance = system.get_core_distance()
                    target.write("%15.5f" % coreDistance)

                if "furthestAtoms" in featureset:
                    dist = cm.analyze.featureCalc.getFurthestAtoms(
                        system.molecules[0], system.molecules[1]
                    )
                    target.write("%15.5f" % dist)
            target.write("\n")

    coreNumber = int(multiprocessing.cpu_count() / 2)
    Parallel(n_jobs=coreNumber)(
        delayed(analyzeDimer)(filename) for filename in fileList
    )
