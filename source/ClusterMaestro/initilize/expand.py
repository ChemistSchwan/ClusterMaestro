import inspect, os, math
from ase.io import read
from scipy.ndimage import rotate
from scipy.linalg import expm
import ClusterMaestro as cm
import numpy as np
import sys, os


def core(fileName, structureName, delete=False, opt=False, level=None):
    """
    Takes the given core-structure and stores it in the database for future usage of initilization.

    .. Attention::
       The Atoms, where the substituents will attach MUST be covered with Flour-atoms. They will be removed, if a molecule is created.
       For a better understanding, look at the structures, that are already in the database.

    .. Warning::
       If you want to use these structures for the creation of Dimers, take attention regarding the orientation of the structure.
       Always doublecheck the structure!

    Parameters
    ----------

    fileName: str
        Name of the xyz-file, which contains the structure
    structureName: str
        Name of the structure, as it will be stored in the database.
    delete: Bool
        Whether or not the old structure with the same name will be removed.
    opt: Bool
        If *True*, optimizes the given structure using GFN2-xTB
    level: str or None
        Must be x, y, z (or type None). Given structure will be leveled to match the x-, y- or z-plane.
    """
    path = inspect.getfile(cm).replace("__init__.py", "") + "initilize/structures/cores"
    if os.path.exists("%s/%s.xyz" % (path, structureName)) and not delete:
        print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("++                     WARNING                     ++")
        print("++  A structure with this name already exists!     ++")
        print("++            Please remove structurefile          ++")
        print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("Structurefile: \n%s/%s.xyz" % (path, structureName))
    else:
        if os.path.exists("%s/%s.xyz" % (path, structureName)):
            os.remove("%s/%s.xyz" % (path, structureName))
        core = cm.System.core.Core(read(fileName))
        if opt:
            core = cm.util.optimize.opt_xtb(core, acc="extreme")
        center = core.center_of_mass()
        core.positions -= center
        if level != None:
            XYZ = core.positions
            p0 = [0, 0, 1, 0]

            def f_min(X, p):
                plane_xyz = p[0:3]
                distance = (plane_xyz * X).sum(axis=1) + p[3]
                return distance / np.linalg.norm(plane_xyz)

            def residuals(params, signal, X):
                return f_min(X, params)

            from scipy.optimize import leastsq

            sol = leastsq(residuals, p0, args=(None, XYZ))[0]
            if level == "x":
                p1 = [1, 0, 0, 0]
            elif level == "y":
                p1 = [0, 1, 0, 0]
            elif level == "z":
                p1 = [0, 0, 1, 0]
            else:
                print("level must be x, y, z, or None!")
            try:
                ang = cm.util.linalg_lib.planeAngle(sol, p1)
                axis = cm.util.linalg_lib.planeIntersect(sol, p1)
                for i, atom in enumerate(core):
                    core[i].position = cm.util.linalg_lib.rotate(
                        atom.position, axis, ang
                    )
            except:
                print(
                    "Molecule can not be rotated, if it is already perfectily alligned."
                )
                print("Please shift a few coordinates a little bit.")
        core.write("%s/%s.xyz" % (path, structureName), format="xyz")


def substituent(fileName, structureName, delete=False, opt=False):
    """
    Takes the given substituent-structure and stores it in the database for future usage of initilization.

    .. Attention::
       The Carbon, that belongs to the framestructure must be the frist atom in the file.
       The atom attached to this carbon must be the second one.

    .. Warning::
       If you want to use these structures for the creation of Dimers, take attention regarding the orientation of the structure.
       **Always doublecheck the structure!**

    Parameters
    ----------

    fileName: str
        Name of the xyz-file, which contains the structure
    structureName: str
        Name of the structure, as it will be stored in the database.
    delete: Bool
        Whether or not the old structure with the same name will be removed.
    opt: Bool
        If *True*, optimizes the given structure using GFN2-xTB
    """
    path = (
        inspect.getfile(cm).replace("__init__.py", "")
        + "initilize/structures/substituents"
    )

    def findVector(vector1, vector2):
        """
        Returns the vector from point "vector1" to point "vector2"
        """
        vector = np.zeros(3)
        for i in range(len(vector)):
            vector[i] = vector1[i] - vector2[i]
        return vector

    def rotateVector(vector1, vector2, vector3):
        """
        Rotates vector 3 to match the direction of vector1. Vector2 is the starting orientation of the fragment.
        """
        c = np.dot(vector1, vector2) / np.linalg.norm(vector1) / np.linalg.norm(vector2)
        angle = np.arccos(np.clip(c, -1, 1))
        cross = np.cross(vector1, vector2)
        M0 = expm(np.cross(np.eye(3), cross / np.linalg.norm(cross) * angle))
        return np.dot(M0, vector3)

    def main(fragmentFile):
        molecule = read(fragmentFile)
        molecule.positions -= molecule.positions[0]
        CX_vec = findVector(molecule.positions[1], molecule.positions[0])
        for i in range(len(molecule.positions)):
            molecule.positions[i] = rotateVector(
                CX_vec, np.array([0, 0, 1]), molecule.positions[i]
            )
        molecule.write(fragmentFile)

    if os.path.exists("%s/%s.xyz" % (path, structureName)) and not delete:
        print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("++                     WARNING                     ++")
        print("++  A structure with this name already exists!     ++")
        print("++            Please remove structurefile          ++")
        print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("Structurefile: \n%s/%s.xyz" % (path, structureName))
    else:
        main(fileName)
        with open(fileName, "r") as source:
            with open("tmp", "w") as target:
                counter = 0
                for line in source:
                    if counter == 0:
                        target.write("%d\n" % (int(line.split()[0]) - 1))
                    elif counter != 2:
                        target.write(line)
                    counter += 1
        os.rename("tmp", "{}/{}.xyz".format(path, structureName))
