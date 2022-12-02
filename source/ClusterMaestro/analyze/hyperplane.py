import ClusterMaestro.initilize.init as init
import ClusterMaestro.util.optimize as opt
from ClusterMaestro.util.general import frange, split
import ClusterMaestro as cm
from ase.io import read
import os, shutil, time, multiprocessing
from multiprocessing import Semaphore
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.mlab as ml
import matplotlib.gridspec as gridspec


def create_structures(
    molecule="AdPh4",
    core=None,
    substituent=None,
    distance=6.0,
    subAngle=None,
    rotAngle=None,
    xshift=0.0,
    yshift=0.0,
    distRange=None,
    rotRange=None,
    xRange=None,
    yRange=None,
    forceConst=0.5,
    dryRun=False,
    subRot=None,
    inverse=True,
):
    """Creates the structures needed for the calculation of an Potential-energy hyperplane.
    A folder "scan" will be created. If this folder already exists, the old one will be deleted
    and a new one will be created, so be careful! For the following optimization, the
    core-atoms will be fixed using "ClusterMaestro.util.optimize.xtb_constrain".

    Parameters
    ----------

    molecule: str
        Choose from fully pre-optimized molecules:
        AdPh4 (default)
    core: str
        Choose core structures See :ref:`Structure parts` for selection.
    substituent: str
        Choose substituents structures, See :ref:`Structure parts` for selection.
    distance: float
        Distance of the two monomers in \AA
    subAngle: float or None
        Defines the rotation-angle of the substituents in the finished structure.
        Very important, if the dimer will be used for a potential-energy hyperplane-scan.
    rotAngle: float or None
        Defines the relative rotation angle between the two monomers. rotAngle=0 is an alternating structure.
    xshift: float
        One dimer gets shifted in x-direction by that amount
    yshift: float
        One dimer gets shifted in y-direction by that amount
    distRange: (1x3) List or None
        If None, the value of "distance" will be used for all structures.
        Listentries: [0]: startvalue, [1]: endvalue, [2] stepsize.
    rotRange: (1x3) List or None
        If None, the value of "distance" will be used for all structures.
        Listentries: [0]: startvalue, [1]: endvalue, [2] stepsize.
    xRange: (1x3) List or None
        If None, the value of "distance" will be used for all structures.
        Listentries: [0]: startvalue, [1]: endvalue, [2] stepsize.
    yRange: (1x3) List or None
        If None, the value of "distance" will be used for all structures.
        Listentries: [0]: startvalue, [1]: endvalue, [2] stepsize.
    forceConst: float
        Force constant from xTB. Defines the force, which fixes the core-atoms of the molecules.
    dryRun: Bool
        Prints the number of resulting structures without creating them.
    subRot: (1x2) List or None
        First value: startangle; second value: endangle
        Rotates the substituents according to the rotation of the dimers. (ONLY to the rotAngle-keyword!)
    inverse: bool
        Whether or not the structure will be inversed.
    """
    if os.path.exists("scan"):
        shutil.rmtree("scan")
    os.mkdir("scan")

    system = init.dimer(
        molecule, core, substituent, distance, subAngle, rotAngle, xshift, yshift
    )
    coreAtoms = system.get_core_tags()

    structCounter = 0

    def _get_ranges(parameter):
        if parameter == None:
            return [0]
        else:
            return frange(parameter[0], parameter[1], parameter[2])

    distCount = 0
    rotCounter = 0
    rotSteps = _get_ranges(rotRange)
    for i in rotSteps:
        rotCounter += 1
    try:
        subRot.append((subRot[1] - subRot[0]) / rotCounter)
    except:
        None
    subAngleList = []
    for i in _get_ranges(subRot):
        subAngleList.append(i)
    for dist in _get_ranges(distRange):
        rotCount = 0
        for rot in _get_ranges(rotRange):
            xCount = 0
            if subRot == None:
                subAngle_use = subAngle
            else:
                subAngle_use = subAngleList[rotCount]
            for xVal in _get_ranges(xRange):
                yCount = 0
                for yVal in _get_ranges(yRange):
                    if not dryRun:
                        system = init.dimer(
                            molecule=molecule,
                            core=core,
                            substituent=substituent,
                            subAngle=subAngle_use,
                            distance=(distance if distRange == None else dist),
                            rotAngle=(rotAngle if rotRange == None else rot),
                            xshift=(xshift if xRange == None else xVal),
                            yshift=(yshift if yRange == None else yVal),
                            inverse=inverse,
                        )
                        print(system.molecules)
                        os.mkdir(
                            "scan/%.3d_%.3d_%.3d_%.3d"
                            % (distCount, rotCount, xCount, yCount)
                        )
                        system.write(
                            "scan/%.3d_%.3d_%.3d_%.3d/start.xyz"
                            % (distCount, rotCount, xCount, yCount)
                        )
                        opt.xtb_constrain(
                            "scan/%.3d_%.3d_%.3d_%.3d/const.inp"
                            % (distCount, rotCount, xCount, yCount),
                            coreAtoms,
                            forceConstant=forceConst,
                        )
                    structCounter += 1
                    yCount += 1
                xCount += 1
            rotCount += 1
        distCount += 1
    with open("hyperplane.dat", "w") as dat:
        dat.write("This is a hyperplane-calculation.\n")
        dat.write(
            "This File gives information needed by the following analysis and plot.\n"
        )
        if distRange != None:
            dat.write("\nDistance was varied\n")
            dat.write(
                "         %10.3f%10.3f%10.3f\n"
                % (distRange[0], distRange[1], distRange[2])
            )
            print(distRange)
        if rotRange != None:
            dat.write("\nRotation was varied\n")
            dat.write(
                "         %10.3f%10.3f%10.3f\n"
                % (rotRange[0], rotRange[1], rotRange[2])
            )
            print(rotRange)
        if xRange != None:
            dat.write("\ny-Range was varied\n")
            dat.write(
                "         %10.3f%10.3f%10.3f\n" % (xRange[0], xRange[1], xRange[2])
            )
            print(xRange)
        if yRange != None:
            dat.write("\nx-Range was varied\n")
            dat.write(
                "         %10.3f%10.3f%10.3f\n" % (yRange[0], yRange[1], yRange[2])
            )
            print(yRange)
        dat.write("\nStructure crestion finished")
        print("%d Structures will be created." % structCounter)


def calc_hyperplane(threads=1, acc="normal", GFN=2):
    """Calculate the hyperplane, crested using create_structures() with an xtb-routine.
    The execution will be done completely automatic with each calculation running on one core.
    The hyperplane will be chopped into chunks of 500 calculations for OS-purposes.

    Parameters
    ----------

    threads: int
        Total Number of accessible threads or cores for the parallelization.
    acc: str
        Level of accuracy. Default: normal Options: crude, sloppy, loose, lax, normal, tight, vtight, extreme
    GFN: int
        GFN-algorithm. Default: 2
    """
    startTime = time.time()
    trueRoot = os.getcwd()

    def get_xtb_string(acc=acc, GFN=GFN):
        xtbString = (
            "xtb start.xyz --opt %s -P 1 -gfn %d --input const.inp > xtb.out"
            % (acc, GFN)
        )
        return xtbString

    xtbString = get_xtb_string(acc=acc, GFN=GFN)

    try:
        with open("hyperplane.dat", "a") as target:
            target.write(
                "\n\n----------------------------------------------------------------------------"
            )
            target.write("\nxTB-command:\n      %s\n" % xtbString)
    except:
        print("No hyperplane.dat file available.\nxTB-command:\n    %s" % xtbString)

    def start_xtb(path, sema, counter, allCount):
        sema.acquire()
        os.chdir(path)
        print("         entering... %s" % os.getcwd())
        os.system(xtbString)
        print("Calculation number %d is finished." % (counter))
        os.chdir(trueRoot)
        sema.release()

    def start_all_number_sema(params=None):
        numberOfThreads = threads
        counter, calcCounter = 0, 0
        sema = Semaphore(numberOfThreads)
        paramChunks = split(params, 500)
        for chunk in paramChunks:
            jobs = []
            for param in chunk:
                path = "./scan/" + param
                if os.path.isfile(path + "xtb.out"):
                    with open(path + "xtb.out", "r") as out:
                        starter = True
                        for line in out:
                            if "finished" in line:
                                starter = False
                    if starter:
                        p = multiprocessing.Process(
                            target=start_xtb, args=(path, sema, counter, len(params))
                        )
                        jobs.append(p)
                        p.start()
                    counter += 1
                else:
                    p = multiprocessing.Process(
                        target=start_xtb, args=(path, sema, counter, len(params))
                    )
                    jobs.append(p)
                    p.start()
                    counter += 1
                print("%d calculations loaded from %d" % (counter, len(params)))
            for p in jobs:
                p.join()

    def find_fails():
        params = os.listdir("./scan")
        goList = []
        for i in params:
            starter = True
            path = "./scan/" + i
            if os.path.isfile(path + "/xtb.out"):
                with open(path + "/xtb.out", "r") as out:
                    for line in out:
                        if "finished" in line:
                            starter = False
            if starter:
                goList.append(i)
        print(len(goList))
        start_all_number_sema(goList)

    find_fails()

    print("\nTime needed for total execution: %.3f s" % (time.time() - startTime))


def create_array():
    """Function for analysis of an energy hyperplane.
    Goes through the folders and returns an array of the energies according to the hyperplane.dat-file.
    Can be used right after the "calc_hyperplane"-function.
    """
    dis, rot, x, y = False, False, False, False
    isOne = None
    with open("hyperplane.dat", "r") as dat:
        for i, line in enumerate(dat):
            if "Distance" in line:
                dis = True
                if isOne == None:
                    isOne = "dis"
            elif "Rotation" in line:
                rot = True
                if isOne == None:
                    isOne = "rot"
            elif "x-Range" in line:
                x = True
                if isOne == None:
                    isOne = "x"
            elif "y-Range" in line:
                y = True
                if isOne == None:
                    isOne = "y"
            if i == 4:
                firstParam = line.split()
                for j in range(len(firstParam)):
                    firstParam[j] = float(firstParam[j])
            if i == 7:
                secondParam = line.split()
                for j in range(len(secondParam)):
                    secondParam[j] = float(secondParam[j])
    #        minVals = getMinima()
    ran1 = int((firstParam[1] - firstParam[0]) / firstParam[2])
    ran2 = int((secondParam[1] - secondParam[0]) / secondParam[2])
    if firstParam[0] < 0 or firstParam[1] < 0:
        ran1 += 1
    if secondParam[0] < 0 or secondParam[1] < 0:
        ran2 += 1
    energyArray = np.zeros((ran1, ran2))
    for i in range(ran1):
        for j in range(ran2):
            string = "%.3d_%.3d_%.3d_%.3d" % (
                i if dis and isOne == "dis" else j if dis else 0,
                i if rot and isOne == "rot" else j if rot else 0,
                i if x and isOne == "x" else j if x else 0,
                i if y and isOne == "y" else j if y else 0,
            )
            try:
                with open("scan/%s/xtbopt.xyz" % string) as source:
                    for line in source:
                        if "energy" in line:
                            energyArray[i, j] = float(line.split()[1])
            except:
                energyArray[i, j] = np.nan
    np.savetxt("arrayFile.dat", energyArray)


def plot(label, mode="plain", zlim=None, show=False):
    """Created a simple plot of the surface, that was calculated in this folder.
    Use after using the "create_array"-function.

    Examples of the outcome can be found in Figure 1.

    .. figure:: ../../doc/source/_static/allplots_lowRes.png
    .. note::
       from ClusterMaestro.analyze import hyperplane

       hyperplane.create_structures(xRange=[-5, 5, 0.5], yRange=[-5, 5, 0.5])

       hyperplane.calc_hyperplane()

       hyperplane.create_array()

       hyperplane.plot(label="contour.png", mode="3D") OR

       hyperplane.plot(label="contour.png", mode="plain")


    Parameters
    ----------

    label: str
        Name of the file of the resulting plot.
    mode: str
        Defines the exact kind of plot. plain: simple 2D-contourplot, 3D: surfaceplot, minVals: plain + red path with minimum values.
        (See Figure 1)
    zlim: None or float
        If a surface has very high parts, it can be helpful to cut off all structure above <zlim>
        and set them to a maximum structure.
    """
    dis, rot, x, y = False, False, False, False
    isOne = None
    with open("hyperplane.dat", "r") as dat:
        for i, line in enumerate(dat):
            if "Distance" in line:
                dis = True
                if isOne == None:
                    isOne = "dis"
            elif "Rotation" in line:
                rot = True
                if isOne == None:
                    isOne = "rot"
                    print("rot")
            elif "x-Range" in line:
                x = True
                if isOne == None:
                    isOne = "x"
            elif "y-Range" in line:
                y = True
                if isOne == None:
                    isOne = "y"
            if i == 4:
                firstParam = line.split()
                for j in range(len(firstParam)):
                    firstParam[j] = float(firstParam[j])
            if i == 7:
                secondParam = line.split()
                for j in range(len(secondParam)):
                    secondParam[j] = float(secondParam[j])
    data = np.genfromtxt("arrayFile.dat")
    data -= np.nanmin(data)
    data *= 2600
    val_xdim, val_ydim = data.shape
    xax = np.zeros(val_xdim)
    yax = np.zeros(val_ydim)
    for i in range(val_xdim):
        xax[i] = firstParam[0] + i * firstParam[2]
    for i in range(val_ydim):
        yax[i] = secondParam[0] + i * secondParam[2]
    xax, yax = np.meshgrid(yax, xax)
    if zlim != None:
        lookarray = data > zlim
        data[lookarray] = zlim
    #        data [data == 0] = zlim
    if dis and isOne == "dis":
        yString = "Distance of core-centers / $\AA$"
    elif rot and isOne == "rot":
        yString = "Rotation / °"
    elif x and isOne == "x":
        yString = "Shift in x-direction / $\AA$"
    elif y and isOne == "y":
        yString = "Shift in y-direction / $\AA$"
    if dis and not isOne == "dis":
        xString = "Distance of core-centers / $\AA$"
    elif rot and not isOne == "rot":
        xString = "Rotation / °"
    elif x and not isOne == "x":
        xString = "Shift in x-direction / $\AA$"
    elif y and not isOne == "y":
        xString = "Shift in y-direction / $\AA$"
    if mode == "plain" or mode == "minVals":
        fig = plt.figure(figsize=(6.5, 5))
        ax = fig.add_subplot()
        levels = 20
        contourf_ = ax.contourf(
            xax, yax, data, vmin=np.nanmin(data), vmax=np.nanmax(data), levels=levels
        )
        ax.set_xlabel(xString)
        ax.set_ylabel(yString)
        if mode == "minVals":
            minVals = getMinima()
            ax.plot(xax[0, minVals[:, 1]], yax[minVals[:, 0], 0], color="red")
        fig.colorbar(contourf_, shrink=1, aspect=20, label="relative energy / kJ/mol")
    elif mode == "3D":
        fig = plt.figure(figsize=(6.5, 5))
        ax = fig.add_subplot(projection="3d")
        ax.plot_wireframe(xax, yax, data, lw=0.2, color="black")
        contourf_ = ax.contourf(
            xax,
            yax,
            data,
            500,
            vmin=np.nanmin(data),
            vmax=(zlim if zlim != None else np.nanmax(data)),
        )
        if zlim != None:
            ax.set_zlim(0, zlim)
            fig.colorbar(
                contourf_,
                shrink=1,
                aspect=20,
                label="relative energy / kJ/mol",
                ticks=range(0, zlim, int(zlim * 0.1)),
            )
        ax.set_xlabel(xString)
        ax.set_ylabel(yString)
        ax.set_zlabel("relative energy / kJ/mol")
        ax.azim, ax.elev = 30, 20
    if show:
        plt.show()
    else:
        plt.savefig(label, format="png", dpi=300)
    if mode == "minVals":
        plt.clf()
        energyProfile = []
        for i in range(len(minVals)):
            energyProfile.append(data[minVals[i, 0], minVals[i, 1]])
        fig, ax = plt.subplots()
        ax.plot(xax[0, :], energyProfile, color="red")
        ax.set_xlabel("Rotation / °")
        ax.set_ylabel("energy / kJ/mol")
        plt.savefig("%s_profile.png" % label[:-4], format="png", dpi=300)


def getMinima():
    """
    Finds the minimum path along the secondary axis for e.g. a scan with: fist axis: distance and second axis: rotation.
    """
    data = np.genfromtxt("arrayFile.dat")
    data -= np.nanmin(data)
    data *= 2600
    minVals = []
    for i in data.transpose():
        minVals.append(
            [
                np.argwhere(data == np.nanmin(i))[0, 0],
                np.argwhere(data == np.nanmin(i))[0, 1],
            ]
        )
    return np.array(minVals)


def writeMinima(fileName="minima.xyz"):
    """
    Writes all the Minimum structures into one xyz-File to be visulized with VMD, Avogadro or similar programs.

    Parameters
    ----------

    fileName: str
        Name of the file with the minia-structures.
    """
    minVals = getMinima()
    minIndex = [0, 0, 0, 0]
    for i in minVals:
        struct = cm.System.structure.Structure(
            read("scan/%.3d_%.3d_%.3d_%.3d/xtbopt.xyz" % (i[0], i[1], 0, 0))
        )
        struct.appendToXYZ(fileName)


def analyzeStuff(mode="rot-z", getMinStruct=False):
    """
    Analyses other important values, such as rotationbarrier and distance-differnce during the rotation.

    Parameters
    ----------

    mode: str
        Which axis were scanned? currently implemented: '*rot-z*' and '*x-y*'.
    getMinStruct: bool
        If True, writes the minimum structure in a file "minimumStructure.xyz"
    """
    data = np.genfromtxt("arrayFile.dat")
    data -= np.nanmin(data)
    data *= 2600
    with open("hyperplane.dat", "r") as dat:
        for i, line in enumerate(dat):
            if i == 4:
                distRange = line.split()
                for j in range(len(distRange)):
                    distRange[j] = float(distRange[j])
            if i == 7:
                rotRange = line.split()
                for j in range(len(rotRange)):
                    rotRange[j] = float(rotRange[j])
    if mode == "rot-z":
        print("\n\n   ***  ANALYSIS OF THE ROTATION-DISTANCE HYPERPLANE  ***\n")
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
        minima = getMinima()
        minDist = min(minima[:, 0])
        maxDist = max(minima[:, 0])
        print(
            "Minimal Distance: %.3f Angstrom" % (minDist * distRange[2] + distRange[0])
        )
        print(
            "Maximal Distance: %.3f Angstrom" % (maxDist * distRange[2] + distRange[0])
        )
        print(
            "Distance-difference during Rotation: %.3f Angstrom\n"
            % (maxDist * distRange[2] - minDist * distRange[2])
        )
        energyList = []
        for i in minima:
            energyList.append(data[i[0], i[1]])
        print("Rotation barrier: %.5f" % max(energyList))
        print(
            "Distance at the fist angle: %.3f"
            % (minima[0, 0] * distRange[2] + distRange[0])
        )
        print(
            "Distance at the middle angle: %.3f"
            % (minima[int(len(minima) / 2), 0] * distRange[2] + distRange[0])
        )
        if getMinStruct:
            minStructure = energyList.index(min(energyList))
            structs = cm.System.Trajectory("minima.xyz")
            minStruct = structs.frames[minStructure]
            minStruct.write("minimumStructure.xyz")
            print("The lowest Structure is number {}".format(minStructure))
    elif mode == "x-y":
        print("\n\n   ***  ANALYSIS OF THE ROTATION-DISTANCE HYPERPLANE  ***\n")
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
        print(
            "Overall energy difference: %.3f kJ/mol (highest point at ( %d | %d ))"
            % (
                data.max(),
                np.where(data == data.max())[0],
                np.where(data == data.max())[1],
            )
        )
