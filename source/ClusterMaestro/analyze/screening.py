"""
Library for the fully-automated screening of different properties of single clusters and cluster dimers.
"""

import os, shutil, time
from copy import deepcopy
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
from ase.io import read
from ase_own.calculators.turbomole import Turbomole
import ClusterMaestro as cm
from ClusterMaestro.analyze.plotFunctions import *


def crest_opt(
    core="Ad",
    substituent="Me",
    coreDist=6,
    crest=False,
    threads=1,
    dftFiles=[],
    startOpt=False,
    optScript="",
    params={},
    justMono=False,
):
    """
    Automated analysis including:
    Crest-conformer analysis,

    Parameters
    ----------

    core: str
        Desired core-structure, See :ref:`Structure parts`.
    substituent: str
        Desired substituent-structure, See :ref:`Structure parts`.
    coreDist: float
        Desired core-core Distance, see ClusterMaestro.initilize.
    crest: bool
        Should a crest run be performed or read in? If the run was already performed and should
        only be read in, the pure runs MUST be in a folder called *"monomer"* and *"dimer"*!
    threads: int
        How many threads should the crest-run use?
    dftFiles: list of strings
        Files, which are needed for the optimization on dft level of theory. They will be copied
        into the *"opt"* folders within the crest-run. The slurm script has to be named
        "*job.slurm*", if the optimization will be started right away.
    startOpt: bool
        Whether or not the optimization will be started. This would mean, that the order
        *"sbatch job.slurm"* is executed.
    optScript: str
        Simple slurm script, which will start the optimization.
    params: dict
        Set of parameters, as beeing used by the *ase.calculator.Turbomole* package.
    justMono: bool
        If True, the dimer structure will not be created or calculated.
    """
    trueRoot = os.getcwd()
    if crest:
        os.mkdir("monomer")
        mono = cm.initilize.init.monomer(core=core, substituent=substituent)
        mono.write("monomer/start.xyz")
        os.chdir("monomer")
        os.system("crest start.xyz -T {} > crest.out".format(threads))
        if not justMono:
            os.chdir(trueRoot)
            os.mkdir("dimer")
            dimer = cm.initilize.init.dimer(
                core=core, substituent=substituent, distance=coreDist
            )
            dimer.write("dimer/start.xyz")
            os.chdir("dimer")
            os.system("crest start.xyz -T {} > crest.out".format(threads))
        os.chdir(trueRoot)

    if justMono:
        folderLists = ["monomer"]
    else:
        folderLists = ["monomer", "dimer"]
    for folder in folderLists:
        os.chdir(folder)
        os.mkdir("opt")
        shutil.copyfile("crest_best.xyz", "opt/start.xyz")
        shutil.copyfile("../{}".format(optScript), "opt/{}".format(optScript))
        for dftfile in dftFiles:
            shutil.copyfile("{}/{}".format(trueRoot, dftfile), "opt/{}".format(dftfile))
        if startOpt:
            os.chdir("opt")
            calc = Turbomole(**params)
            mol = read("start.xyz")
            mol.set_calculator(calc)
            calc.initialize()
            os.system("chmod u+x {}".format(optScript))
            os.system("./{}".format(optScript))
        os.chdir(trueRoot)


def bindingInfo(
    structureFile="coord",
    startBind=False,
    dftFiles=[],
    spScript="",
    coreAtoms=[],
    params={},
    fromXYZ=False,
    analyzedDimer=None,
):
    """
    Prepares the calculation of the binding energies for the dimer structure.

    Parameters
    ----------

    structureFile: str
        Name of the structure file of the optimized structure in the *dimer* folder.
    startBind: bool
        Whether or not the singlepoint calculations will be started. This would mean, that
        the order *"sbatch jobBind.slurm"* is executed.
    spScript: str
        Startscript for the singlepoint calculations. Files will be copied to each folder for
        the execution.
    params: dict
        Set of parameters, as beeing used by the *ase.calculator.Turbomole* package.
    fromXYZ: bool
        Whether or not the structures are already prepared.
    analyzedDimer: Dimer
        Dimer, which was already analyzed.
    """
    trueRoot = os.getcwd()
    print(trueRoot)
    if not fromXYZ:
        os.mkdir("bindingEnergies")
        inputFile = "../dimer/opt/{}".format(structureFile)
        os.chdir("bindingEnergies")
        if not "xyz" in inputFile:
            cm.analyze.bonding.initStructureParts(
                readInput=inputFile,
                inputType="turbomole",
                coreAtoms=coreAtoms,
                analyzedDimer=analyzedDimer,
            )
        else:
            cm.analyze.bonding.initStructureParts(
                readInput=inputFile, analyzedDimer=analyzedDimer
            )
    #    os.chdir("..")
    else:
        os.chdir("bindingEnergies")
    if startBind:
        folderList = [
            "dimer",
            "mono1",
            "mono2",
            "core",
            "core1",
            "core2",
            "substituent",
            "substituent1",
            "substituent2",
        ]
        for slurmfile in dftFiles:
            for i in folderList:
                shutil.copyfile(
                    "{}/{}".format(trueRoot, slurmfile), "{}/{}".format(i, slurmfile)
                )
        for i in folderList:
            os.chdir(i)
            shutil.copyfile("../../{}".format(spScript), spScript)
            mol = cm.System.Structure(read("start.xyz"))
            elemList = list(set(mol.get_chemical_symbols()))
            for i in range(len(elemList)):
                elemList[i] = elemList[i].lower()
            deleteECP = True
            params_2 = deepcopy(params)
            for i in params_2:
                if "ecp" in i:
                    try:
                        if params_2[i].split()[0] in elemList:
                            deleteECP = False
                    except:
                        for j in params_2[i]:
                            if j.split()[0] in elemList:
                                deleteECP = False
            if deleteECP:
                try:
                    params_2.pop("ecp name")
                    params_2.pop("ecp file")
                    params_2.pop("basis set atom")
                    params_2.pop("basis set file")
                except:
                    pass
            try:
                calc = Turbomole(**params_2)
                mol.set_calculator(calc)
                calc.initialize()
            except:
                params_2["use redundant internals"] = False
                calc = Turbomole(**params_2)
                mol.set_calculator(calc)
                calc.initialize()
            os.system("chmod u+x {}".format(spScript))
            os.system("./{}".format(spScript))
            os.chdir("..")
    os.chdir(trueRoot)
    time.sleep(2)
    print(trueRoot)


def bindingAnalysis(calctype="tm", plot=True):
    """
    Analyses the binding energy, writes the datafile and plots a simple plot.
    Keep in mind to delete the *slurm.out* files with a find function for full automation.

    Parameters
    ----------

    calctype: str
        See *ClusterMaestro.analyze.bonding.analyseSinglePoint*
    """
    try:
        trueRoot = os.getcwd()
        print(trueRoot)
        os.chdir("bindingEnergies")
        cm.analyze.bonding.analyseSinglePoint(calcType=calctype, writeFile=True)
        if plot:
            cm.analyze.bonding.plotBondshare(width=5, height=3.2, legendloc=None)
        os.chdir(trueRoot)
        time.sleep(2)
    except:
        pass
    os.chdir(trueRoot)


def crestStart(structureFile="dimer/crest_best.xyz", crestScript=""):
    """
    Prepares the calculation for a total crest-analysis.
    Be aware, that the energy window (*-ewin*) should be set manually to e.g. 30 kcal.

    Parameters
    ----------

    structureFile: str
        Which structure should be used?
    crestScript: str
        Startscript for the crest calculations. Files will be copied to each folder for the
        execution.
    """
    os.mkdir("crestAnalysis")
    shutil.copyfile(structureFile, "crestAnalysis/start.xyz")
    shutil.copyfile(crestScript, "crestAnalysis/{}".format(crestScript))
    os.chdir("crestAnalysis")
    os.system("chmod u+x {}".format(crestScript))
    os.system("./{}".format(crestScript))
    os.chdir("..")


def crestAnalysis(maxRMSD=1.0, symTolerance=0.5):
    """
    Analyses the mass-crest routine. If you want to change the parameters, delete the folders
    *crestAnalysis/rmsdGroups* and *crestAnalysis/symmetryGroups* before restarting the analysis.

    Parameters
    ----------

    maxRMSD: float
        See *ClusterMaestro.analyze.crestutils.sortCRESTrmsd*
    symTolerance: float
        See *ClusterMaestro.analyze.crestutils.sortCRESTsym*
    """
    os.chdir("crestAnalysis")
    if not os.path.exists("allStructures"):
        cm.analyze.crestutils.splitter()
    if not os.path.exists("rmsdGroups"):
        cm.analyze.crestutils.sortCRESTrmsd(maxRMSD=maxRMSD)
    if not os.path.exists("symmetryGroups"):
        cm.analyze.crestutils.sortCRESTsym(tolerance=symTolerance)
    cm.analyze.crestutils.plotSymEng()
    cm.analyze.crestutils.dimerBondEng(data="sym", monomer="../monomer/crest_best.xyz")
    cm.analyze.crestutils.dimerSymData()
    cm.analyze.crestutils.plotDimerSymData(
        height=3, width=4, xRange=(5, 12), legend=False
    )
    try:
        data = np.genfromtxt("bondlength.data")
    except:
        data = cm.analyze.crestutils.getCoreBondLengths()
        np.savetxt("bondlength.data", data)
    plt.clf()
    plt.rcParams.update({"figure.figsize": (4, 3), "figure.dpi": 100})
    plt.hist(data, bins=50)
    plt.gca().set(xlabel="bondlength / $\AA$")
    mu, std = norm.fit(data)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, "k", linewidth=2)
    title = "Fit results: mu = {:.3f},  std = {:.5f}".format(mu, std)
    plt.xlim(mu - 0.1, mu + 0.1)
    plt.title(title)
    plt.tight_layout()
    plt.yticks([], [])
    plt.savefig("bondLengthDist.pdf")
    os.chdir("..")


def getAllInfo(graphics=False, optmode="tm", core_cutoff=25):
    """
    Collects all information generated by the functions: *crest_opt*, *bindingInfo*,
    *bindingAnalysis*, *crestStart* and *crestAnalysis*. It is recommended to run all
    the mentioned functions FIRST and run this function afterwords to get a small
    summary of the calculations.

    A new folder will be created with the most important Information nice and compact, all in
    one place.

    Parameters
    ----------

    graphics: bool
        Whether or not the graphics are copied into the new folder as well, or not.
    optmode: str
        Wich program was used for the optimizations. Currently available: *tm*: turbomole
    core_cutoff: int
        See *ClusterMaestro.System.monomer*
    """

    def getED_tm(path, copy=True, name=False):
        if not name:
            for i in os.listdir(path):
                if ".out" in i or ".log" in i:
                    dataFile = i
                    dataPath1 = "{}/{}".format(path, dataFile)
                    with open(dataPath1, "r") as testFile:
                        for line in testFile:
                            if "RUNNING PROGRAM" in line:
                                dataPath = "{}/{}".format(path, dataFile)
                                break
        else:
            dataPath = path
        energy, dispersion = None, None
        with open(dataPath, "r") as source:
            for line in source:
                if "total energy" in line:
                    energy = line.split()[4]
                if "empirical dispersive energy correction" in line:
                    dispersion = line.split()[-1]
        if copy:
            shutil.copyfile(dataPath, "archive/{}/{}.out".format(systemName, path))
        if type(energy) != str:
            return None
        else:
            return [float(energy), float(dispersion)]

    if optmode == "tm":
        dimerMin = cm.System.Oligomer(read("dimer/opt/coord", format="turbomole"))
        monoMin = cm.System.Oligomer(read("monomer/opt/coord", format="turbomole"))
        dimerMin.energy = getED_tm("dimer/opt", copy=False)[0] * 0.036749308137
        monoMin.energy = getED_tm("monomer/opt", copy=False)[0] * 0.036749308137
    dimerMin.analyze(core_cutoff)
    print("+++++++++++++++++++++++++++++++++++++++")
    print(dimerMin)
    print(dimerMin.energy)
    print(dimerMin.get_core_distance())
    print(monoMin)
    print(monoMin.energy)
    print(
        "Dissociation energy: {:.5f} kJ/mol".format(
            (dimerMin.energy - 2 * monoMin.energy) * 2625.4996427573383
        )
    )


def moreInfo(dftFiles=[], getOrbitals=False, params={}):
    """
    Calculates additional information using turbomole. The calculations must be written in the
    *more.sh*-file, so this program can find it. An optimization must be performed beforehand!

    Parameters
    ----------

    dftFiles: list of strings
        Files, which are needed for the optimization on dft level of theory. They will be
        copied into the *"opt"* folders within the crest-run. The slurm script has to be named
        "*job.slurm*", if the optimization will be started right away.
    params: dict
        Set of parameters, as beeing used by the *ase.calculator.Turbomole* package.
    """
    trueRoot = os.getcwd()
    params["scf energy convergence"] = 8
    for i in ["monomer", "dimer"]:
        os.mkdir("{}/opt/more".format(i))
        shutil.copyfile("more.sh", "{}/opt/more/more.sh".format(i))
        for dftfile in dftFiles:
            shutil.copyfile(
                "{}/{}".format(trueRoot, dftfile), "{}/opt/more/{}".format(i, dftfile)
            )
        system = read("{}/opt/coord".format(i))
        system.write("{}/opt/more/start.xyz".format(i), format="xyz")
        os.chdir("{}/opt/more".format(i))

        calc = Turbomole(**params)
        system.set_calculator(calc)
        calc.initialize()

        os.system("chmod u+x more.sh")
        os.system("./more.sh")
        os.chdir(trueRoot)


def createRotationStructs(
    system,
    mode="tm",
    writeXYZ=True,
    rotrange=[0, 360, 10],
    dftFiles=[],
    params={},
    threads=1,
    substituent=0,
):
    """
    Creates the structures for a scan of the rotation of the first substituent in a molecule.
    Creates a folder for the final structures and initilized the structures using the
    ASE-*calculator* object. I recomment my modified ASE-package for that, as constrained
    optimizations are not included in the normal ASE-package.

    For the actual calculation, a start file named *rotation.sh* has to be in the startfolder.
    This will be handled just like all the other calculations (more, crest...).

    Parameters
    ----------

    system: ClusterMaestro.System.Molecule
        System, where the substituent will be rotated.
    mode: str
        Currently only *tm*
    calculator: ASE.Calculator or None
        Defines the calculation parameters using the ASE-Calculator object.
    writeXYZ: bool
        Whether or not to write the created structures into a xyz-file for easy inspection.
    rotrange: 1x3 list
        Start-, end and step value for the rotation angle.
    dftFiles: list of str
        Additional files needed for the DFT calculation, e.g. pseudopotential.
    threads: int
        Number of parallelization used. This defines the numbers of optimization calculations run
        simultaniously.
    substituent: int
        Number of the substituent, that is beeing rotated. This is the number, in which the
        substituent is written in the sytem.substituents list.
    """
    import multiprocessing
    from multiprocessing import Semaphore

    def getAngleAndAtoms(system):
        tagList = np.array([0, 0, 0, 0])
        print(system.substituents[substituent][0])
        tagList[2] = system.substituents[substituent][0].tag
        connectivity = system.substituents[substituent].connectivity_matrix()
        for index, i in enumerate(connectivity[0]):
            if i == 1 and system.substituents[substituent][index].tag not in tagList:
                tagList[3] = system.substituents[substituent][index].tag
                break

        connectMatrix = system.connectivity_matrix()
        for index, i in enumerate(connectMatrix[tagList[2]]):
            if i == 1 and system[index].tag not in tagList:
                tagList[1] = system[index].tag
                break

        for index, i in enumerate(connectMatrix[tagList[1]]):
            if i == 1 and system[index].tag not in tagList:
                tagList[0] = system[index].tag
                break

        angleArray = [
            system[tagList[0]].position,
            system[tagList[1]].position,
            system[tagList[2]].position,
            system[tagList[3]].position,
        ]
        angleArray = np.array(angleArray)
        angle = cm.util.linalg_lib.dihedralAngle(
            angleArray[0], angleArray[1], angleArray[2], angleArray[3]
        )
        tagList += 1
        print(tagList)
        print(angle)
        return [angle, tagList]

    def startCalculation(path, sema):
        sema.acquire()
        root = os.getcwd()
        os.chdir(path)
        os.system("./rotation.sh")
        os.chdir(root)
        sema.release()

    if os.path.exists("RotationScreening"):
        shutil.rmtree("RotationScreening")
    os.mkdir("RotationScreening")
    os.chdir("RotationScreening")
    system.write("rotation.xyz")
    torsAngle, torsAtoms = getAngleAndAtoms(system)
    params["idef"] = "f tors {} {} {} {}".format(
        torsAtoms[0], torsAtoms[1], torsAtoms[2], torsAtoms[3]
    )
    if mode == "tm":
        import ase_own
    system.substituents[substituent].connection = system.substituents[
        substituent
    ].get_connectedAtom(system.core)
    for i in range(rotrange[0], rotrange[1], rotrange[2]):
        os.mkdir("{:03}".format(i))
        os.chdir("{:03}".format(i))
        for j in dftFiles:
            shutil.copyfile("../../{}".format(j), "{}".format(j))
        system.substituents[substituent].torsion(rotrange[2])
        for i in system:
            for k in system.substituents[substituent]:
                if i.tag == k.tag:
                    i.position = k.position
        if mode == "tm":
            system_2 = ase_own.Atoms(system)
            calculator = Turbomole(**params)
            system_2.set_calculator(calculator)
            calculator.initialize()
            shutil.copyfile("../../rotation.sh", "rotation.sh")
            os.system("chmod u+x rotation.sh")
        system.write("start.xyz")
        os.chdir("..")
        if writeXYZ:
            system.appendToXYZ("rotation.xyz")
    if mode == "tm":
        numberOfThreads = threads
        sema = Semaphore(numberOfThreads)
        counter = rotrange[0]
        for i in range(rotrange[0], rotrange[1], rotrange[2]):
            jobs = []
            p = multiprocessing.Process(
                target=startCalculation, args=("{:03}".format(counter), sema)
            )
            jobs.append(p)
            counter += rotrange[2]
            p.start()
        for p in jobs:
            p.join()

    os.chdir("..")


def analyseRotation(outFile="fsp.out", plot=True):
    """
    Analyses the rotational scan created by the *createRotationStructs* function. The
    final single-point calculation must be named *fsp.out* or *sp.out*.

    Parameters
    ----------

    outFile: str
        Name of the final single point calculation of the different screening-steps.
    plot: bool
        Whether or not a simple plot is beeing created.
    """
    folders = os.listdir("RotationScreening")
    folderList = []
    for i in folders:
        if cm.util.general.isInt(i):
            folderList.append(int(i))
    folderList.sort()
    energyArray = []
    optStructures = []
    for folder in folderList:
        with open("RotationScreening/{:03}/{}".format(folder, outFile), "r") as source:
            for line in source:
                if "total energy" in line and not "convergence" in line:
                    energyArray.append([folder, float(line.split()[4])])
        appendSystem = cm.System.Structure(
            read("RotationScreening/{:03}/coord".format(folder))
        )
        optStructures.append(appendSystem)
    for i in optStructures:
        i.appendToXYZ("RotationScreening/rotation_opt.xyz")
    energyArray = np.array(energyArray)
    energyMin = min(energyArray[:, 1])
    with open("RotationScreening/data.dat", "w") as target:
        target.write(
            "#{:>9}{:>20}{:>20}{:>20}\n".format(
                "Angle", "energy", "relative energy", "rel. energy / kJ/mol"
            )
        )
        for i in energyArray:
            target.write(
                "{:10.0f}{:20.8f}{:20.8}{:20.8f}\n".format(
                    i[0],
                    i[1],
                    i[1] - energyMin,
                    (i[1] - energyMin) * 2625.4996427573383,
                )
            )
    if plot:
        try:
            fig = plt.figure(figsize=(4, 3))
            energyArray[:, 1] -= energyMin
            energyArray[:, 1] *= 2625.4996427573383
            plt.plot(energyArray[:, 0], energyArray[:, 1])
            maxEnergy = max(energyArray[:, 1])
            print(
                "Rotation barrier for one substituent: {:.5f} kJ/mol".format(maxEnergy)
            )
            with open("RotationScreening/rotationBarrier", "w") as target:
                target.write("{}".format(max(energyArray[:, 1])))
            plt.tight_layout()
            plt.subplots_adjust(bottom=0.15, left=0.16)
            plt.xlabel("Rotation angle / °")
            plt.ylabel("relative Energy / kJ/mol")
            plt.savefig("RotationScreening/rotationBarrier.pdf")
        except:
            pass


def createEcxitation(
    params={}, exMode="rpas", exString="a 1", dftFiles=[], dimer=False
):
    """
    Creates a TDDFT-calculation from the optimized monomer structure. The caclulation will
    be carried out using the *excitation.sh* script in the base-directory.

    Parameters
    ----------
    params: dict
        Set of parameters, as beeing used by the *ase.calculator.Turbomole* package.
    exMode: str
        Kind of the excitation: *rpas* for singulett and *rpat* for triplett.
    exString: str
        Specify the kind of excitation, that should be calculated. Example: *a 20* calculates
        the first 20 excitation states of the a mode. I recommend to calculate without symmetry,
        so there are only *a* states and no *b* or higher states. Or specify the excitation modes
        by hand and ignore this function.
    dimer: bool
        If dimer is True, the orbitals of the dimer are being calculated. Otherwise, the monomer
        will be calculated.
    """
    trueRoot = os.getcwd()
    if dimer:
        folder = "dimer"
    else:
        folder = "monomer"
    os.mkdir("{}/opt/excitation".format(folder))
    os.chdir("{}/opt/excitation".format(folder))
    for j in dftFiles:
        shutil.copyfile("../../../{}".format(j), "{}".format(j))

    shutil.copyfile("../coord", "coord")
    system = read("coord", format="turbomole")
    params["excitation"] = [exMode, exString]

    calc = Turbomole(**params)
    system.set_calculator(calc)
    calc.initialize()
    shutil.copyfile("../../../excitation.sh", "excitation.sh")
    os.system("chmod u+x excitation.sh")
    os.system("./excitation.sh")
    os.chdir(trueRoot)


def analyseExcitation(outFile="exspectrum", plot=True):
    """
    Analyse the excitation calculation, which was executed in the *monomer/opt/excitation* folder.

    Parameters
    ----------

    outFile: str
        Name of the result of the excitation calculation. Default from Turbomole is *exspectrum*
    plot: bool
        Whether or not a simple plot is beeing created.
    """
    try:
        exSpectrum(
            "monomer/opt/excitation/{}".format(outFile),
            "monomer/opt/excitation/exPlot.pdf",
        )
    except:
        pass


def exportOrbitals(params={}, dftFiles=[], dimer=False):
    """
    Export orbitals as cube-files. Before this Function is beeing executed, the "more" function
    MUST be run before with the export of the MO-Orbital, as the "eiger.out"-File WILL be read-in
    in this function in order to write the cube-files.

    Parameters
    ----------

    params: dict
        Set of parameters, as beeing used by the *ase.calculator.Turbomole* package.
    dftFiles: list of str
        Additional files needed for the DFT calculation, e.g. pseudopotential.
    dimer: bool
        If dimer is True, the orbitals of the dimer are being calculated. Otherwise, the
        monomer will be calculated.
    """
    readOrbitals = False
    if dimer:
        folder = "dimer"
    else:
        folder = "monomer"
    with open("{}/opt/more/eiger.out".format(folder), "r") as source:
        for line in source:
            if readOrbitals:
                splitLine = line.split()
                if len(splitLine) == 8:
                    LUMO = int(splitLine[1])
                elif len(splitLine) == 9:
                    HOMO = int(splitLine[1])
                    break
            if "Orbital" in line:
                readOrbitals = True
    system = read("{}/opt/coord".format(folder), format="turbomole")
    if dimer:
        os.mkdir("Orbitals_dimer")
        os.chdir("Orbitals_dimer")
    else:
        os.mkdir("Orbitals")
        os.chdir("Orbitals")
    for j in dftFiles:
        shutil.copyfile("../{}".format(j), "{}".format(j))
    calc = Turbomole(**params)
    system.set_calculator(calc)
    calc.initialize()
    with open("control", "r") as source:
        with open("control_alt", "w") as target:
            for line in source:
                if not "$end" in line:
                    target.write(line)
                else:
                    target.write("$pointval mo {}-{} fmt=cub\n".format(HOMO, LUMO))
                    target.write("$end")
    shutil.move("control_alt", "control")
    shutil.copyfile("../sp.sh", "sp.sh")
    os.system("chmod u+x sp.sh")
    os.system("./sp.sh")
    shutil.move("{}a.cub".format(LUMO), "LUMO.cub")
    shutil.move("{}a.cub".format(HOMO), "HOMO.cub")
    os.chdir("..")


def getIE_EA(params={}, dftFiles=[]):
    """
    Calculate ionization energy and electron affinity. For relaxation, the *opt.sh* script
    will be used.

    Parameters
    ----------

    params: dict
        Set of parameters, as beeing used by the *ase.calculator.Turbomole* package.
    dftFiles: list of str
        Additional files needed for the DFT calculation, e.g. pseudopotential.
    """

    def getEnergy(file):
        with open(file, "r") as source:
            for line in source:
                if "|  total energy " in line:
                    return float(line.split()[4])

    trueRoot = os.getcwd()
    system = read("monomer/opt/coord", format="turbomole")
    os.mkdir("monomer/opt/charges")
    os.chdir("monomer/opt/charges")
    os.mkdir("IE")
    os.mkdir("EA")
    for j in dftFiles:
        shutil.copyfile("../../../{}".format(j), "IE/{}".format(j))
        shutil.copyfile("../../../{}".format(j), "EA/{}".format(j))
    shutil.copyfile("../../../opt.sh", "IE/opt.sh")
    shutil.copyfile("../../../opt.sh", "EA/opt.sh")
    os.chdir("IE")
    params["total charge"] = 1
    calc = Turbomole(**params)
    system.set_calculator(calc)
    calc.initialize()
    os.system("chmod u+x opt.sh")
    os.system("./opt.sh")
    os.chdir("../EA")
    params["total charge"] = -1
    calc = Turbomole(**params)
    system.set_calculator(calc)
    calc.initialize()
    os.system("chmod u+x opt.sh")
    os.system("./opt.sh")
    os.chdir(trueRoot)
    orgEnergy = getEnergy("monomer/opt/final_sp.out")
    IEenergy = getEnergy("monomer/opt/charges/IE/final_sp.out")
    EAenergy = getEnergy("monomer/opt/charges/EA/final_sp.out")
    with open("charges", "w") as target:
        target.write(
            "Ionization energy:   {} eV\n".format(
                (IEenergy - orgEnergy) * 27.2113966413
            )
        )
        target.write(
            "Electron Affinity:   {} eV".format((EAenergy - orgEnergy) * 27.2113966413)
        )


def sumAllInfo(
    folder,
    summaryFolder="ScreeningSummary",
    moreFiles=["eiger.out", "aoforce.out"],
    plotFiles=["monomer.pml", "dimer.pml"],
    plotStructs=False,
):
    """
    Method to analyse the screening of a chemical composition. This screening must be
    performed with the included screening functions.

    Parameters
    ----------

    folder: str
        Name of the folder with all the screening calculations.
    summaryFolder: str
        Folder, where the analysis will be saved to.
    moreFiles: list of str
        The screening includes a *moreInfo* function. This function will calculate some
        additional info, the user can choose. To store the important files from this function,
        write the name of these files into the *moreFiles*-list. This also functions as a list
        for additional analysis, namely the plotting of the MO-diagramm (if *eiger.out* is in
                the list) and the vibration spectrum (if *aoforce.out* is in the list.)
    plotFiles: list of str
        Extermal scripts for plotting the geometric structures. I use pymol, but others are
        possible as well, if scriptable.
    plotStructs: bool
        Whether or not the structures will be plottet with pymol. For the filtering on the
        cluster, I recommend to leave this on *False*, if working on a cluster
        for data-filtering and use a small script to plot it all afterwards.
    """
    rootDir = os.getcwd()
    if not os.path.exists(summaryFolder):
        os.mkdir(summaryFolder)
    if not os.path.exists("{}/{}".format(summaryFolder, folder)):
        os.mkdir("{}/{}".format(summaryFolder, folder))
    if os.path.exists("{}/Orbitals".format(folder)):
        monomer = read("{}/Orbitals/coord".format(folder))
    else:
        monomer = read("{}/monomer/opt/coord".format(folder))
    dimer = read("{}/dimer/opt/coord".format(folder))
    try:
        if not os.path.exists("{}/{}/monomer.xyz".format(summaryFolder, folder)):
            monomer.write("{}/{}/monomer.xyz".format(summaryFolder, folder))
        if not os.path.exists("{}/{}/dimer.xyz".format(summaryFolder, folder)):
            dimer.write("{}/{}/dimer.xyz".format(summaryFolder, folder))
    except:
        print("{} failed. (xyz-Files)".format(folder))
    os.chdir(rootDir)

    try:
        if not os.path.exists("{}/bindingEnergies/dataFile.dat".format(folder)):
            os.chdir(folder)
            bindingAnalysis()
            os.chdir("..")
            shutil.copyfile(
                "{}/bindingEnergies/bondShare.pdf".format(folder),
                "{}/{}/bondShare.pdf".format(summaryFolder, folder),
            )
            shutil.copyfile(
                "{}/bindingEnergies/dataFile.dat".format(folder),
                "{}/{}/bondData.dat".format(summaryFolder, folder),
            )
        elif not os.path.exists("{}/{}/bondShare.pdf".format(summaryFolder, folder)):
            os.chdir("{}/bindingEnergies".format(folder))
            cm.analyze.bonding.plotBondshare(width=5, height=3.2, legendloc=None)
            os.chdir("../..")
            shutil.copyfile(
                "{}/bindingEnergies/bondShare.pdf".format(folder),
                "{}/{}/bondShare.pdf".format(summaryFolder, folder),
            )
            shutil.copyfile(
                "{}/bindingEnergies/dataFile.dat".format(folder),
                "{}/{}/bondData.dat".format(summaryFolder, folder),
            )
    except:
        print("{} failed. (Binding energies)".format(folder))
    os.chdir(rootDir)

    try:
        for outFile in moreFiles:
            if not os.path.exists(f"{summaryFolder}/{folder}/monomer_{outFile}"):
                shutil.copyfile(
                    f"{folder}/monomer/opt/more/{outFile}",
                    f"{summaryFolder}/{folder}/monomer_{outFile}",
                )
            if not os.path.exists(f"{summaryFolder}/{folder}/dimer_{outFile}"):
                shutil.copyfile(
                    f"{folder}/dimer/opt/more/{outFile}",
                    f"{summaryFolder}/{folder}/dimer_{outFile}",
                )
        monomerGAP, dimerGAP = -9999, -9999
        monomerGAP = MO_Diagram(
            f"{summaryFolder}/{folder}/monomer_eiger.out",
            f"{summaryFolder}/{folder}/monomer_MO.pdf",
        )
        dimerGAP = MO_Diagram(
            f"{summaryFolder}/{folder}/dimer_eiger.out",
            f"{summaryFolder}/{folder}/dimer_MO.pdf",
        )
    except:
        print(f"{folder} failed. (HOMO- LUMO Gap)")
    os.chdir(rootDir)

    try:
        for outFiles in moreFiles:
            if outFiles == "eiger.out":
                if not os.path.exists(f"{summaryFolder}/{folder}/monomer_MO.pdf"):
                    monomerGAP = MO_Diagram(
                        f"{summaryFolder}/{folder}/monomer_eiger.out",
                        f"{summaryFolder}/{folder}/monomer_MO.pdf",
                    )
                    dimerGAP = MO_Diagram(
                        f"{summaryFolder}/{folder}/dimer_eiger.out",
                        f"{summaryFolder}/{folder}/dimer_MO.pdf",
                    )
            elif outFiles == "vibspectrum":
                if not os.path.exists(
                    f"{summaryFolder}/{folder}/monomer_Vibspectrum.pdf"
                ):
                    vibSpectrum(
                        f"{summaryFolder}/{folder}/monomer_vibspectrum",
                        f"{summaryFolder}/{folder}/monomer_Vibspectrum.pdf",
                    )
                if not os.path.exists(
                    f"{summaryFolder}/{folder}/dimer_Vibspectrum.pdf"
                ):
                    vibSpectrum(
                        f"{summaryFolder}/{folder}/dimer_vibspectrum",
                        f"{summaryFolder}/{folder}/dimer_Vibspectrum.pdf",
                    )
    except:
        print(f"{folder} failed. (Vibration spectrum)")
    os.chdir(rootDir)

    try:
        if os.path.exists(f"{folder}/monomer/opt/excitation"):
            if os.path.exists(f"{folder}/monomer/opt/excitation/eiger.out"):
                monomerGAP_tddft = MO_Diagram(
                    f"{folder}/monomer/opt/excitation/eiger.out",
                    f"{summaryFolder}/{folder}/monomer_MO_tddft.pdf",
                )
            if not os.path.exists(f"{folder}/monomer/opt/excitation.pdf"):
                os.chdir(folder)
                analyseExcitation()
                os.chdir("..")
            if not os.path.exists(f"{summaryFolder}/{folder}/monomer_excitation.pdf"):
                shutil.copyfile(
                    f"{folder}/monomer/opt/excitation/exPlot.pdf",
                    f"{summaryFolder}/{folder}/monomer_excitation.pdf",
                )
    except:
        print(f"{folder} failed. (Excitation spectrum)")
    os.chdir(rootDir)

    try:
        rotBarrier = np.nan
        if os.path.exists(f"{folder}/RotationScreening/rotationBarrier"):
            with open(f"{folder}/RotationScreening/rotationBarrier"):
                rotBarrier = float(target.readline())
        rot_exists = os.path.exists(f"{folder}/RotationScreening")
        rotsub_exists = os.path.exists(
            f"{summaryFolder}/{folder}/substituent_rotation.pdf"
        )
        if rot_exists and not rotsub_exists:
            if not os.path.exists(f"{folder}/RotationScreening/rotationBarrier.pdf"):
                os.chdir(folder)
                analyseRotation()
                os.chdir("..")
            shutil.copyfile(
                f"{folder}/RotationScreening/rotationBarrier.pdf"
                f"{summaryFolder}/{folder}/substituent_rotation.pdf"
            )
            shutil.copyfile(
                f"{folder}/RotationScreening/rotation.xyz"
                f"{summaryFolder}/{folder}/rotation_structuresSTART.xyz"
            )
            shutil.copyfile(
                f"{folder}/RotationScreening/rotation_opt.xyz"
                f"{summaryFolder}/{folder}/rotation_structuresOPT.xyz"
            )
        totBondEng = 0
    except:
        print(f"{folder} failed. (Rotiationbarrier)")
    os.chdir(rootDir)

    try:
        orb_exists = os.path.isdir(f"{folder}/Orbitals")
        cube_exists = os.path.exists(f"{summaryFolder}/{folder}/HOMO.cube")
        if orb_exists and not cube_exists:
            shutil.copyfile(
                f"{folder}/Orbitals/HOMO.cub", f"{summaryFolder}/{folder}/HOMO.cube"
            )
            shutil.copyfile(
                f"{folder}/Orbitals/LUMO.cub", f"{summaryFolder}/{folder}/LUMO.cube"
            )
    except:
        print(f"{folder} failed. (HOMO and LUMO cube files)")
    os.chdir(rootDir)

    try:
        with open(f"{summaryFolder}/{folder}/bondData.dat", "r") as source:
            for index, line in enumerate(source):
                if index == 1:
                    totBondEng = float(line.split()[2]) * 2625.4996427573383
                    coreBondEng = float(line.split()[3])
                    subBondEng = float(line.split()[4])
                elif index == 2:
                    totBondEngD = float(line.split()[2]) * 2625.4996427573383
                    coreBondEngD = float(line.split()[3])
                    subBondEngD = float(line.split()[4])
        if plotStructs:
            for i in plotFiles:
                shutil.copyfile(f"{i}", f"{summaryFolder}/{folder}/{i}")
                os.chdir(f"{summaryFolder}/{folder}")
                os.system(f"pymol {i}".format(i))
                os.chdir("../..")
    except:
        print(f"{folder} failed. (Binding data)")
    os.chdir(rootDir)

    try:
        IEenergy, EAenergy = 0, 0
        if os.path.exists(f"{folder}/charges"):
            with open(f"{folder}/charges", "r") as source:
                for line in source:
                    if "Ionization" in line:
                        IEenergy = float(line.split()[2])
                    if "Affinity" in line:
                        EAenergy = float(line.split()[2])
    except:
        print(f"{folder} failed. (charged calculations)")

    try:
        with open("{}/{}/data.csv".format(summaryFolder, folder), "w") as target:
            target.write("property,value,unit\n")
            target.write("Monomer HOMO LUMO Gap,{:.3f},eV\n".format(monomerGAP))
            target.write("Dimer HOMO LUMO Gap,{:.3f},eV\n".format(dimerGAP))
            target.write("Total bond energy,{:.3f},kJ/mol\n".format(totBondEng))
            target.write(
                "$E_{bond, core} / E_{bond, sub}$,{:.3f},\n".format(
                    coreBondEng / subBondEng
                )
            )
            target.write(
                "Total dispersive bond energy,{:.3f},kJ/mol\n".format(totBondEngD)
            )
            target.write(
                "$E_{disp, core} / E_{disp, sub}$,{:.3f},\n".format(
                    coreBondEngD / subBondEngD
                )
            )
            target.write(
                "Substituent rotation barrier,{:.3f},kJ/mol\n".format(rotBarrier)
            )
            target.write("Ionization energy,{:.5f},eV\n".format(IEenergy))
            target.write("Electron affinity,{:.5f},eV\n".format(EAenergy))
    except:
        print(f"{folder} failed. (data summary tabular)")
    os.chdir(rootDir)


def plotStructs(
    plotFiles=["monomer.pml", "dimer.pml"], summaryFolder="ScreeningSummary"
):
    """
    If the primary analysis of the screening was performed on a cluster without the possibility
    of on-the-fly rendering of high-quality pictures,
    This function allows for an easy creation of these pictures from the summary-folder.

    Parameters
    ----------

    plotFiles: list of str
        pymol-scripts for the different plots. These scripts will be executed with the pymol-order.
    summaryFolder: str
        Folder, where the analysis will be saved to.
    """
    for folder in os.listdir("ScreeningSummary"):
        if not os.path.exists(f"ScreeningSummary/{folder}/monomer.png"):
            for i in plotFiles:
                shutil.copyfile(f"{i}", f"{summaryFolder}/{folder}/{i}")
                os.chdir(f"{summaryFolder}/{folder}")
                os.system(f"pymol {i}")
                os.chdir("../..")


def plotOrbitals(
    plotFiles=["monomer_HOMO.pml", "monomer_LUMO.pml"],
    summaryFolder="ScreeningSummary",
    allNew=False,
):
    """
    If the primary analysis of the screening was performed on a cluster without the possibility
    of on-the-fly rendering of high-quality pictures,
    This function allows for an easy creation of these pictures from the summary-folder.

    Parameters
    ----------

    plotFiles: list of str
        pymol-scripts for the different plots. These scripts will be executed with the pymol-order.
    summaryFolder: str
        Folder, where the analysis will be saved to.
    allNew: bool
        Whether or not all the structures should be rendered new.
    """
    for folder in os.listdir("ScreeningSummary"):
        print(folder)
        if not os.path.exists(f"ScreeningSummary/{folder}/HOMO.png") or allNew:
            print(folder)
            for i in plotFiles:
                shutil.copyfile(f"{i}", f"{summaryFolder}/{folder}/{i}")
                os.chdir("{}/{}".format(summaryFolder, folder))
                if os.path.exists("HOMO.cub"):
                    os.system("mv HOMO.cub HOMO.cube")
                if os.path.exists("LUMO.cub"):
                    os.system("mv LUMO.cub LUMO.cube")
                os.system(f"pymol {i}")
                os.chdir("../..")


def exchange_opt(
    sourceFolder, substituteAtoms, params={}, dftFiles=[], optScript="", exchange="both"
):
    """
    Uses an already optimized structure with a different chemical composition as a basis
    for the optimization of the new sysystem.

    Parameters
    ----------

    sourceFolder: str
        This function will go BACK one instance (to the foldersystem with the many systems)
        and goes to the given folder, extracting the structure for the monomer and the dimer,
        exchanges the atoms and uses these as a basis for the optimization.
    substituteAtoms: tuple
        Exchanges the (first tuple-instance) atom with the (second instance-tuple) atom to
        match the wanted system.
    params: dict
        Set of parameters, as beeing used by the *ase.calculator.Turbomole* package.
    dftFiles: list of strings
        Files, which are needed for the optimization on dft level of theory. They will be copied
        into the *"opt"* folders within the crest-run. The slurm script has to be named
        "*job.slurm*", if the optimization will be started right away.
    optScript: str
        Simple slurm script, which will start the optimization.
    exchange: str
        "both", "monomer", or "dimer". Which will be taken from the other folder and optimized
        again. Usually there is no need to optimize the monomer as well.
    """
    if exchange == "both":
        folders = ["monomer", "dimer"]
    elif exchange == "dimer":
        folders = ["dimer"]
    elif exchange == "monomer":
        folders = ["monomer"]
    trueRoot = os.getcwd()
    for folder in folders:
        if not os.path.isdir(folder):
            os.mkdir(folder)
            os.mkdir(f"{folder}/opt")
        else:
            if not os.path.isdir(f"{folder}/opt"):
                os.mkdir(f"{folder}/opt")
    baseSystems = []
    if exchange in ("monomer", "both"):
        baseSystems.append(
            read(f"../{sourceFolder}/monomer/opt/coord", format="turbomole")
        )
    if exchange in ("dimer", "both"):
        baseSystems.append(
            read(f"../{sourceFolder}/dimer/opt/coord", format="turbomole")
        )

    for system in baseSystems:
        for atom in system:
            if atom.symbol == substituteAtoms[0]:
                atom.symbol = substituteAtoms[1]
    for i, folder in enumerate(folders):
        baseSystems[i].write(f"{folder}/opt/start.xyz")

    for folder in folders:
        os.chdir(folder)
        shutil.copyfile(f"../{optScript}", f"opt/{optScript}")
        for dftfile in dftFiles:
            shutil.copyfile(f"{trueRoot}/{dftfile}", f"opt/{dftfile}")
        os.chdir("opt")
        calc = Turbomole(**params)
        mol = read("start.xyz")
        mol.set_calculator(calc)
        calc.initialize()
        os.system(f"chmod u+x {optScript}")
        os.system(optScript)
        os.chdir(trueRoot)
