import os
import shutil
from copy import deepcopy

import ClusterMaestro as cm
import matplotlib.pyplot as plt
import numpy as np
from ClusterMaestro.initilize import init
from ase.io import read
from ase_own.calculators.turbomole import Turbomole

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


def initStructureParts(
    readInput=None,
    molecule="AdPh4",
    core=None,
    substituent=None,
    distance=6.0,
    subAngle=None,
    rotAngle=None,
    xshift=0.0,
    yshift=0.0,
    core_cutoff=25,
    threads=1,
    inputType="xyz",
    coreAtoms=[],
    analyzedDimer=None,
):
    """
    Initilizes a dimer with the built-in initilize-function, optimizes it on the highest xtb-level *(GFN2, acc=extreme)* and
    seperates the two molecules, the cores and substituents. All will be put into a different, new created folder. Other
    folders with the same name will be removed.

    Parameters
    ----------

    molecule: str
        Choose from fully pre-optimized molecules:
        AdPh4 (default)
    core: str
        Choose core structures from: See :ref:`Structure parts`
    substituent: str
        Choose substituents from: See :ref:`Structure parts`
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
    core_cutoff: int
        The N densest atoms to keep after the initial KDE. See find_cores for more info.
    inputType: str
        Type of the inputstructure. (see ASE)
    analyzedDimer: Dimer
    """
    for i in folderList:
        if os.path.exists(i):
            shutil.rmtree(i)
        os.mkdir(i)
    if analyzedDimer == None:
        if readInput == None:
            initDimer = init.dimer(
                molecule=molecule,
                core=core,
                substituent=substituent,
                subAngle=subAngle,
                distance=distance,
                rotAngle=rotAngle,
                xshift=xshift,
                yshift=yshift,
            )
            dimer = cm.System.Oligomer(
                cm.util.optimize.opt_xtb(initDimer, acc="normal", threads=threads)
            )
        else:
            dimer = cm.System.Oligomer(read(readInput))  # , format=inputType))
        dimer.write("dimer/start.xyz", format="xyz")
        if core == None or substituent == None:
            dimer.analyze(
                core=core,
                substituent=substituent,
                core_cutoff=core_cutoff,
                elements=coreAtoms,
            )
        else:
            dimer.matchPattern(core=core, substituent=substituent)
    else:
        dimer = analyzedDimer
        print(dimer)
        dimer.write("dimer/start.xyz", format="xyz")
    print(dimer)
    for i in dimer.molecules:
        print(i)
        print(i.core)
        for j in i.substituents:
            print(j)

    dimer.molecules[0].write("mono1/start.xyz", format="xyz")
    dimer.molecules[1].write("mono2/start.xyz", format="xyz")
    dimer.molecules[0].core.saturateCH()
    dimer.molecules[1].core.saturateCH()
    cores = cm.System.molecule.Molecule(deepcopy(dimer.molecules[0].core))
    cores += deepcopy(dimer.molecules[1].core)
    cores.saturateCH()
    cores.write("core/start.xyz", format="xyz")
    dimer.molecules[0].core.write("core1/start.xyz", format="xyz")
    dimer.molecules[1].core.write("core2/start.xyz", format="xyz")
    if len(dimer.molecules[0].substituents) == 0:
        for i in dimer.molecules:
            subs = i.find_substituents()

    def mergeSubs(subs):
        allSubs = cm.System.molecule.Molecule(subs[0])
        for i in range(len(subs)):
            if i != 0:
                allSubs += subs[i]
        return allSubs

    for i in dimer.molecules:
        for j in i.substituents:
            j.saturateCH()
    sub1 = mergeSubs(dimer.molecules[0].substituents)
    sub1.write("substituent1/start.xyz", format="xyz")
    sub2 = mergeSubs(dimer.molecules[1].substituents)
    sub2.write("substituent2/start.xyz", format="xyz")
    subs = deepcopy(sub1)
    subs += sub2
    subs.write("substituent/start.xyz", format="xyz")


def calcSinglePoint(threads=1):
    """
    Goes throught the folders using the folderList and performes a single-point calculation on xTB-level.

    Parameters
    ----------

    threads: int
        number of cpu-threads used for the xTB-calculation.
    """
    trueRoot = os.getcwd()
    for i in folderList:
        os.chdir(i)
        os.system("xtb start.xyz -P %s > xtb.out" % threads)
        os.chdir(trueRoot)


def analyseSinglePoint(
    calcType="xTB", writeFile=True, systemName="not_spec", archive=False
):
    """
    Goes throught the folders in folderList and collocts the total energy of the systems as well as the dispersion energy.
    These values will be calculated with each other to bet a good overview over the bondingenergies in the system.

    Parameters
    ----------

    calcType: str
        Currently only xTB and gaussian, more will follow.
    writeFile: Bool
        Whether or not the data is written into a datafile. This file is needed for the built-in plot-function.
    systemName: str
        If a dataFile is written, a name will be handy for identification. Must be ONE STRING, WITHOUT BLANKS.
    archive: Bool
        If true: go into archive and re-calculate the energy.
    """
    if not archive:
        if not os.path.exists("archive"):
            os.mkdir("archive")
        if os.path.exists("archive/%s" % systemName):
            shutil.rmtree("archive/%s" % systemName)
        os.mkdir("archive/%s" % systemName)

    def getED_xtb(path, copy=True, name=False):
        if not name:
            for i in os.listdir(path):
                if ".out" in i or ".log" in i:
                    dataFile = i
                    dataPath = "%s/%s" % (path, dataFile)
        else:
            dataPath = path
        with open(dataPath, "r") as source:
            for line in source:
                if "TOTAL ENERGY" in line:
                    energy = line.split()[3]
                if "dispersion" in line:
                    dispersion = line.split()[3]
        if copy:
            shutil.copyfile(path + "/xtb.out", "archive/%s/%s.out" % (systemName, path))
            shutil.copyfile(
                path + "/start.xyz", "archive/%s/%s.xyz" % (systemName, path)
            )
        return [float(energy), float(dispersion)]

    def getED_gau(path, copy=True, name=False):
        if not name:
            for i in os.listdir(path):
                if ".out" in i or ".log" in i:
                    dataFile = i
                    dataPath = "%s/%s" % (path, dataFile)
        else:
            dataPath = path
        with open(dataPath, "r") as source:
            for line in source:
                if "SCF Done" in line:
                    energy = line.split()[4]
                if " Dispersion energy" in line:
                    dispersion = line.split()[-2].replace("D", "E")
        if copy:
            shutil.copyfile(
                "%s/%s" % (path, dataFile), "archive/%s/%s.out" % (systemName, path)
            )
        return [float(energy), float(dispersion)]

    def getED_tm(path, copy=True, name=False):
        if not name:
            for i in os.listdir(path):

                if (".out" in i or ".log" in i) and not "ASE" in i:
                    dataFile = i
                    dataPath1 = "%s/%s" % (path, dataFile)
                    with open(dataPath1, "r") as testFile:
                        for line in testFile:
                            # if "RUNNING PROGRAM" in line:
                            dataPath = "%s/%s" % (path, dataFile)
                            break
        else:
            dataPath = path
        energy, dispersion = None, None
        with open(dataPath, "r") as source:
            for line in source:
                if "total energy" in line and not "scf" in line:
                    energy = line.split()[4]
                if (
                    "Energy contribution" in line
                    or "empirical dispersive energy correction" in line
                ):
                    dispersion = line.split()[-1]
        if copy:
            shutil.copyfile("%s" % (dataPath), "archive/%s/%s.out" % (systemName, path))
            shutil.copyfile(
                "{}/start.xyz".format(dataPath.split("/")[0]),
                "archive/%s/%s.xyz" % (systemName, path),
            )
        if type(energy) != str:
            return None
        else:
            return [float(energy), float(dispersion)]

    def getBondE(dimer, mono1, mono2):
        return dimer - mono1 - mono2

    energyArray = []
    if archive:
        for i in folderList:
            if calcType == "xTB":
                energyArray.append(
                    getED_xtb(
                        "archive/%s/%s.out" % (systemName, i), copy=False, name=True
                    )
                )
            if calcType == "gau":
                energyArray.append(
                    getED_gau(
                        "archive/%s/%s.out" % (systemName, i), copy=False, name=True
                    )
                )
            if calcType == "tm":
                energyArray.append(
                    getED_tm(
                        "archive/%s/%s.out" % (systemName, i), copy=False, name=True
                    )
                )
    else:
        for i in folderList:
            if calcType == "xTB":
                energyArray.append(getED_xtb(i))
            elif calcType == "gau":
                energyArray.append(getED_gau(i))
            elif calcType == "tm":
                energyArray.append(getED_tm(i))
    energyArray = np.array(energyArray)
    print("\n\n   ***  ANALYSIS OF THE BONDING IN THE GIVEN MOLECULE  ***\n")
    EBondList = []
    for i in range(0, 9, 3):
        print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
        print("         Bonding energy of: %s:" % folderList[i].upper())
        EBondList.append(
            [getBondE(energyArray[i, 0], energyArray[i + 1, 0], energyArray[i + 2, 0])]
        )
        print("Bondenergy (total): %.5f" % EBondList[-1][-1])
        EBondList[-1].append(
            getBondE(energyArray[i, 1], energyArray[i + 1, 1], energyArray[i + 2, 1])
        )
        print("Bondenergy (dispersion): %.5f" % EBondList[-1][-1])
    EBondList = np.array(EBondList)
    print("\nRelative parts of the TOTAL BONDING ENERGY:")
    print(
        "Relative part of the core-core-interaction (total): %.3f %%"
        % ((EBondList[1, 0] / EBondList[0, 0]) * 100)
    )
    print(
        "Relative part of the substituent-substituent-interaction (total): %.3f %%"
        % ((EBondList[2, 0] / EBondList[0, 0]) * 100)
    )
    print(
        "Relative part of the core-substituent-interaction (total): %.3f %%"
        % (
            (
                1
                - (EBondList[2, 0] / EBondList[0, 0])
                - (EBondList[1, 0] / EBondList[0, 0])
            )
            * 100
        )
    )
    print("\nRelative parts of the DISPERSIVE BONDING ENERGY:")
    print(
        "Relative part of the core-core-interaction (dispersion): %.3f %%"
        % ((EBondList[1, 1] / EBondList[0, 1]) * 100)
    )
    print(
        "Relative part of the substituent-substituentinteraction (dispersion): %.3f %%"
        % ((EBondList[2, 1] / EBondList[0, 1]) * 100)
    )
    print(
        "Relative part of the core-substituent-interaction (dispersion): %.3f %%"
        % (
            (
                1
                - (EBondList[2, 1] / EBondList[0, 1])
                - (EBondList[1, 1] / EBondList[0, 1])
            )
            * 100
        )
    )
    if writeFile:
        if not os.path.exists("dataFile.dat"):
            with open("dataFile.dat", "w") as target:
                target.write(
                    "#%19s%20s%20s%20s%20s\n"
                    % (
                        "system name",
                        "complete binding",
                        "core-core / %",
                        "sub-sub / %",
                        "core-sub / %",
                    )
                )
        with open("dataFile.dat", "a") as target:
            target.write(
                "%20s%20f%20f%20f%20f\n"
                % (
                    systemName + " total",
                    EBondList[0, 0],
                    EBondList[1, 0],
                    EBondList[2, 0],
                    (EBondList[0, 0] - EBondList[1, 0] - EBondList[2, 0]),
                )
            )
            target.write(
                "%20s%20f%20f%20f%20f\n"
                % (
                    systemName + " dispersion",
                    EBondList[0, 1],
                    EBondList[1, 1],
                    EBondList[2, 1],
                    (EBondList[0, 1] - EBondList[1, 1] - EBondList[2, 1]),
                )
            )


def plotBondshare(
    name="bondShare",
    part="both",
    width=7,
    height=5,
    dataFile="dataFile.dat",
    boxRange=None,
    totalBesides=False,
    legendloc="lower right",
    ylims=None,
):
    """
    Plots the bondshares, which were calculated and written in *dataFile.dat* as a bar-diagram. The Result may look something
    like the example in Figure 2.

    .. figure:: ../../doc/source/_static/bonding_plot.png


    Parameters
    ----------

    name: str
        Name of the resulting figure, WITHOUT the suffix.
    part: str
        Plot either the total binding energy (*total*), just the dispersive part (*dispersion*), or both next to each other.
    width: float
        Width of resulting image in inches.
    height: float
        Height of resulting image in inches.
    dataFile: str
        File containing data. Preferrebly produced by the *analyseSinglePoint*-function.
    boxRange: int
        Boxes are produced for a better optical seperation of the datapoints.
    totalBesides: bool
        Whether the total binding energy is beeing plotted as well or not.
    legendloc: str
        Where should the legend be? None removes the legend.
    ylims: tuple or None
        What should the limits for the y-axis be?
    """
    labelList = []
    totalArray = []
    dispArray = []
    with open(dataFile, "r") as source:
        for line in source:
            if "total" in line:
                splitLine = line.split()
                labelList.append(splitLine[0])
                totalArray.append(
                    [
                        float(splitLine[2]),
                        float(splitLine[3]),
                        float(splitLine[4]),
                        float(splitLine[5]),
                    ]
                )
            elif "dispersion" in line:
                splitLine = line.split()
                dispArray.append(
                    [
                        float(splitLine[2]),
                        float(splitLine[3]),
                        float(splitLine[4]),
                        float(splitLine[5]),
                    ]
                )
    totalArray, dispArray = np.array(totalArray), np.array(dispArray)
    totalArray = totalArray * 2625.4996427573383
    dispArray = dispArray * 2625.4996427573383
    for i in totalArray:
        print(i[0])
    fig, ax = plt.subplots(figsize=(width, height))
    plt.grid(axis="y", zorder=10)
    for array in [totalArray, dispArray]:
        for i in array[0]:
            array[0, 3] = 0 if array[0, 3] > 0 else array[0, 3]
    if part == "both":
        x = np.arange(len(labelList))
        width = 0.28 if totalBesides else 0.4
        offset = width / 2 if totalBesides else 0
        ax.bar(
            x - width / 2 + offset,
            totalArray[:, 2],
            color=cm.util.constants.colorList[0],
            label="substituent-substituent-interaction",
            width=width,
        )
        ax.bar(
            x - width / 2 + offset,
            totalArray[:, 1],
            bottom=totalArray[:, 2],
            color=cm.util.constants.colorList[1],
            label="core-core-interaction",
            width=width,
        )
        ax.bar(
            x - width / 2 + offset,
            totalArray[:, 3],
            bottom=totalArray[:, 2] + totalArray[:, 1],
            color="#f7a204",
            label="core-substituent-interaction",
            width=width,
        )
        ax.bar(
            x + width / 2 + offset,
            dispArray[:, 2],
            color=cm.util.constants.colorList[0],
            hatch="///",
            edgecolor="black",
            width=width,
        )
        ax.bar(
            x + width / 2 + offset,
            dispArray[:, 1],
            bottom=dispArray[:, 2],
            color=cm.util.constants.colorList[1],
            hatch="///",
            edgecolor="black",
            width=width,
        )
        ax.bar(
            x + width / 2 + offset,
            dispArray[:, 3],
            bottom=dispArray[:, 2] + dispArray[:, 1],
            color="#f7a204",
            hatch="///",
            edgecolor="black",
            width=width,
        )
        ax.bar(0, 0, color="white", hatch="///", edgecolor="black", label="dispersion")
        plt.ylabel("Binding energy / kJ/mol")
    if part == "total":
        width = 0.4 if totalBesides else 0.8
        offset = width / 2 if totalBesides else 0
        x = np.arange(len(labelList))
        ax.bar(
            x + offset,
            totalArray[:, 2],
            color=cm.util.constants.colorList[0],
            label="substituent$-$substituent interaction",
            width=width,
        )
        ax.bar(
            x + offset,
            totalArray[:, 1],
            bottom=totalArray[:, 2],
            color=cm.util.constants.colorList[1],
            label="core$-$core interaction",
            width=width,
        )
        ax.bar(
            x + offset,
            totalArray[:, 3],
            bottom=totalArray[:, 2] + totalArray[:, 1],
            color="#f7a204",
            label="core$-$substituent interaction",
            width=width,
        )
        plt.ylabel("Binding energy / kJ/mol")
    if part == "dispersion":
        width = 0.4 if totalBesides else 0.8
        offset = width / 2 if totalBesides else 0
        x = np.arange(len(labelList))
        ax.bar(
            x + offset,
            dispArray[:, 2],
            color=cm.util.constants.colorList[0],
            label="substituent-substituent-interaction",
            width=width,
        )
        ax.bar(
            x + offset,
            dispArray[:, 1],
            bottom=dispArray[:, 2],
            color=cm.util.constants.colorList[1],
            label="core-core-interaction",
            width=width,
        )
        ax.bar(
            x + offset,
            dispArray[:, 3],
            bottom=dispArray[:, 2] + dispArray[:, 1],
            color="#f7a204",
            label="core-substituent-interaction",
            width=width,
        )
        plt.ylabel("Dispersive part of binding energy / kJ/mol")
    if totalBesides:
        x = np.arange(len(labelList))
        width = 0.28 if part == "both" else 0.4
        offset = 0.2 if part != "both" else 0.28
        ax.bar(
            x - offset,
            totalArray[:, 0],
            color=cm.util.constants.colorList[3],
            label="total binding energy",
            width=width,
        )
    limits = plt.ylim()
    handles, labels = plt.gca().get_legend_handles_labels()
    if not legendloc == None:
        if totalBesides:
            order = [4, 0, 1, 2, 3]
            ax.legend(
                [handles[idx] for idx in order],
                [labels[idx] for idx in order],
                loc="lower right",
            )
        else:
            ax.legend(loc=legendloc)

    if boxRange != None:
        boxList = []
        for i in range(0, len(labelList), boxRange * 2):
            boxList.append((i - 0.5, boxRange))
        ax.broken_barh(
            boxList,
            (limits[0] - 100, limits[1] + 1000),
            facecolors="#c0f1f8",
            zorder=0.1,
        )
    ax.set_axisbelow(True)
    plt.xticks(range(len(labelList)), labelList, rotation=45)
    if ylims == None:
        plt.ylim(limits[0] + 0.1 * limits[0], limits[1])
    else:
        plt.ylim(ylims)
    plt.xlim(-0.5, len(labelList) - 0.5)
    plt.subplots_adjust(right=0.98, top=0.98, left=0.08, bottom=0.20)
    plt.tight_layout()
    plt.savefig(name + ".pdf", format="pdf")
    plt.savefig(name + ".png", format="png", dpi=300)


def fullAnalysis(
    readInput=None,
    molecule="AdPh4",
    core=None,
    substituent=None,
    distance=6.0,
    subAngle=None,
    rotAngle=None,
    xshift=0.0,
    yshift=0.0,
    core_cutoff=25,
    threads=1,
    calcType="xTB",
    writeFile=True,
    systemName="not_spec",
):
    """
    Performes a full creation, calculation and analysis *(bonding.initStructureParts, bonding.calcSinglePoint,
    bonding.analyseSinglePoint)* of a given system. One line of code to get the full bonding-energy (on an xTB-GFN2 level).

    For a description of the parameters, see the other functions in this module package.
    """
    initStructureParts(
        readInput,
        molecule,
        core,
        substituent,
        distance,
        subAngle,
        rotAngle,
        xshift,
        yshift,
        core_cutoff,
        threads,
    )
    calcSinglePoint(threads)
    analyseSinglePoint(calcType, writeFile, systemName)


def getBindingEnergy(
    dimer,
    dftFiles=[],
    params={},
    calculate=True,
    calculateDimer=False,
    decompose=False,
    spScript="",
):
    """
    Seperates the dimer into two molecules, initializes the calculation and stores everything in folders.

    Parameters
    ----------
    dftFiles: list of strings
        Files, which are needed for the optimization on dft level of theory. They will be copied into the
        *"opt"* folders within the crest-run. The slurm script has to be named "*job.slurm*", if the
        optimization will be started right away.
    params: dict
        Set of parameters, as beeing used by the *ase.calculator.Turbomole* package.
    calculate: bool
        Whether or not the calculation will just be started.
    calculateDimer: bool
        Whether or not the dimer will be calculated again.
    decompose: bool
        Whether or not the decomposed binding energy will be calculated.
    spScript: str
        Startscript for the singlepoint calculations. Files will be copied to each folder for the execution.
        If none is given, *ridft* will be executed.
    """
    if not decompose:
        for i in range(len(dimer.molecules)):
            os.mkdir("mol_{}".format(i))
            os.chdir("mol_{}".format(i))
            dimer.molecules[i].write("start.xyz")
            for dftfile in dftFiles:
                shutil.copyfile("../{}".format(dftfile), "{}".format(dftfile))
            mol = read("start.xyz")
            if calculate:
                calc = Turbomole(**params)
                mol.set_calculator(calc)
                calc.initialize()
                os.system("ridft > final_sp.out")
            os.chdir("..")
        if calculateDimer:
            dimer.write("start.xyz")
            calc = Turbomole(**params)
            mol = read("start.xyz")
            mol.set_calculator(calc)
            calc.initialize()
            os.system("ridft > final_sp.out")
        with open("final_sp.out", "r") as source:
            for line in source:
                if "|  total energy" in line:
                    dimerEnergy = float(line.split()[4])
        with open("mol_0/final_sp.out", "r") as source:
            for line in source:
                if "|  total energy" in line:
                    mol0Energy = float(line.split()[4])
        with open("mol_1/final_sp.out", "r") as source:
            for line in source:
                if "|  total energy" in line:
                    mol1Energy = float(line.split()[4])
        eBind = dimerEnergy - mol0Energy - mol1Energy
        with open("bindingEnergy", "w") as target:
            target.write("Dimer Energy: {}  Ha\n".format(dimerEnergy))
            target.write("monomer0 Energy: {}  Ha\n".format(mol0Energy))
            target.write("monomer1 Energy: {}  Ha\n\n".format(mol1Energy))
            target.write("Binding Energy: {} Ha\n".format(eBind))
            target.write("Binding Energy: {} kJ/mol".format(eBind * 2625.499))
    else:
        trueRoot = os.getcwd()
        initStructureParts(analyzedDimer=dimer)
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
                    "%s/%s" % (trueRoot, slurmfile), "%s/%s" % (i, slurmfile)
                )
        for i in folderList:
            os.chdir("%s" % i)
            if spScript != "":
                shutil.copyfile("../../%s" % spScript, spScript)
            mol = cm.System.Structure(read("start.xyz"))
            elemList = list(set(mol.get_chemical_symbols()))
            for i in range(len(elemList)):
                elemList[i] = elemList[i].lower()
            deleteECP = True
            params_2 = deepcopy(params)
            for i in params_2:
                if "ecp" in i:
                    if params_2[i].split()[0] in elemList:
                        deleteECP = False
            if deleteECP:
                try:
                    params_2.pop("ecp name")
                    params_2.pop("ecp file")
                    params_2.pop("basis set atom")
                    params_2.pop("basis set file")
                except:
                    None
            try:
                calc = Turbomole(**params_2)
                mol.set_calculator(calc)
                calc.initialize()
            except:
                params_2["use resolution of identity"] = False
                calc = Turbomole(**params_2)
                mol.set_calculator(calc)
                calc.initialize()
            if spScript != "":
                os.system("chmod u+x %s" % spScript)
                os.system("./%s" % spScript)
            else:
                os.system("ridft > sp.out")
            os.chdir("..")
        os.chdir(trueRoot)
        print(trueRoot)
