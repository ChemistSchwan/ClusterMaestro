import ClusterMaestro as cm
from ase.io import read
from statistics import mean
import numpy as np
import os, shutil, linecache, random
import pymatgen as pmg
import matplotlib.pyplot as plt


def splitter(dataFile="crest_conformers.xyz"):
    """
    Splits all structures contained in the "crest_conformers" file into seperate files.

    Parameters
    ----------

    dataFile: str
        File, which contains the crest-structures.
    """
    files = []
    firstLine = ""
    xyzString = ""
    if os.path.exists("allStructures"):
        shutil.rmtree("allStructures")
    os.mkdir("allStructures")

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
        files.append(xyzString)
    for i, string in enumerate(files):
        with open("allStructures/%.3d.xyz" % i, "w") as target:
            target.write(string)
    return len(files)


def sortCRESTrmsd(
    maxRMSD=1.0,
    dataFile="crest_conformers.xyz",
    core=False,
    core_cutoff=25,
    split=False,
):
    """
    Sorts the structures resulting from crest into groups depending on their RMSD.

    Parameters
    ----------
    maxRMSD: float
        RMSD value to decide, if the structures are "similar" ot not.
    dataFile: str
        File, which contains the crest-structures.
    core: bool
        True: only the core-structures are compared, False: the whole structures are comared.
    core_cutoff: int
         The N densest atoms to keep after the initial KDE. See find_cores for more info.
    split: bool
        Whether or not the 'crest_conformers.xyz' will be splitted. If False: takes the structures in 'allStructures'
    """
    molecules, rmsds, sortedMols, groups = [], [], [], []
    if split:
        number = splitter(dataFile=dataFile)
    else:
        number = len(os.listdir("allStructures"))
    print("There are %d different molecules to be sorted" % number)

    if os.path.exists("rmsdGroups"):
        shutil.rmtree("rmsdGroups")
    os.mkdir("rmsdGroups")

    fileList = os.listdir("allStructures")
    for i in fileList:
        molecules.append(cm.System.Oligomer(read("allStructures/%s" % i)))
        molecules[-1].energy = float(linecache.getline("allStructures/%s" % i, 2))
    groupCounter = -1
    for i, mol1 in enumerate(molecules):
        if i not in sortedMols:
            groupCounter += 1
            if core:
                mol1.analyze(core_cutoff)
            os.mkdir("rmsdGroups/%.3d" % groupCounter)
            shutil.copyfile(
                "allStructures/%s" % fileList[i],
                "rmsdGroups/%.3d/%s" % (groupCounter, fileList[i]),
            )
            sortedMols.append(i)
            for j, mol2 in enumerate(molecules):
                if j not in sortedMols:
                    if core:
                        mol2.analyze(core_cutoff)
                        cores1 = mol1.molecules[0].core + mol1.molecules[1].core
                        cores2 = mol2.molecules[0].core + mol2.molecules[1].core
                        rmsd = cm.util.rmsdFAST.combineAllRMSD(cores1, cores1)
                    else:
                        rmsd = cm.util.rmsdFAST.combineAllRMSD(mol1, mol2)
                    if rmsd < maxRMSD:
                        shutil.copyfile(
                            "allStructures/%s" % fileList[j],
                            "rmsdGroups/%.3d/%s" % (groupCounter, fileList[j]),
                        )
                        sortedMols.append(j)
    print("     Molecules successfully sorted for RMSD")


def sortCRESTsym(
    dataFile="crest_conformers.xyz",
    tolerance=0.3,
    eigen_tolerance=0.01,
    matrix_tol=0.1,
    split=False,
):
    """
    Takes a crest-output and sorts the resulting structures according to their symmetry into different folders.
    This function is using the pymatgen.symmetry.analyzer module.

    Explanation of the resulting symmetry: https://en.wikipedia.org/wiki/Point_group

    Parameters
    ----------

    dataFile: str
        File, which contains the crest-structures.
    tolerance: float
        Distance tolerance to consider sites as symmetrically equivalent.
    eigen_tolerance: float
        Tolerance to compare eigen values of the inertia tensor.
    matrix_tol: float
        Tolerance used to generate the full set of symmetry operations of the point group.
    split: bool
        Whether or not the 'crest_conformers.xyz' will be splitted. If False: takes the structures in 'allStructures'
    """
    molecules, rmsds, detectedSyms, groups = [], [], [], []
    if split:
        number = splitter(dataFile=dataFile)
    else:
        number = len(os.listdir("allStructures"))
    print("There are %d different molecules to be sorted" % number)

    fileList = os.listdir("allStructures")
    for i in fileList:
        molecules.append(cm.System.molecule.Molecule(read("allStructures/%s" % i)))
        try:
            molecules[-1].energy = float(linecache.getline("allStructures/%s" % i, 2))
        except:
            molecules[-1].energy = np.nan
    #        os.remove("crest_split_%.3d.xyz" % i)

    if os.path.exists("symmetryGroups"):
        shutil.rmtree("symmetryGroups")
    os.mkdir("symmetryGroups")
    for i, mol in enumerate(molecules):
        molecule = pmg.core.structure.Molecule(
            species=mol.symbols, coords=mol.positions
        )
        symmetry = pmg.symmetry.analyzer.PointGroupAnalyzer(
            molecule, tolerance, eigen_tolerance, matrix_tol
        ).sch_symbol
        if symmetry == "S2":
            symmetry = "Ci"
        print("Molecule %.3d has symmetry %s" % (i, symmetry))
        if symmetry not in detectedSyms:
            detectedSyms.append(symmetry)
            os.mkdir("symmetryGroups/%s" % symmetry)
        mol.write("symmetryGroups/%s/%s_org" % (symmetry, fileList[i]), format="xyz")
        with open("symmetryGroups/%s/%s_org" % (symmetry, fileList[i]), "r") as source:
            with open("symmetryGroups/%s/%s" % (symmetry, fileList[i]), "w") as target:
                for j, line in enumerate(source):
                    if j == 1:
                        line = "   %f\n" % mol.energy
                        target.write(line)
                    else:
                        target.write(line)
        os.remove("symmetryGroups/%s/%s_org" % (symmetry, fileList[i]))
    print("     Molecules successfully sorted for symmetry")


def plotSymEng(name="symEnergy", width=7, height=5):
    """
    Plots the energy (according to the value in the crest-File) vs the point-group. Run *sortCRESTsym* first!

    Parameters
    ----------
    name: str
        Name of the resulting figure, WITHOUT the suffix.
    width: float
       Width of resulting image in inches.
    height: float
        Height of resulting image in inches.
    """
    symList = os.listdir("symmetryGroups")
    energyList = []
    firstVal = None

    for i, sym in enumerate(symList):
        energyList.append([])
        for j, filename in enumerate(os.listdir("symmetryGroups/%s" % sym)):
            # if firstVal == None:
            #    firstVal = float(linecache.getline("symmetryGroups/%s/%s" % (sym, filename), 2))
            energyList[-1].append(
                float(linecache.getline("symmetryGroups/%s/%s" % (sym, filename), 2))
            )

    minVal = min(energyList)[0]
    for i in range(len(energyList)):
        for j in range(len(energyList[i])):
            energyList[i][j] = (float(energyList[i][j]) - minVal) * 2625.4996427573383

    fig, ax = plt.subplots(figsize=(width, height))
    for i, eList in enumerate(energyList):
        xarray = np.zeros(len(eList))
        xarray += i
        ax.scatter(xarray, eList, color="black")

    plt.xlabel("Point group")
    plt.ylabel("Relative energy / kJ/mol")
    plt.xticks(range(len(symList)), symList)

    plt.savefig(name + ".pdf", format="pdf")


def dimerBondEng(
    core=None,
    substituent=None,
    match=True,
    data="rmsd",
    name="BondEng",
    monomer="monomer.xyz",
    core_cutoff=25,
):
    """
    Takes the outcome of either the *sortCRESTrmsd* or *sortCRESTsym* and sorts them accordingly to their groups.
    Also plots them as I did in the paper from Hanau *et al.*

    Parameters
    ----------

    core: str
        Core structure. See :ref:`Structure parts`
    substituent: str
        Substituent structure. See :ref:`Structure parts`
    data: str
        Either *rmsd* or *sym*. Sorts data according to eigther the rmsd or the symmetry of the structures.
    name: str
        Filename of resulting plot WITHOUT SUFFIX!
    monomer: str
        XYZ-file containing the monomer-minimum structure.
    core_cutoff: int
        See *ClusterMaestro.System.dimer.Dimer.analyze()*
    """

    monomerEnergy = float(linecache.getline(monomer, 2))
    if data == "rmsd":
        rootFolder = "rmsdGroups"
    elif data == "sym":
        rootFolder = "symmetryGroups"
    folders = os.listdir(rootFolder)
    # if os.path.exists("data_structures"):
    #    shutil.rmtree("data_structures")
    # os.mkdir("data_structures")

    moleculeList = []
    for i in folders:
        subfiles = os.listdir("%s/%s" % (rootFolder, i))
        for j in range(len(subfiles)):
            subfiles[j] = int(subfiles[j][:-4])  # .replace(".xyz", "")
        moleculeList.append(
            [
                i,
                cm.System.Oligomer(
                    read(
                        "%s/%s/%.3d.xyz" % (rootFolder, i, min(subfiles)), format="xyz"
                    )
                ),
            ]
        )
        #   shutil.copyfile("%s/%s/%.3d.xyz" % (rootFolder, i, min(subfiles)), "data_structures/%s.xyz" % moleculeList[-1][0])
        moleculeList[-1][1].energy = float(
            linecache.getline("%s/%s/%.3d.xyz" % (rootFolder, i, min(subfiles)), 2)
        )

    monomerEnergy = float(linecache.getline(monomer, 2))

    dataArray, labelList = [], []
    with open("data.dat", "w") as target:
        target.write(
            "#%9s%20s%20s%20s\n" % ("molecule", "Energy", "Distance", "Diss. Energy")
        )
        for i, mol in enumerate(moleculeList):
            mol[1].analyze(core=core, substituent=substituent, core_cutoff=core_cutoff)
            target.write(
                "%10s%20.5f%20.5f%20.5f\n"
                % (
                    str(mol[0]),
                    mol[1].energy,
                    mol[1].get_core_distance(),
                    mol[1].energy - 2 * monomerEnergy,
                )
            )
            labelList.append(str(mol[0]))
            dataArray.append(
                [
                    mol[1].energy,
                    mol[1].get_core_distance(),
                    mol[1].energy - 2 * monomerEnergy,
                ]
            )
    dataArray = np.array(dataArray)
    dataArray[:, 0] *= 2625.5002
    dataArray[:, 2] *= 2625.5002

    fig, ax = plt.subplots()
    ax.scatter(range(len(labelList)), dataArray[:, 1])
    plt.xticks(range(len(labelList)), labelList)
    plt.ylabel("Core-Core distance / $\mathrm{\AA}$")
    plt.savefig(
        "%sBondEng.pdf" % rootFolder.replace("Groups", ""), dpi=500, transparent=True
    )

    plt.clf()
    fig, ax = plt.subplots()
    ax.scatter(dataArray[:, 1], dataArray[:, 2])
    plt.xlabel("Core-Core distance / $\mathrm{\AA}$")
    plt.ylabel("Dissociation Energy / kJ/mol")
    plt.savefig(
        "%sBondEDiss.pdf" % rootFolder.replace("Groups", ""), dpi=500, transparent=True
    )


def dimerSymData(core_cutoff=25, distanceLimit=None):
    """
    Collocts data about the symmetry, energy and core-core distance about all systems and saved them in a file, so it can
    be plotted later.

    Can also handle Oligomers. In this case, the smallest core-core distance will be written in

    Parameters
    ----------

    core_cutoff: int
        See *ClusterMaestro.System.dimer.Dimer.analyze()*
    distanceLimit: None or float
        See *ClusterMaestro.System.dimer.Dimer.get_core_distance()*
    """
    #    monomerEnergy = float(linecache.getline(monomer, 2))
    rootFolder = "symmetryGroups"
    folders = os.listdir(rootFolder)
    #    if os.path.exists("data_structures"):
    #        shutil.rmtree("data_structures")
    #    os.mkdir("data_structures")
    moleculeList, symList, numberList = [], [], []
    for i in folders:
        subfiles = os.listdir("%s/%s" % (rootFolder, i))
        moleculeList.append([])
        symList.append(i)
        for j in subfiles:
            numberList.append(j[:-4])
            moleculeList[-1].append(
                cm.System.Oligomer(read("%s/%s/%s" % (rootFolder, i, j), format="xyz"))
            )
            moleculeList[-1][-1].energy = float(
                linecache.getline("%s/%s/%s" % (rootFolder, i, j), 2)
            )
    #    monomerEnergy = float(linecache.getline(monomer, 2))
    dataArray, labelList = [], []
    counter, minimum = 0, 0
    with open("data.dat", "w") as target:
        target.write(
            "#%9s%10s%20s%20s%20s\n"
            % ("molecule", "symmetry", "Energy", "min. Distance", "avg. Distance")
        )
        for i, molList in enumerate(moleculeList):
            dataArray.append([])
            for j, mol in enumerate(molList):
                mol.analyze(core_cutoff)
                if len(mol.molecules) >= 2:
                    distance = mol.get_core_distance(limit=distanceLimit)
                else:
                    distance = np.nan
                try:
                    target.write(
                        "%10s%10s%20.5f%20.5f%20.5f\n"
                        % (
                            numberList[counter],
                            symList[i],
                            mol.energy,
                            distance,
                            distance,
                        )
                    )
                except:
                    target.write(
                        "%10s%10s%20.5f%20.5f%20.5f\n"
                        % (
                            numberList[counter],
                            symList[i],
                            mol.energy,
                            min(distance),
                            mean(distance),
                        )
                    )
                if mol.energy < minimum:
                    minimum = mol.energy
                counter += 1


def plotDimerSymData(name="dimerAnalysis", xRange=None, width=7, height=5, legend=True):
    """
    Plots the data generated from *dimerSymData*

    Parameters
    ----------

    name: str
        Filename of resulting plot WITHOUT SUFFIX!
    xRange: (float, float)
        Range of x-axis. If nothing is given, will default into auto range.
    width: float
       Width of resulting image in inches.
    height: float
        Height of resulting image in inches.
    legend: bool
        Whether or not a legend is shown.
    """
    symmetryList = [
        "C1",
        "C2",
        "C3",
        "C4",
        "Td",
        "D2",
        "D3",
        "D2h",
        "D2d",
        "D3d",
        "Cs",
        "Ci",
        "C*v",
        "C2v",
        "C3v",
        "S4",
    ]
    symColorList = [
        "#005b9a",
        "#15b434",
        "#c33232",
        "#972b9a",
        "#46c1c0",
        "#c07f2b",
        "#007ea0",
        "#941651",
        "#ffb000",
        "#696966",
        "#60ff05",
        "#ff5c1c",
        "#c802f9",
        "#0f6c75",
        "#fcf80c",
        "#fc0c78",
    ]
    dataArray, symList, minimum = [], [], 0
    with open("data.dat", "r") as source:
        for i, line in enumerate(source):
            if i == 1:
                splitLine = line.split()
                symList.append(splitLine[1])
                dataArray.append([])
                dataArray[-1].append(
                    [float(splitLine[2]) * 2625.5002, float(splitLine[3])]
                )
            elif i > 1:
                splitLine = line.split()
                if splitLine[1] != symList[-1]:
                    dataArray.append([])
                    symList.append(splitLine[1])
                dataArray[-1].append(
                    [float(splitLine[2]) * 2625.5002, float(splitLine[3])]
                )
            if i >= 1:
                if float(splitLine[2]) * 2625.5002 < minimum:
                    minimum = float(splitLine[2]) * 2625.5002
    fig, ax = plt.subplots(figsize=(width, height))
    for i, dat in enumerate(dataArray):
        dat = np.array(dat)
        ax.scatter(
            dat[:, 1],
            dat[:, 0] - minimum,
            color=symColorList[symmetryList.index(symList[i])],
            label=symList[i],
        )
    if not xRange == None:
        plt.xlim(xRange)

    if legend:
        plt.legend(loc="lower right")
    plt.ylabel("relative energy / kJ/mol")
    plt.xlabel("Core-Core distance / $\mathrm{\AA}$")
    plt.subplots_adjust(right=0.98, top=0.98, left=0.15, bottom=0.20)
    plt.savefig(name + ".pdf", dpi=500, transparent=True)


def plotDimerSymDataLegend(symList):
    """
    Plots just the legend of the *plotDimerSymData*- function.

    Parameters
    ----------

    symList: list
        List with all symmetry strings, that should be in the legend.
    """
    symmetryList = [
        "C1",
        "C2",
        "C3",
        "C4",
        "Td",
        "D2",
        "D3",
        "D2h",
        "D3d",
        "Cs",
        "Ci",
        "C*v",
        "C2v",
        "C3v",
        "S4",
    ]
    symColorList = [
        "#005b9a",
        "#15b434",
        "#c33232",
        "#972b9a",
        "#46c1c0",
        "#c07f2b",
        "#007ea0",
        "#941651",
        "#ffb000",
        "#696966",
        "#60ff05",
        "#ff5c1c",
        "#c802f9",
        "#0f6c75",
        "#fcf80c",
        "#fc0c78",
    ]
    for i in symList:
        plt.scatter([0], [0], label=i, color=symColorList[symmetryList.index(i)])
    plt.legend(ncol=int(len(symList) / 2))
    plt.savefig("legend.pdf", dpi=500)


def getSampleStructures(number=20):
    """
    Collects a number of structures from a pre-sorted CREST run.
    This CREST run must be sorted using the *sortCRESTrmsd* and *sortCRESTsym* functions od this module.

    The samples will be collected by the following rules:

    The best structure from each symmetry group is taken first as well as the best structure from each rmsd group.

    If this is not enough to fulfill the required number, the rest will be selected randomly.

    The structures are then copied to a new folder *sampleStructures* in a seperate folder each, for easier calculation afterwards.

    Parameters
    ----------

    number: int
        Number of molecules, which are chosen for the sample set. MUST be smaller or equal to the total set of structures.
    """
    structureNumbers = []
    if os.path.exists("sampleStructures"):
        shutil.rmtree("sampleStructures")
    os.mkdir("sampleStructures")
    symGroups = os.listdir("symmetryGroups")
    rmsdGroups = os.listdir("rmsdGroups")
    molSymGroups = []
    for i in symGroups:
        molSymGroups.append([])
        for j in os.listdir("symmetryGroups/%s" % i):
            molSymGroups[-1].append(int(j[:-4]))
    for i in molSymGroups:
        i.sort()
        if not len(structureNumbers) >= number:
            structureNumbers.append(i[0])
    molrmsdGroups = []
    for i in rmsdGroups:
        molrmsdGroups.append([])
        for j in os.listdir("rmsdGroups/%s" % i):
            molrmsdGroups[-1].append(int(j[:-4]))
    for i in molrmsdGroups:
        i.sort()
        if i[0] not in structureNumbers and not len(structureNumbers) >= number:
            structureNumbers.append(i[0])
    allList = os.listdir("allStructures")
    for i in range(len(allList)):
        allList[i] = int(allList[i][:-4])
    if len(allList) < number:
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("+ The total number of molecules is smaller than the number variable. +")
        print("+       Sampleset will be smaller and contain every molecule.        +")
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    if len(allList) >= number:
        while len(structureNumbers) < number:
            randomNumber = random.randint(0, len(allList) - 1)
            if randomNumber not in structureNumbers:
                structureNumbers.append(randomNumber)
    else:
        for i in range(len(allList)):
            if i not in structureNumbers:
                structureNumbers.append(i)
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("The chosen structures for the sample group are: ")
    print(structureNumbers)
    for i in structureNumbers:
        os.mkdir("sampleStructures/%.3d" % i)
        shutil.copyfile(
            "allStructures/%.3d.xyz" % i, "sampleStructures/%.3d/start.xyz" % i
        )


def getCoreBondLengths(core_cutoff=25):
    """
    Returns a distribution of the bondlengths within the core-structure of the *allStructures* folder.

    Parameters
    ----------

    core_cutoff: int
        See dimer.analyze.
    """
    rootFolder = "allStructures"
    files = os.listdir(rootFolder)
    averageList, maxList, minList, totalList = [], [], [], []
    for filename in files:
        dimer = cm.System.Oligomer(read("%s/%s" % (rootFolder, filename)))
        try:
            dimer.analyze(core_cutoff)
            for mol in dimer.molecules:
                bondList = mol.core.getCoreBondlength()
                averageList.append(bondList.sum() / len(bondList))
                maxList.append(bondList.max())
                minList.append(bondList.min())
                for i in bondList:
                    totalList.append(i)
        except:
            print(
                "Could not resolve structure. Structure %s will not be considered."
                % filename
            )
    totalList = np.array(totalList)
    return totalList
