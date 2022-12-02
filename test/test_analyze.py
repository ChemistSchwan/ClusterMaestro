import ClusterMaestro as cm
import os, shutil


def test_bonding():
    trueRoot = os.getcwd()
    os.chdir("cases/bondingTest")
    print(os.getcwd())
    os.mkdir("initTest")
    os.chdir("initTest")
    cm.analyze.bonding.initStructureParts(readInput="../../SiSPh4_dimer.xyz", core="SiS", substituent="Ph")
    os.chdir("..")
    cm.analyze.bonding.analyseSinglePoint(calcType="tm", systemName="SiPh4", archive=True)
    cm.analyze.bonding.analyseSinglePoint(calcType="gau", systemName="AdMe4", archive=True)
    cm.analyze.bonding.analyseSinglePoint(systemName="AdPh", archive=True)
    cm.analyze.bonding.plotBondshare()
    shutil.rmtree("initTest")
    for file in ["bondShare.pdf", "bondShare.png", "dataFile.dat"]: os.remove(file)
    os.chdir(trueRoot)


def test_crestutils():
    trueRoot = os.getcwd()
    os.chdir("cases/crest")
    cm.analyze.crestutils.splitter()
    cm.analyze.crestutils.sortCRESTrmsd()
    cm.analyze.crestutils.dimerBondEng(core="Ad", substituent="Me")
    cm.analyze.crestutils.getCoreBondLengths()
    for folder in ["allStructures", "rmsdGroups"]: shutil.rmtree(folder)
    for file in ["data.dat", "rmsdBondEDiss.pdf", "rmsdBondEng.pdf"]: os.remove(file)
    os.chdir(trueRoot)

def test_hyperplane():
    trueRoot = os.getcwd()
    os.chdir("cases")
    os.mkdir("hyperplane")
    os.chdir("hyperplane")
    cm.analyze.hyperplane.create_structures(core="Ad", substituent="Me", distRange=[4, 6.01, 1], rotRange=[0, 60, 30], subAngle=-20)#, dryRun=True)
    cm.analyze.hyperplane.calc_hyperplane(threads=4, acc="crude")
    cm.analyze.hyperplane.create_array()
    cm.analyze.hyperplane.plot(label="contour3D.png", mode="3D", zlim=200)
    cm.analyze.hyperplane.plot(label="contour2D.png", mode="plain", zlim=200)
    cm.analyze.hyperplane.analyzeStuff()
    cm.analyze.hyperplane.getMinima()
    cm.analyze.hyperplane.writeMinima()
    os.chdir("..")
    shutil.rmtree("hyperplane")
    os.chdir(trueRoot)

def test_uspexutils():
    trueRoot = os.getcwd()
    os.chdir("cases/uspex")
    cm.analyze.uspexutils.splitter(showDensity=True)
    arrays = cm.analyze.uspexutils.readInput()
    cm.analyze.uspexutils.plot(arrays, yLims=[(-0.05, 0.1), None, (-0.03, 0.1)], xLims=[None, None, (0.95, 1.025)])
    cm.analyze.uspexutils.plot(arrays, yLims=[(-0.05, 0.5), None, (-0.03, 2)])
    cm.analyze.uspexutils.plot(arrays)
    shutil.rmtree("poscar_split")
    for file in ["Density-Enthalpy_ID.pdf", "Enthalpy_Density.pdf", "Density_ID.pdf", "Enthalpy_ID.pdf"]: os.remove(file)
    os.chdir(trueRoot)
