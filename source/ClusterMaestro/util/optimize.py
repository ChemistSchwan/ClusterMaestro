import os, re
import shutil
import ClusterMaestro as cm
from ase.io import read

"""On-the-fly optimization of chemical systems.

    Further analysis has to be performed seperately after the optimization.
"""


def opt_xtb(
    self,
    acc="normal",
    threads=1,
    GFN=2,
    output=False,
    showString=False,
    folder=None,
    inputFile=None,
):
    """Optimization using the xTB-GFN programm in version by Grimme et al.

    prerequisite: working xtb-version in bin-directory. (Download: https://github.com/grimme-lab/xtb/releases)
    Also view the xtb-Documentation for more precize information. (https://xtb-docs.readthedocs.io/en/latest/contents.html).

    Please cite the original xtb-publications when using. These can vary depending on the xTB-version you are using.

    Parameters
    ----------

    acc: str
        Level of accuracy. Default: normal
        Options: crude, sloppy, loose, lax, normal, tight, vtight, extreme
    threads: int
        Number of threads the optimization is running on. Default: 1
    GFN: inti or "ff"
        GFN-algorithm. Default: 2
    output: Bool
        Whether or not the file "xtbopt.out" is beeing written and stored. Defalt: False
    showString: Bool
        Whether or not the string for the xtbinput is beeing printed.
    folder: str or None
        If not None: keep folder, where the optimization is performed. str is the name of that folder
    inputFile: str or None
        If not None: File with detailed input (see xTB Documentation)
    """

    trueRoot = os.getcwd()
    if folder == None:
        folder = ".opt_tmp"
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)
    os.chdir(folder)
    if inputFile == None:
        if GFN == "ff":
            xtbString = "xtb --gfnff __file__ --opt __acc__ -P __threads__ __output__"
        else:
            xtbString = (
                "xtb __file__ --opt __acc__ -P __threads__ --gfn __GFN__ __output__"
            )
            xtbString = re.sub("__GFN__", str(GFN), xtbString)
    else:
        shutil.copyfile("../{}".format(inputFile), inputFile)
        if GFN == "ff":
            xtbString = "xtb --gfnff __file__ --opt __acc__ -P __threads__ __output__ --input __inputFile__"
        else:
            xtbString = "xtb __file__ --opt __acc__ -P __threads__ --gfn __GFN__ __output__ --input __inputFile__"
            xtbString = re.sub("__GFN__", str(GFN), xtbString)
            xtbString = re.sub("__inputFile__", str(inputFile), xtbString)

    xtbString = re.sub("__file__", "start.xyz", xtbString)
    xtbString = re.sub("__acc__", acc, xtbString)
    xtbString = re.sub("__threads__", str(threads), xtbString)
    xtbString = re.sub("__output__", " > xtbopt.out", xtbString)
    if showString:
        print(xtbString)
    self.write("start.xyz", format="xyz")
    os.system(xtbString)
    if output:
        shutil.copyfile("xtbopt.out", "../xtbopt.out")
    molecule = type(self)(read("xtbopt.xyz"))
    with open("xtbopt.out", "r") as outFile:
        for line in outFile:
            if "TOTAL ENERGY" in line:
                splitLine = line.split()
                for word in splitLine:
                    if cm.util.general.isFloat(word):
                        molecule.energy = float(word)
                        break
    os.chdir(trueRoot)
    if folder == ".opt_tmp":
        shutil.rmtree(".opt_tmp")

    return molecule


def xtb_constrain(filename, atomtags, forceConstant=1.0):
    """Create a constrain-file for an xtb-optimization.

    Parameters
    ----------

    filename: str
        Name of the constrain-file
    atomtags: List
        List of the atoms, whose positions should be constrained.
    forceConstant: float
        Force constant of the constrains.
    """
    constString = (
        "$constrain\n force constant=%.3f\n atoms: __atoms__\n$end" % forceConstant
    )
    atomString = ""
    for i in atomtags:
        atomString += "%d," % (i + 1)
    atomString = atomString[:-1]
    constString = re.sub("__atoms__", atomString, constString)
    with open(filename, "w") as target:
        target.write(constString)
