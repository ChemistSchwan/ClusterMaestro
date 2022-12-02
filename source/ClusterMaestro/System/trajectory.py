"""The Trajectory-class can handle the analysis of an MD simulation with Dimer or Solid Objects.
"""

from copy import deepcopy
import numpy as np
import os, shutil, ase
from ase.io import read, write
import matplotlib.pyplot as plt
import ClusterMaestro as cm


class Trajectory(list):
    """
    Class for handling of trajectories returned from e.g. simulations with Lammps.
    Contains essentially a list of the different timesteps, where each listentry
    contains a Dimer- or Solid-object.
    """

    def __init__(
        self, dataFile=None, frames=[], box=None, timeStep=1, temperature=None
    ):
        if box == None:
            self.box = []
        else:
            self.box = box
        self.frames = frames
        if dataFile != None:
            self.dataFile = dataFile
            self._read(dataFile)
        self.timeStep = timeStep
        self.temperature = np.full([len(self.frames)], np.nan)
        try:
            self._getTemp()
        except:
            print("Cound not fill temperature-values.")

    def __repr__(self):
        return_string = f"CMTrajectory: Frames: {len(self.frames)}"
        return_string += f", Atoms per Frame: {len(self.frames[0])}"
        return_string += f", Timestep per Frame: {self.timeStep} ps"
        return return_string

    def _read(self, fileName):
        files = []
        firstLine = ""
        xyzString = ""
        with open(fileName, "r") as source:
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
            with open(".tmp_traj.xyz", "w") as target:
                target.write(string)
            if len(self.box) != 0:
                self.frames.append(
                    cm.System.solid.Solid(read(".tmp_traj.xyz", format="xyz"))
                )
                self.frames[-1].cell = self.box
            else:
                self.frames.append(
                    cm.System.Oligomer(read(".tmp_traj.xyz", format="xyz"))
                )

    def _getTemp(self):
        """
        Writes the temperature from the log.lammps-file.

        .. Warning::
           Only reliable, if the dump-command and the thermo command have the same value
           for output. If this is not the case, this function will interpolate the
           missing values evenly. The best would be to give dump and thermo the same value.
        """
        if os.path.exists("log.lammps"):
            with open("log.lammps", "r") as source:
                reader = False
                tempList = []
                for line in source:
                    if "dump" in line:
                        dumpStep = int(line.split()[4])
                    if "thermo" in line:
                        thermoStep = int(line.split()[1])
                    if "Loop" in line:
                        reader = False
                    if reader and len(line.split()) == 6:
                        tempList.append(float(line.split()[1]))
                    if "Step Temp E_pair E_mol TotEng Press" in line:
                        reader = True
            for i in range(0, len(self.frames), int(round(thermoStep / dumpStep, 0))):
                self.temperature[i] = tempList[int(i / round(thermoStep / dumpStep, 0))]
        nans, newTempList = cm.util.linalg_lib.nan_helper(self.temperature)
        self.temperature[nans] = np.interp(
            newTempList(nans), newTempList(~nans), self.temperature[~nans]
        )

    def write(self, fileName):
        """
        Writes the trajectory in a file in xyz-format.

        Parameters
        ----------

        fileName: str
            Name of resulting file.
        """
        if os.path.exists(fileName):
            os.remove(fileName)
        for i, system in enumerate(self.frames):
            with open(fileName, "a") as target:
                target.write("%d\n" % len(system))
                target.write(
                    "frame: %d of %d, timestep: %s ps, box: %.3f %.3f %.3f\n"
                    % (
                        i,
                        len(self.frames),
                        self.timeStep,
                        system.cell[0, 0],
                        system.cell[1, 1],
                        system.cell[2, 2],
                    )
                )
                for j in system:
                    target.write(
                        "%2s%17.8f%17.8f%17.8f%9d\n"
                        % (j.symbol, j.position[0], j.position[1], j.position[2], j.tag)
                    )

    def shiftMols(self, atomsPerMol):
        """
        Cleans the remaining molecule-jumps from the *cleanTrajectory* Function.

        Parameters
        ----------

        atomsPerMol: int
            How many atoms do we have per dedecated molecule in the system?
        """
        for i in range(len(self.frames)):
            if i != 0:
                if cm.util.rmsd.simpleRMSD(self.frames[i], self.frames[i - 1]) > 2:
                    print("\n")
                    print("Shift of a molecule from timestep %d to %d" % (i - 1, i))
                    print("Running analyze-function...")
                    sys1 = deepcopy(self.frames[i])
                    sys2 = deepcopy(self.frames[i - 1])
                    sys1.analyze(atomsPerMol)
                    sys2.analyze(atomsPerMol)
                    for j in range(len(sys1.molecules)):
                        if (
                            cm.util.rmsdFAST.simpleRMSD(
                                sys1.molecules[j], sys2.molecules[j]
                            )
                            > 1
                        ):
                            center_1 = sys1.molecules[j].center_of_mass()
                            center_2 = sys2.molecules[j].center_of_mass()
                            centerOfMasses = center_1 - center_2
                            for k, val in enumerate(centerOfMasses):
                                if abs(val) > 5:
                                    if val > 0:
                                        for l in range(i, len(self.frames)):
                                            self.frames[l].positions[
                                                j * atomsPerMol : (j + 1) * atomsPerMol,
                                                k,
                                            ] -= self.box[k]
                                    if val < 0:
                                        for l in range(i, len(self.frames)):
                                            self.frames[l].positions[
                                                j * atomsPerMol : (j + 1) * atomsPerMol,
                                                k,
                                            ] += self.box[k]
                                    print("Shift of molecule %d was corrected." % j)

    def shiftMolsFAST(self, atomsPerMol):
        """
        Cleans the remaining molecule-jumps from the *cleanTrajectory* Function.

        Compiled code using Cython.

        Parameters
        ----------

        atomsPerMol: int
            How many atoms do we have per dedecated molecule in the system?
        """
        self = cm.lib.solidFuncs.shiftMols(self, atomsPerMol)

    def getRMSD(self):
        """
        Returns an 1D array containing the RMSD values for each timestep. It is highly
        recommended to use this function on the trajectory file, which was created by
        the *cleanTrajectory*-function.

        Parameters
        ----------

        trajFile: str
            Trajectoryfile, which will be analysed.
        """
        rmsdList = []
        for i, frame in enumerate(self.frames):
            if i == 0:
                firstSystem = deepcopy(frame)
            rmsdList.append(cm.util.rmsd.simpleRMSD(firstSystem, frame))
        return np.array(rmsdList)

    def getRMSDfast(self):
        """
        Returns an 1D array containing the RMSD values for each timestep. It is highly
        recommended to use this function on the trajectory file, which was created by
        the *cleanTrajectory*-function.

        .. note::
           Same code as *getRMSD*, but compiled using Cython, resulting in about 30\%
           better performance.

        Parameters
        ----------

        trajFile: str
            Trajectoryfile, which will be analysed.
        """
        rmsdList = []
        for i, frame in enumerate(self.frames):
            if i == 0:
                firstSystem = deepcopy(frame)
            rmsdList.append(cm.util.rmsdFAST.simpleRMSD(firstSystem, frame))
        return np.array(rmsdList)

    def plotRMSD(self, filename="RMSD.pdf", rmsdList=None, show=False):
        """
        Simple plot-function for lists returned from the getRMSD-function.

        Parameters
        ----------

        filename: str
            Name of the resulting imagefile.
        rmsdList: array
            List with the rmsd-values. Will be generated automatically, if not given.
        show: bool
            Whether or not the image will just be shown, or saved.
        """
        import matplotlib.pyplot as plt

        if self.timeStep == 0:
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
            print("++                     WARNING                     ++")
            print("++            timeStep is still set to 0           ++")
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
        plt.clf()
        if rmsdList == None:
            rmsdList = self.getRMSD()
        xarray = np.arange(0, len(rmsdList), 1)
        xarray = xarray * self.timeStep

        plt.plot(xarray, rmsdList)
        plt.ylabel("RMSD / $\AA$")
        plt.xlabel("time in ps")
        if show:
            plt.show()
        else:
            plt.savefig(filename)
        plt.clf()

    def plotTEMP(self, filename="TEMP.pdf", tempList=None, show=False):
        """
        Simple plot-function for lists returned from the getRMSD-function.

        Parameters
        ----------

        filename: str
            Name of the resulting imagefile.
        tempList: array
            List with the temperature-values. Will be generated automatically, if not given.
        show: bool
            Whether or not the image will just be shown, or saved.
        """
        if self.timeStep == 0:
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
            print("++                     WARNING                     ++")
            print("++            timeStep is still set to 0           ++")
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
        import matplotlib.pyplot as plt

        plt.clf()
        tempList = deepcopy(self.temperature)
        nans, newTempList = cm.util.linalg_lib.nan_helper(tempList)
        tempList[nans] = np.interp(
            newTempList(nans), newTempList(~nans), tempList[~nans]
        )
        xarray = np.arange(0, len(tempList), 1)
        xarray = xarray * (self.timeStep)
        plt.plot(xarray, tempList)
        plt.ylabel("Temperature / K")
        plt.xlabel("time in ps")
        if show:
            plt.show()
        else:
            plt.savefig(filename)
        plt.clf()

    def plotRT(
        self, filename="RMSD_TEMP.pdf", tempList=None, rmsdList=None, show=False
    ):
        """
        Simple plot-function for lists returned from the getRMSD-function.

        Parameters
        ----------

        filename: str
            Name of the resulting imagefile.
        tempList: array
            List with the temperature-values. Will be generated automatically, if not given.
        rmsdList: array
            List with the rmsd-values. Will be generated automatically, if not given.
        show: bool
            Whether or not the image will just be shown, or saved.
        """
        if self.timeStep == 0:
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
            print("++                     WARNING                     ++")
            print("++            timeStep is still set to 0           ++")
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
        if rmsdList == None:
            rmsdList = self.getRMSD()
        import matplotlib.pyplot as plt

        plt.clf()
        xarray = np.arange(0, len(rmsdList), 1)
        xarray = xarray * self.timeStep

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        color1 = "#005b9a"
        color2 = "#c33232"

        ax1.plot(xarray, rmsdList, color=color1)
        ax2.plot(xarray, self.temperature, color=color2)
        ax1.set_ylabel("RMSD / $\AA$", color=color1)
        ax2.set_ylabel("Temperature / K", color=color2)
        ax1.set_xlabel("time in ps")
        ax1.set_xlim(xmin=0)
        ax2.set_xlim(xmin=0)
        ax1.set_ylim(ymin=0)
        ax2.set_ylim(ymin=0)
        ax1.tick_params(axis="y", labelcolor=color1)
        ax2.tick_params(axis="y", labelcolor=color2)
        ax2.set_yticks(
            np.linspace(
                ax2.get_yticks()[0], ax2.get_yticks()[-1], len(ax1.get_yticks())
            )
        )
        ax1.set_yticks(
            np.linspace(
                ax1.get_yticks()[0], ax1.get_yticks()[-1], len(ax2.get_yticks())
            )
        )
        ax1.grid(axis="y", zorder=10)
        plt.subplots_adjust(left=0.1, right=0.85, top=0.98)
        if show:
            plt.show()
        else:
            plt.savefig(filename)

    def getBox(self, dataFile):
        """
        Extracts the box from a lammps-input file.

        Parameters
        ----------

        dataFile: str
            Lammps-inputfile (endung: '.dat')
        """
        box = []
        with open(dataFile, "r") as source:
            for line in source:
                if "lo" in line:
                    splitLine = line.split()
                    box.append(float(splitLine[1]))
        self.box = box

    def getDensities(self):
        """
        Returns a list (and also prints) all densities at each individual frame.
        """
        denList = []
        if type(self.frames[0]) == cm.System.solid.Solid:
            for i, frame in enumerate(self.frames):
                dens = frame.calcDensity()
                denList.append(dens)
                print("Frame: %d, density: %f g/cm^3" % (i, dens))
        elif type(self.frames[0]) == cm.System.Oligomer:
            for i, frame in enumerate(self.frames):
                try:
                    dens = frame.density
                    denList.append(dens)
                except:
                    system = frame.analyzeParticle()
                    dens = system.density
                    denList.append(dens)
                print("Frame: %d, density: %f g/cm^3" % (i, dens))

        return np.array(denList)

    def plotDensities(self):
        """
        Creates a PDF called "LMP_densities.pdf" with all the densities according to the
        timesteps.
        """
        plt.clf()
        denList = self.getDensities()
        plt.plot(
            cm.util.general.frange(0, (len(denList) * self.timeStep), self.timeStep),
            denList,
        )
        plt.xlabel("time / ps")
        plt.ylabel("Density / g/cm$^3$")
        plt.savefig("LMP_densities.pdf")

    def makeDense(self, atomsPerMolecule=None, core=None, substituent=None):
        """
        Lammps often returns the structure in a very un-dense manner. In order to extract
        e.g. Dimers from this structure, the molecules need to be shifted into the box,
        which is done with this function.

        Parameters
        ----------

        atomsPerMolecule: int
            see Solid.analyze
        core: Core
            see Solid.analyze
        substituent: Substituent
            see Solid.analyze
        """

        def pointInBox(point, box):
            returnVal = []
            for i, val in enumerate(point):
                if val < 0 or val > box[i, i]:
                    returnVal.append(i)
            if len(returnVal) == 0:
                return [True, None]
            else:
                return [False, returnVal]

        def shiftMolVal(point, box, vals):
            newPoint = np.zeros((3))
            for i in vals:
                if point[i] < box[i, i]:
                    newPoint[i] = +box[i, i]
                elif point[i] > box[i, i]:
                    newPoint[i] = -box[i, i]
            return newPoint

        def shiftMolInBox(molecule, shiftVal):
            for i in molecule:
                i.position += shiftVal
            return molecule

        for i, frame in enumerate(self.frames):
            print("Frame number {} is being analyzed.".format(i))
            frame.analyze(
                atomsPerMolecule=atomsPerMolecule, core=core, substituent=substituent
            )
            for j, molecule in enumerate(frame.molecules):
                inBox = pointInBox(molecule.center_of_mass(), frame.cell)
                if not inBox[0]:
                    shiftVal = shiftMolVal(
                        molecule.center_of_mass(), frame.cell, inBox[1]
                    )
                    self.frames[i].molecules[j] = shiftMolInBox(molecule, shiftVal)
            frame._renew_positions_mols()

    def makeDenseParallel(
        self, atomsPerMolecule=None, core=None, substituent=None, threads=None
    ):
        """
        Lammps often returns the structure in a very un-dense manner. In order to extract e.g.
        Dimers from this structure, the molecules need to be shifted into the box, which is
        done with this function.

        Parameters
        ----------

        atomsPerMolecule: int
            see Solid.analyze
        core: Core
            see Solid.analyze
        substituent: Substituent
            see Solid.analyze
        threads: int or None
            How many threads should be used for parallelization? Default: output of
            multiprocessing.cpu_count() / 2
        """
        from joblib import Parallel, delayed
        from multiprocessing import cpu_count

        def pointInBox(point, box):
            returnVal = []
            for i, val in enumerate(point):
                if val < 0 or val > box[i, i]:
                    returnVal.append(i)
            if len(returnVal) == 0:
                return [True, None]
            else:
                return [False, returnVal]

        def shiftMolVal(point, box, vals):
            newPoint = np.zeros((3))
            for i in vals:
                if point[i] < box[i, i]:
                    newPoint[i] = +box[i, i]
                elif point[i] > box[i, i]:
                    newPoint[i] = -box[i, i]
            return newPoint

        def shiftMolInBox(molecule, shiftVal):
            for i in molecule:
                i.position += shiftVal
            return molecule

        def getMolsright(frame):
            frame.analyze(
                atomsPerMolecule=atomsPerMolecule, core=core, substituent=substituent
            )
            for j, molecule in enumerate(frame.molecules):
                inBox = pointInBox(molecule.center_of_mass(), frame.cell)
                if not inBox[0]:
                    shiftVal = shiftMolVal(
                        molecule.center_of_mass(), frame.cell, inBox[1]
                    )
                    frame.molecules[j] = shiftMolInBox(molecule, shiftVal)
            frame._renew_positions_mols()
            return frame

        if threads != None:
            cores = threads
        else:
            cores = cpu_count() / 2
        allFrames = Parallel(n_jobs=int(cores))(
            delayed(getMolsright)(frame) for frame in self.frames
        )
        for i in range(len(self.frames)):
            self.frames[i] = allFrames[i]

    def extractDimers(
        self,
        core="Ad",
        substituent="Ph",
        optimize=False,
        acc="normal",
        threads=1,
        GFN=2,
        maxCoreDist=12,
        frameNumbers="all",
        ignoreOuter=False,
    ):
        """
        Uses the *ClusterMaestro.System.solid.Solid.extractDimers* function to
        extract dimers from EACH frame of the system. Those will be renamed accordingly.

        Parameters
        ----------

        core: str
            See :ref:`Structure parts`
        substituent: str
            See :ref:`Structure parts`
        optimize: bool
            Whether or not the extracted dimer structure will be optimized
            (using xTB and their default parameters).
        acc: str
            See *ClusterMaestro.util.optimize*.
        threads: int
            See *ClusterMaestro.util.optimize*.
        GFN: int, "ff" or list
            See *ClusterMaestro.util.optimize*. If *GFN* is a List, a sequence will be made,
            where the *GFN* parameters are executed one by one.
        maxCoreDist: float
            Maximum accepted core-core distance, which is still acceptible as a dimer.
        frameNumbers: list or "all"
            Which frames should be analysed? Default value is "all", which analyses all
            frames. If this function is used, it is NOT mandatory, that the dimers are
            in the right order!
        ignoreOuter: bool
            Whether or not outer molecules are beeing ignored. This can only be used for
            *Dimer*-type frames, as *Solid*-type frames have periodic boundary conditions
            and therefore no outer molecules.
        """
        if not os.path.exists("allDimers"):
            os.mkdir("allDimers")
            dimerCounter = 1
        else:
            fileListAD = os.listdir("allDimers/")
            if len(fileListAD) == 0:
                shutil.rmtree("allDimers")
                os.mkdir("allDimers")
                dimerCounter = 1
            else:
                currentNumber = int(fileListAD[-1][:-4])
                dimerCounter = currentNumber
        if frameNumbers == "all":
            for frameNumber, i in enumerate(self.frames):
                for j in i.molecules:
                    j.matchPattern(core=core, substituent=substituent)
                if ignoreOuter:
                    dimerSystem = i.analyzeParticle()
                    dimerSystem.extractDimers(
                        optimize=optimize,
                        acc=acc,
                        threads=threads,
                        GFN=GFN,
                        maxCoreDist=maxCoreDist,
                    )
                else:
                    i.extractDimers(
                        optimize=optimize,
                        acc=acc,
                        threads=threads,
                        GFN=GFN,
                        maxCoreDist=maxCoreDist,
                    )
                fileList = os.listdir("dimers")
                fileList.sort()
                for fileName in fileList:
                    shutil.copyfile(
                        "dimers/%s" % fileName, "allDimers/%.5d.xyz" % (dimerCounter)
                    )
                    dimerCounter += 1
                shutil.rmtree("dimers")
        else:
            for frameNumber, i in enumerate(self.frames):
                if frameNumber not in frameNumbers:
                    continue
                for j in i.molecules:
                    j.matchPattern(core=core, substituent=substituent)
                if ignoreOuter:
                    dimerSystem = i.analyzeParticle()
                    dimerSystem.extractDimers(
                        optimize=optimize,
                        acc=acc,
                        threads=threads,
                        GFN=GFN,
                        maxCoreDist=maxCoreDist,
                    )
                else:
                    i.extractDimers(
                        optimize=optimize,
                        acc=acc,
                        threads=threads,
                        GFN=GFN,
                        maxCoreDist=maxCoreDist,
                    )
                fileList = os.listdir("dimers")
                fileList.sort()
                print(fileList)
                for fileName in fileList:
                    shutil.copyfile(
                        "dimers/%s" % fileName, "allDimers/%.5d.xyz" % (dimerCounter)
                    )
                    dimerCounter += 1
                shutil.rmtree("dimers")

    def getNextNeighbors(
        self, rMax=12.0, frames="all", plot=False, plotName="NN_Hist.pdf"
    ):
        """
        Analyses a trajectory in regards to the next neighbors of each molecule. The
        trajectory must be analysed.

        Parameters
        ----------

        rMax: float
            Every molecule in a radius of r <= rMax is considered a neighbor. Please choose
            this function wisely by using the *pairDistFuncMolecules* function.
        frames: "all" or range
            Which frames are beeing taken into account? For longer trajectories, it might be
            wise to ignore the first couple frames.
        plot: bool
            If True, a histogram plot will be created.
        plotName: str
            Name of resulting histogram plot.
        """
        resultNNarray = []  # np.zeros(len(self.frames[0].molecules))
        counter = 0
        if frames == "all":
            for i in range(len(self.frames)):
                for k in self.frames[i].getNextNeighbors(rMax=rMax):
                    resultNNarray.append(k)
                counter += 1
                print(counter)
        else:
            frameRange = frames
            for i in frameRange:
                for k in self.frames[i].getNextNeighbors(rMax=rMax):
                    resultNNarray.append(k)
                counter += 1
                print(counter)
        resultNNarray = np.array(resultNNarray)
        if plot:
            n_bins = 15
            fig, ax = plt.subplots()
            ax.hist(resultNNarray, bins=n_bins, align="left", rwidth=1, range=(0, 15))
            plt.xlim(0, 15)
            plt.xticks(range(0, 21, 2), range(0, 21, 2))
            plt.xlabel("Number of next neighbors")
            plt.ylabel("Appearances of this number of NN")
            plt.savefig(plotName)
        return resultNNarray

    def pairDistFuncMolecules(
        self,
        frames=[],
        rMax=10,
        rMin=0,
        dr=0.1,
        plotName="pairDistFuncMolecule_traj.pdf",
    ):
        """
        This function plots the pair distribution function of the given System with
        respect to whole molecules and NOT atoms. It is based on
        *https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation*.

        Returns also the return arrays from the *pairCorrelationFunction_3D*.

        The solid needs to be analysed already!

        Parameters
        ----------

        frames: list
            Which frames will be plotted.
        rMax: float
            Maximum radius to which the PDF will be created.
        plotName: str
            Name of the resulting plot.
        """
        import matplotlib.pyplot as plt
        from matplotlib import style

        style.use("seaborn")
        plt.figure(figsize=(4, 3))
        colorList = cm.util.constants.colorList_2
        counter = 0
        for frame in frames:
            positions = np.empty((len(self.frames[frame].molecules), 3))
            for i, molecule in enumerate(self.frames[frame].molecules):
                positions[i] = molecule.center_of_mass()
            dummy = ase.Atoms()
            for i in positions:
                dummy += ase.Atom(symbol="Xe", position=i)
            newSolid = cm.System.solid.Solid(dummy)
            newSolid.cell = self.frames[frame].cell

            while 4 * rMax > newSolid.cell[0, 0]:
                print("exapnding cell for pairDistFuncMolecules ...")
                newSolid = newSolid.expandSystem()

            g_r, r, reference_indices = cm.util.linalg_lib.pairCorrelationFunction_3D(
                newSolid.positions[:, 0],
                newSolid.positions[:, 1],
                newSolid.positions[:, 2],
                newSolid.cell[0, 0],
                rMax,
                dr,
            )

            # Plot defining region:
            plt.plot(r, g_r, label=frame, color=colorList[counter])
            counter += 1
        plt.xlabel("Radius / ${\AA}$")
        plt.ylabel("G(r)")
        plt.xlim((rMin, rMax))
        plt.ylim((0, 1.05 * g_r.max()))
        plt.tight_layout()
        plt.legend()
        plt.savefig(plotName)

        return g_r, r, reference_indices

    def viewNeighborPositions(
        self,
        reference,
        core="Ad",
        substituent="Ph",
        maxCoreDist=10,
        write=True,
        frameNumbers="all",
    ):
        """
        Uses *rotateDimerOnRef* to fit all two-molecule systems of the *Dimer* and created
        dummy-helium molecules at the positions of their center of mass. If the folder
        *allDimers* with extracted dimer structures already
        exists, this function will use the structures from that folder.

        Parameters
        ----------

        reference: Dimer
            Reference-molecule, which will be used for the rotation.
        maxCoreDist: float
            Maximum accepted core-core distance, which is still acceptible as a dimer.
        write: bool
            Whether or not the files "centers.xyz" and "referenceMol.xyz" are written.
        frameNumbers: list or "all"
            Which frames should be analysed? Default value is "all", which analyses all
            frames. If this function is used, it is NOT mandatory, that the dimers are
            in the right order!
        """
        if not os.path.exists("allDimers"):
            self.extractDimers(
                threads=4, maxCoreDist=maxCoreDist, frameNumbers=frameNumbers
            )
        centers = ase.Atoms()
        counter = 0
        for dimerNum in os.listdir("allDimers"):
            print(counter)
            dimer = cm.System.Dimer(read("allDimers/{}".format(dimerNum)))
            dimer.matchPattern(core=core, substituent=substituent)
            rotDimer = dimer.rotateDimerOnRef(reference)
            center = ase.Atom(
                symbol="H", position=rotDimer.molecules[1].center_of_mass()
            )
            centers += center
            counter += 1
        referenceMol = reference.molecules[0]

        return centers, referenceMol

    def getDistributionPlot(
        self, reference, maxCoreDist=10, core="Ad", substituent="Ph", frameNumbers="all"
    ):
        """
        Creates a 2D plot for the distribution of other molecules around a
        reference molecule. If the function *viewNeighborPositions* was executed before and
        the return objects were saved as "centers.xyz" and "referenceMol.xyz", this function
        will not perform this again, but read the two files.

        Parameters
        ----------

        maxCoreDist: float
            Maximum accepted core-core distance, which is still acceptible as a dimer.
        reference: Dimer
            Reference-molecule, which will be used for the rotation.
        frameNumbers: list or "all"
            Which frames should be analysed? Default value is "all", which analyses all
            frames. If this function is used, it is NOT mandatory, that the dimers are
            in the right order!
        """
        import matplotlib.pyplot as plt
        from matplotlib import colors

        if os.path.exists("centers.xyz") and os.path.exists("referenceMol.xyz"):
            surroundings = read("centers.xyz")
            ref = cm.System.Molecule(read("referenceMol.xyz"))
        else:
            surroundings, ref = self.viewNeighborPositions(
                reference,
                core=core,
                substituent=substituent,
                maxCoreDist=maxCoreDist,
                frameNumbers=frameNumbers,
            )
            surroundings.write("centers.xyz")
            ref.write("referenceMol.xyz")
        vals = [[], []]
        distances = []
        for i in surroundings:
            distances.append(np.linalg.norm(i.position))
        distances = np.array(distances)
        for i in surroundings:
            test = cm.util.linalg_lib.Point(i.position[0], i.position[1], i.position[2])

            sph = test.toSpherical()

            sphCoords = sph.degrees()
            vals[0].append(sphCoords[1])
            vals[1].append(sphCoords[2])
        c = "C0"
        for i, val in enumerate(distances):
            if val > maxCoreDist:
                vals[0][i] = np.nan
                vals[1][i] = np.nan
        fig, axs = plt.subplots(figsize=(4, 3), tight_layout=True)
        plt.scatter(
            vals[0], vals[1], s=2, c=distances, cmap=plt.get_cmap("cividis"), alpha=0.7
        )
        plt.colorbar(label="core-core distance / $\AA$")
        plt.xlabel(r"polar angle $\theta$ / °")
        plt.ylabel(r"azimut angle $\varphi$ / °")
        plt.xlim(0, 180)
        plt.ylim(0, 360)
        plt.savefig("angleDistrib.pdf")
        plt.savefig("angleDistrib.png", dpi=500)

    def getBondLengthDist(self, part="core", atoms=["C", "C"], frames="all"):
        """
        Returns a distribution of the bondlengths within the core-structure of the system.
        The system needs to be analyzed.

        Parameters
        ----------

        part: str
            Currently: "core" or "substituent". This will only take one of the
            structure parts into account.
        atoms: list
            Only bonds between these two atoms types are considered.
        frames: list or "all"
            Which frames should be considered?
        """
        if frames == "all":
            frames = range(len(self.frames))
        totalList = []
        for i in frames:
            newList = self.frames[i].getBondLengthDist(part=part, atoms=atoms)
            totalList = np.concatenate((totalList, newList))
        totalList = np.array(totalList)
        return totalList

    def getBondAngles(self, part="core", atoms=["C", "C", "C"], frames="all"):
        """
        Returns a distribution of the bondlengths within the core-structure of the system.
        The system needs to be completely analyzed.

        Parameters
        ----------

        part: str
            Currently: "core" or "substituent". This will only take one of the structure
            parts into account.
        atoms: list
            Only bonds between these two atoms types are considered.
        frames: list or "all"
            Which frames should be considered?
        """
        if frames == "all":
            frames = range(len(self.frames))
        totalList = []
        for i in frames:
            newList = self.frames[i].getBondAngles(part=part, atoms=atoms)
            totalList = np.concatenate((totalList, newList))
        totalList = np.array(totalList)
        return totalList

    def getBondLengthCoreSub(self, frames):
        """
        Returns a distribution of the bondlengths of the core-substituent
        bond of the system. The system needs to be completely analyzed.

        Parameters
        ----------

        frames: list or "all"
            Which frames should be considered?
        """
        if frames == "all":
            frames = range(len(self.frames))
        totalList = []
        for i in frames:
            newList = self.frames[i].getBondLengthCoreSub()
            totalList = np.concatenate((totalList, newList))
        totalList = np.array(totalList)
        return totalList
