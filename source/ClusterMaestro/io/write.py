import ClusterMaestro.util.constants


def poscar(name, system, box=[]):
    """Function to write a VASP-style PASCAR-file from objects of the ClusterMaestro.System or its subclasses.

    Parameters
    ----------

    name: str
        Name of the resulting file
    system: Atoms
        Molecule, which should be written in the POSCAR-file
    box: 1x3 array
        Edge-length of the CUBIC box
    """
    atoms = []
    atomsNumbers = []
    positions = []
    for i in system:
        if i.symbol not in atoms:
            atoms.append(i.symbol)
            atomsNumbers.append(1)
            positions.append([i.position])
        else:
            for j, symbol in enumerate(atoms):
                if i.symbol == symbol:
                    atomsNumbers[j] += 1
                    positions[j].append(i.position)
    if box == []:
        box = system.cell
    with open(name, "w") as target:
        target.write("Created with ClusterMaestro\n   1.00000000000\n")
        if box.shape == (3,):
            target.write("%20.10f%15.10f%15.10f\n" % (box[0], 0.0, 0.0))
            target.write("%20.10f%15.10f%15.10f\n" % (0.0, box[0], 0.0))
            target.write("%20.10f%15.10f%15.10f\n" % (0.0, 0.0, box[0]))
        else:
            cell = system.cell
            target.write(
                "%20.10f%15.10f%15.10f\n" % (cell[0, 0], cell[0, 1], cell[0, 2])
            )
            target.write(
                "%20.10f%15.10f%15.10f\n" % (cell[1, 0], cell[1, 1], cell[1, 2])
            )
            target.write(
                "%20.10f%15.10f%15.10f\n" % (cell[2, 0], cell[2, 1], cell[2, 2])
            )
        for i in atoms:
            target.write("%5s" % i)
        target.write("\n")
        for i in atomsNumbers:
            target.write("%5d" % i)
        target.write("\nCartesian\n")
        for i in positions:
            for j in i:
                target.write("%20.12f%20.12f%20.12f\n" % (j[0], j[1], j[2]))


def tmPeriodic(name, system, box=[]):
    """Function to write a turbomole-style input-file from objects of the ClusterMaestro.System or its subclasses.
    Can also be used for xtb-calculations.

    Parameters
    ----------

    name: str
        Name of the resulting file
    system: Atoms
        Molecule, which should be written in the POSCAR-file
    box: 1x3 array
        Edge-length of the CUBIC box
    """
    with open(name, "w") as target:
        if len(box) == 3:
            target.write(
                "$periodic 3\n$cell angs\n   %8.3f%8.3f%8.3f  90  90  90\n$coord angs\n"
                % (box[0], box[1], box[2])
            )
        else:
            target.write(
                "$periodic 3\n$cell angs\n   %8.3f%8.3f%8.3f  90  90  90\n$coord angs\n"
                % (system.cell[0, 0], system.cell[1, 1], system.cell[2, 2])
            )
        for i in system:
            target.write(
                " %9.5f%9.5f%9.5f%3s\n"
                % (i.position[0], i.position[1], i.position[2], i.symbol)
            )
        target.write("$end")


def lammps(name, system, box, atom_style="charge", atomsPerMolecule=66):
    """Function to write a lammps data-file from objects of the ClusterMaestro.System or its subclasses.
    Can only be used with the 'atom_style   charge'.

    Parameters
    ----------

    name: str
        Name of the resulting file
    system: Atoms
        Molecule, which should be written in the POSCAR-file
    box: 1x3 array
        Edge-length of the CUBIC box
    atom_style: str
        Atom style in lammps-data file. Currently supported: charge and full.
    atomsPerMolecule: int
        See *ClusterMaestro.System.solid.Solid.analyze*. Only needed for atom_style *full* or *molecular*
    """

    def getAtomTypes(molecule):
        counter = 0
        typeList = []
        typeDict = {}
        for i in molecule:
            if i.symbol not in typeList:
                typeList.append(i.symbol)
                typeDict.update({i.symbol: counter + 1})
                counter += 1
        return [counter, typeList, typeDict]

    atomTypes, typeList, typeDict = getAtomTypes(system)
    if atom_style == "charge":
        with open(name, "w") as target:
            target.write("#Creaed with ClusterMaestro\n\n")
            target.write("%d atoms\n%d atom types\n\n" % (len(system), atomTypes))
            target.write(
                "%.1f %.1f xlo xhi\n%.1f %.1f ylo yhi\n%.1f %.1f zlo zhi\n\nMasses\n\n"
                % (0.0, box[0], 0.0, box[1], 0.0, box[2])
            )
            for i, symbol in enumerate(typeList):
                target.write(
                    "%d %f\n" % (i + 1, ClusterMaestro.util.constants.atomMasses[symbol])
                )
            target.write("\nAtoms\n\n")
            print(system[0].tag)
            for i, atom in enumerate(system):
                target.write(
                    "%4d%4d%5.1f%10.5f%10.5f%10.5f\n"
                    % (
                        atom.tag + 1,
                        typeDict[atom.symbol],
                        0.0,
                        atom.position[0],
                        atom.position[1],
                        atom.position[2],
                    )
                )
    elif atom_style == "full":
        with open(name, "w") as target:
            target.write("#Creaed with ClusterMaestro\n\n")
            target.write("%d atoms\n%d atom types\n\n" % (len(system), atomTypes))
            target.write(
                "%.1f %.1f xlo xhi\n%.1f %.1f ylo yhi\n%.1f %.1f zlo zhi\n\nMasses\n\n"
                % (0.0, box[0], 0.0, box[1], 0.0, box[2])
            )
            for i, symbol in enumerate(typeList):
                target.write(
                    "%d %f\n" % (i + 1, ClusterMaestro.util.constants.atomMasses[symbol])
                )
            target.write("\nAtoms\n\n")
            system.analyze(atomsPerMolecule=atomsPerMolecule)
            atomCounter = 0
            for i, molecule in enumerate(system.molecules):
                for j, atom in enumerate(molecule):
                    atomCounter += 1
                    target.write(
                        "%4d%4d%4d%5.1f%10.5f%10.5f%10.5f\n"
                        % (
                            atomCounter,
                            i + 1,
                            typeDict[atom.symbol],
                            atom.charge,
                            atom.position[0],
                            atom.position[1],
                            atom.position[2],
                        )
                    )


def xyzTrajectory(name, system):
    """
    Writes the pure xyz-coordinates into a simple xyz-trajectory file. No PBC are considered here!

    Parameters
    ----------

    name: str
        Name of the resulting file
    system: Trajectory
        Trajectory-Object, which will be written.
    """
    with open(name, "w") as target:
        for i, frame in enumerate(system.frames):
            target.write("%d\n" % len(frame))
            target.write("step: %d\n" % i)
            for j in frame:
                target.write(
                    "%3s%10.5f%10.5f%10.5f\n"
                    % (j.symbol, j.position[0], j.position[1], j.position[2])
                )


def lmpTrajectory(name, system):
    """
    Writes the Trajectory as .lmp-file in the format: *dump d1 all     custom 100000 traj.lmp element xu yu zu*

    Parameters
    ----------

    name: str
        Name of the resulting file.
    system: Trajectory
        Trajectory-Object, which will be written.
    """
    with open(name, "w") as target:
        for i, frame in enumerate(system.frames):
            target.write("ITEM: TIMESTEP\n%d\n" % (i * system.timeStep))
            target.write(
                "ITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n" % (len(frame))
            )
            target.write(
                "0  %.15f\n0  %.15f\n0  %.15f\n"
                % (frame.cell[0, 0], frame.cell[1, 1], frame.cell[2, 2])
            )
            target.write("ITEM: ATOMS element xu yu zu\n")
            for atom in frame:
                target.write(
                    "%2s%10.5f%10.5f%10.5f\n"
                    % (
                        atom.symbol,
                        atom.position[0],
                        atom.position[1],
                        atom.position[2],
                    )
                )
