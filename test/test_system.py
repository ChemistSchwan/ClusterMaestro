import ClusterMaestro as cm
from ase.io import read
from copy import deepcopy
import os

def test_analyse():
    molecule1 = cm.System.molecule.Molecule(read("cases/AdPh4.xyz", format="xyz"))
    assert len(molecule1) == 66, "Correct number of atoms."
    molecule1.find_core(25)
    molecule1.find_substituents()
    assert len(molecule1.core) == 22, "Correct number of atoms in core."
    assert len(molecule1.substituents) == 4, "Correct number of substituents."
    assert len(molecule1.substituents[0]) == 11, "Correct number of atoms in substituents."


    molecule2 = cm.System.molecule.Molecule(read("cases/SnSNp4.xyz", format="xyz"))
    assert len(molecule2) == 78, "Correct number of atoms."
    molecule2.analyze(core="SnS", substituent="Np")
    #molecule2.find_core(14)
    #molecule2.find_substituents()
    assert len(molecule2.core) == 10, "Correct number of atoms in core."
    assert len(molecule2.substituents) == 4, "Correct number of substituents."
    assert len(molecule2.substituents[0]) == 17, "Correct number of atoms in substituents."

    dimer1 = cm.System.Oligomer(read("cases/SiSPh4_dimer.xyz", format="xyz"))
    assert len(dimer1) == 108, "Correct number of atoms."
    dimer1.analyze(core="SiS", substituent="Ph")
    assert len(dimer1.molecules) == 2, "Correct number of molecules."
    assert len(dimer1.molecules[0].core) == 10 and len(dimer1.molecules[1].core) == 10, "Correct number of atoms in core."
    assert len(dimer1.molecules[0].substituents) == 4 and len(dimer1.molecules[0].substituents), "Correct number of substituents."
    assert len(dimer1.molecules[0].substituents[0]) == 11 and len(dimer1.molecules[0].substituents[1]), "Correct number of atoms in substituents."

def test_core_sub():
    mol = cm.System.molecule.Molecule(read("cases/AdPh4.xyz", format="xyz"))
    mol.find_core(25)
    mol.core.volume()
    mol.core.getCoreBondlength()

def test_structure():
    mol = cm.System.molecule.Molecule(read("cases/AdPh4.xyz", format="xyz"))
    molDel = deepcopy(mol)
    molDel.removeAtom(0)
    mol.find_core(25)
    mol.find_substituents()
    mol.substituents[0].get_connection(mol.core)
    mol.connectivityMatrix()
    mol.torsion(mol.substituents[0], 90)
    #mol.plot()
    mol.appendToXYZ("cases/dummy.xyz")
    mol.aromatic_carbons()
    mol.saturateCH()

def test_molecule():
    mol = cm.System.molecule.Molecule(read("cases/AdPh4.xyz", format="xyz"))
    mol.find_core(25)
    mol.find_substituents()
    mol.aromatic_carbons()
    mol.replace_core("Si")

def test_dimer():
    dimer1 = cm.System.Oligomer(read("cases/SiSPh4_dimer.xyz", format="xyz"))
    dimer1.analyze(core="SiS", substituent="Ph")
    dimer1.rm_outer_ph()
    dimer1.get_core_distance()
    dimer1.set_core_distance(9)
    dimer1.mean_pairwise_angle()
    dimer1.bonding_substituents()
    dimer1.get_core_tags()
    dimer1.matchPattern(core="SiS", substituent="Ph")

def test_solid():
    solid = cm.System.solid.Solid(read("cases/solid.vasp", format="vasp"))
    solid.analyze()
    solid.analyze(atomsPerMolecule=38)
    solid.expandSystem()
    solid.seperateMolecules(190)
    solid.completeMolecules_OLD(38)
    solid.expandMolecules()

def test_traj():
    traj = cm.System.trajectory.Trajectory(dataFile="cases/traj.xyz")
    traj.getBox("cases/traj.dat")
    traj.timeStep = 1
    print(traj)
    traj.write("cases/traj_write.xyz")
#    traj.cleanTrajectory(core="Ad", substituent="Ph")
    for i, frame in enumerate(traj.frames):
        frame.matchPattern(core="Ad", substituent="Ph")
        traj.frames[i] = frame.analyzeParticle()
    traj.shiftMols(66)
    traj.shiftMolsFAST(66)
    traj.getRMSD()
    traj.getRMSDfast()
    traj.plotRMSD()
    traj.getDensities()
    os.remove("cases/traj_write.xyz")
    os.remove("RMSD.pdf")



