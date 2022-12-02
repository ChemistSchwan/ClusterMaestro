import ClusterMaestro as cm
import numpy.linalg as npl
import os, time

def test_init_mono():
    molecule1 = cm.initilize.init.monomer()
    assert len(molecule1) == 66, "Correct number of atoms."
    assert len(molecule1.core) == 22, "Correct number of atoms in core."
    assert len(molecule1.substituents) == 4, "Correct number of substituents."
    assert len(molecule1.substituents[0]) == 11, "Correct number of atoms in substituents."

    molecule2 = cm.initilize.init.monomer(core="SnS", substituent="Np")
    assert len(molecule2) == 78, "Correct number of atoms."
    assert len(molecule2.core) == 10, "Correct number of atoms in core."
    assert len(molecule2.substituents) == 4, "Correct number of substituents."
    assert len(molecule2.substituents[0]) == 17, "Correct number of atoms in substituents."

def test_init_di():
    dimer1 = cm.initilize.init.dimer(subAngle=-30)
    assert len(dimer1) == 132, "Correct number of atoms"
    assert len(dimer1.molecules) == 2, "Correct number of molecules"
    assert len(dimer1.molecules[0]) == 66, "Correct number of atoms in molecule"
    assert dimer1.get_core_distance() > 5.99 and dimer1.get_core_distance() < 6.01, "dimer-core-distances are correct"

    dimer2 = cm.initilize.init.dimer(core="SnS", substituent="Np", subAngle=-30)
    assert len(dimer2) == 156, "Correct number of atoms"
    assert len(dimer2.molecules) == 2, "Correct number of molecules"
    assert len(dimer2.molecules[0]) == 78, "Correct number of atoms in molecule"
    assert dimer2.get_core_distance() > 5.99 and dimer2.get_core_distance() < 6.01, "dimer-core-distances are correct"

    dimer3 = cm.initilize.init.dimer(core="SnS", substituent="Np", subAngle=-30, rotAngle=60)
    assert len(dimer3) == 156, "Correct number of atoms"
    assert len(dimer3.molecules) == 2, "Correct number of molecules"
    assert len(dimer3.molecules[0]) == 78, "Correct number of atoms in molecule"
    assert dimer3.get_core_distance() > 5.99 and dimer3.get_core_distance() < 6.01, "dimer-core-distances are correct"

    dimer4 = cm.initilize.init.dimer(core="SiS", substituent="Ph", subAngle=-30, xshift=2)
    assert len(dimer4) == 108, "Correct number of atoms"
    assert len(dimer4.molecules) == 2, "Correct number of molecules"
    assert len(dimer4.molecules[0]) == 54, "Correct number of atoms in molecule"
    assert dimer4.get_core_distance() > 6.01, "dimer-core-distances are correct"


def test_init_solid():
    solid1 = cm.initilize.init.amorphSolid(core="Ad", substituent="Me", density=0.6, number=5)
    solid1.analyzeFAST(len(cm.initilize.init.monomer(core="Ad", substituent="Me")))
    density = solid1.calcDensity()
    assert density > 0.59 and density < 0.61, "Correct density"
    cm.io.write.poscar("cases/solid.vasp", solid1)
    cm.io.write.tmPeriodic("cases/solid.tm", solid1)
    time.sleep(0.01)
    os.remove("cases/solid.tm")

