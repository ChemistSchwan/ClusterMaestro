import ClusterMaestro as cm
import pymatgen as pmg
import numpy as np
from ase.io import read
import os, shutil, time


def symmetryScanMono(core="Ad", substituent="Ph"):
    system = cm.initilize.init.monomer(core=core, substituent=substituent)
    for i in system.substituents:
        i.get_connection(system.core)
        i.get_dihedral(system.core)
        i.torsion(0 - i.dihedral)
        i.get_connection(system.core)
    system._renew_positions(structure=["all"])

    if os.path.exists("symScanstructs"):
        shutil.rmtree("symScanstructs")
    os.mkdir("symScanstructs")

    symmetryList = []
    molecule = pmg.core.structure.Molecule(
        species=system.symbols, coords=system.positions
    )
    symmetry = pmg.symmetry.analyzer.PointGroupAnalyzer(molecule).sch_symbol
    symmetryList.append(symmetry)
    system.write("symScanstructs/%s.xyz" % symmetry)
    counter = 0
    for val0, i in enumerate(range(0, 180, 30)):
        system.torsion(system.substituents[0], i)
        for val1, j in enumerate(range(0, 180, 30)):
            system.torsion(system.substituents[1], j)
            for val2, k in enumerate(range(0, 180, 30)):
                system.torsion(system.substituents[2], k)
                for val3, l in enumerate(range(0, 180, 30)):
                    system.torsion(system.substituents[3], l)
                    molecule = pmg.core.structure.Molecule(
                        species=system.symbols, coords=system.positions
                    )
                    symmetry = pmg.symmetry.analyzer.PointGroupAnalyzer(
                        molecule
                    ).sch_symbol
                    print(symmetry)
                    if symmetry not in symmetryList:
                        symmetryList.append(symmetry)
                        system.write("symScanstructs/%s.xyz" % symmetry)
                    counter += 1
    system.write("test.xyz")
    print(counter)
