from ClusterMaster.initilize import init
import ClusterMaster as cm

solid = init.amorphSolid(number=10, density=1)

cm.util.write.poscar("solid.vasp", solid, box=[solid.cell[0,0], solid.cell[1,1], solid.cell[2,2]])

dimer = init.dimer()

dimer.write("dimer.xyz", format="xyz")

mono = init.monomer()

mono.write("mono.xyz", format="xyz")
