import ClusterMaestro as cm
import os

def test_core():
    cm.initilize.expand.core("cases/benzol.xyz", "benzol", delete=True, level="z")
    benzol = cm.initilize.init.dimer(core="benzol", substituent="Ph", subAngle=0)
    benzol = cm.initilize.init.monomer(core="benzol", substituent="Ph")

    benzol.core.write("test.xyz")
    benzol.write("full.xyz")
    print(len(benzol.substituents))
    print(benzol.core)
    for i in range(len(benzol.substituents)):
        print(len(benzol.substituents[i]))
        benzol.substituents[i].write("%d.xyz" % i)
    print(benzol)
    ref_1 = cm.initilize.init.dimer(core="SnTe", substituent="Me", subAngle=0)
    ref = cm.initilize.init.dimer(core="Ad", substituent="Ph", subAngle=0)
    for file in ["0.xyz", "1.xyz", "2.xyz", "full.xyz", "test.xyz"]: os.remove(file)

