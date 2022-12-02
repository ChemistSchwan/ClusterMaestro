import ClusterMaestro as cm
from ase.io import read, write
import rmsd, time, os
import numpy as np


def test_rmsdMethods():
    dimer1 = cm.System.Oligomer(read("cases/rmsd_0.xyz", format="xyz"))
    dimerPerm = cm.System.Oligomer(read("cases/rmsd_1.xyz", format="xyz"))
    dimerRot = cm.System.Oligomer(read("cases/rmsd_2.xyz", format="xyz"))
    dimerRotPerm = cm.System.Oligomer(read("cases/rmsd_3.xyz", format="xyz"))

    dimerList = [dimer1, dimerPerm, dimerRot, dimerRotPerm]
    labelList = ["dimer", "Perm", "Rot", "PermRot"]

    result = np.zeros((4,4), dtype="float64")
    print("RMSD using the simple RMSD method")
    print("%20s%10s%10s%10s" % ("dimer", "Perm", "Rot", "PermRot"))
    start = time.time()
    for i, dim1 in enumerate(dimerList):
        for j, dim2 in enumerate(dimerList):
            result[i,j] =  cm.util.rmsd.simpleRMSD(dim1, dim2)
    for i, vals in enumerate(result):
        print("%10s%10.3f%10.3f%10.3f%10.3f" % (labelList[i], vals[0], vals[1], vals[2], vals[3]))
    print("Time needed for total execution: %.5f\n" % (time.time() - start))


    result = np.zeros((4,4), dtype="float64")
    print("RMSD using the KABSCH RMSD method")
    print("%20s%10s%10s%10s" % ("dimer", "Perm", "Rot", "PermRot"))
    start = time.time()
    for i, dim1 in enumerate(dimerList):
        for j, dim2 in enumerate(dimerList):
            result[i,j] =  cm.util.rmsd.kabschRMSD(dim1, dim2)
    for i, vals in enumerate(result):
        print("%10s%10.3f%10.3f%10.3f%10.3f" % (labelList[i], vals[0], vals[1], vals[2], vals[3]))
    print("Time needed for total execution: %.5f\n" % (time.time() - start))

    result = np.zeros((4,4), dtype="float64")
    print("RMSD using the HUNGARIAN RMSD method")
    print("%20s%10s%10s%10s" % ("dimer", "Perm", "Rot", "PermRot"))
    start = time.time()
    for i, dim1 in enumerate(dimerList):
        for j, dim2 in enumerate(dimerList):
            result[i,j] =  cm.util.rmsd.hungarianRMSD(dim1, dim2)
    for i, vals in enumerate(result):
        print("%10s%10.3f%10.3f%10.3f%10.3f" % (labelList[i], vals[0], vals[1], vals[2], vals[3]))
    print("Time needed for total execution: %.5f\n" % (time.time() - start))


    result = np.zeros((4,4), dtype="float64")
    print("RMSD using the combined KABSCH and HUNGARIAN RMSD method")
    print("%20s%10s%10s%10s" % ("dimer", "Perm", "Rot", "PermRot"))
    start = time.time()
    for i, dim1 in enumerate(dimerList):
        for j, dim2 in enumerate(dimerList):
            result[i,j] =  cm.util.rmsd.kabschHungRMSD(dim1, dim2)
            assert result[i,j] < 0.01, "Kabsch-hungarian is right"
    for i, vals in enumerate(result):
        print("%10s%10.3f%10.3f%10.3f%10.3f" % (labelList[i], vals[0], vals[1], vals[2], vals[3]))
    print("Time needed for total execution: %.5f\n" % (time.time() - start))


    result = np.zeros((4,4), dtype="float64")
    print("RMSD using the QUARTERNIONEN RMSD method")
    print("%20s%10s%10s%10s" % ("dimer", "Perm", "Rot", "PermRot"))
    start = time.time()
    for i, dim1 in enumerate(dimerList):
        for j, dim2 in enumerate(dimerList):
            result[i,j] =  cm.util.rmsd.quarternionRMSD(dim1, dim2)
    for i, vals in enumerate(result):
        print("%10s%10.3f%10.3f%10.3f%10.3f" % (labelList[i], vals[0], vals[1], vals[2], vals[3]))
    print("Time needed for total execution: %.5f\n" % (time.time() - start))


    result = np.zeros((4,4), dtype="float64")
    print("RMSD using the REORDER_DISTANCE RMSD method")
    print("%20s%10s%10s%10s" % ("dimer", "Perm", "Rot", "PermRot"))
    start = time.time()
    for i, dim1 in enumerate(dimerList):
        for j, dim2 in enumerate(dimerList):
            result[i,j] =  cm.util.rmsd.reorderRMSD(dim1, dim2)
    for i, vals in enumerate(result):
        print("%10s%10.3f%10.3f%10.3f%10.3f" % (labelList[i], vals[0], vals[1], vals[2], vals[3]))
    print("Time needed for total execution: %.5f\n" % (time.time() - start))


    result = np.zeros((4,4), dtype="float64")
    print("RMSD using the combined REORDER_DISTANCE-KABSCH RMSD method")
    print("%20s%10s%10s%10s" % ("dimer", "Perm", "Rot", "PermRot"))
    start = time.time()
    for i, dim1 in enumerate(dimerList):
        for j, dim2 in enumerate(dimerList):
            result[i,j] =  cm.util.rmsd.reorderKabschRMSD(dim1, dim2)
    for i, vals in enumerate(result):
        print("%10s%10.3f%10.3f%10.3f%10.3f" % (labelList[i], vals[0], vals[1], vals[2], vals[3]))
    print("Time needed for total execution: %.5f\n" % (time.time() - start))


    result = np.zeros((4,4), dtype="float64")
    print("RMSD using the COMBINED RMSD method")
    print("%20s%10s%10s%10s" % ("dimer", "Perm", "Rot", "PermRot"))
    start = time.time()
    for i, dim1 in enumerate(dimerList):
        for j, dim2 in enumerate(dimerList):
            result[i,j] =  cm.util.rmsd.combineAllRMSD(dim1, dim2)
            assert result[i,j] < 0.01, "Combination approach is right"
    for i, vals in enumerate(result):
        print("%10s%10.3f%10.3f%10.3f%10.3f" % (labelList[i], vals[0], vals[1], vals[2], vals[3]))
    print("Time needed for total execution: %.5f\n" % (time.time() - start))


    result = np.zeros((4,4), dtype="float64")
    print("RMSD using the OLD COMBINED RMSD method")
    print("%20s%10s%10s%10s" % ("dimer", "Perm", "Rot", "PermRot"))
    start = time.time()
    for i, dim1 in enumerate(dimerList):
        for j, dim2 in enumerate(dimerList):
            result[i,j] =  cm.util.rmsd.combineOLDRMSD(dim1, dim2)
    for i, vals in enumerate(result):
        print("%10s%10.3f%10.3f%10.3f%10.3f" % (labelList[i], vals[0], vals[1], vals[2], vals[3]))
    print("Time needed for total execution: %.5f\n" % (time.time() - start))
