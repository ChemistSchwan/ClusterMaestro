# Original source code proveded via BSD-2-clause license from:
# Copyright (c) 2013, Jimmy Charnley Kromann <jimmy@charnley.dk> & Lars Bratholm
# All rights reserved.
# Modifications are limited to the inclusion of the "combineOLDRMSD" and "combineAllRMSD" functions.

import numpy as np
import rmsd
from copy import deepcopy


def sortList(list_1, sorting_list):
    sorted_list = []
    for i in sorting_list:
        sorted_list.append(list_1[i])
    return np.array(sorted_list)


def combineOLDRMSD(molecule1, molecule2):
    rmsd_List = []
    rmsd_List.append(rmsd.kabsch_rmsd(molecule1.positions, molecule2.positions))  # 1
    moleculeOrder = rmsd.reorder_hungarian(
        np.array(molecule1.symbols),
        np.array(molecule2.symbols),
        molecule1.positions,
        molecule2.positions,
    )
    molecule2_Positions = deepcopy(molecule2.positions)
    molecule2_Positions = sortList(molecule2_Positions, moleculeOrder)
    rmsd_List.append(rmsd.kabsch_rmsd(molecule1.positions, molecule2_Positions))  # 2
    molecule_1 = deepcopy(molecule1.positions)
    molecule_2 = deepcopy(molecule2.positions)
    molecule_2 = rmsd.kabsch_rotate(molecule_2, molecule_1)
    moleculeOrder = rmsd.reorder_hungarian(
        np.array(molecule1.symbols), np.array(molecule2.symbols), molecule_1, molecule_2
    )
    molecule_2 = sortList(molecule_2, moleculeOrder)
    rmsd_List.append(rmsd.kabsch_rmsd(molecule_1, molecule_2))  # 3
    return min(rmsd_List)


def combineAllRMSD(
    molecule1,
    molecule2,
):
    """Different approaches for the calculation of the RMSD are performed and the smallest one is returned,
    as that one is assumed to be the most accurate.

    The calculation of the RMSD of teo molecules is surprisingly difficult, if the order of the atom AND the
    rotation shifts at the same time. Please keep in mind, that in these cases the return value might not be trusted.

    Parameters
    ----------

    molecule1: Atoms-object (or derivatives)
        First molecule to analyse.
    molecule2: Atoms-object (or derivative)
        Second molecule to analyse.
    """
    rmsdList = []
    rmsdList.append(simpleRMSD(molecule1, molecule2))
    rmsdList.append(kabschRMSD(molecule1, molecule2))
    rmsdList.append(hungarianRMSD(molecule1, molecule2))
    rmsdList.append(kabschHungRMSD(molecule1, molecule2))
    rmsdList.append(quarternionRMSD(molecule1, molecule2))
    rmsdList.append(reorderRMSD(molecule1, molecule2))
    rmsdList.append(reorderKabschRMSD(molecule1, molecule2))
    return min(rmsdList)


def simpleRMSD(molecule1, molecule2):
    """
    Returnes the RMSD of two molecules calculated using the simple method.

    Parameters
    ----------

    molecule1: Atoms-object (or derivatives)
        First molecule to analyse.
    molecule2: Atoms-object (or derivative)
        Second molecule to analyse.
    """
    return rmsd.rmsd(molecule1.positions, molecule2.positions)


def kabschRMSD(molecule1, molecule2):
    """
    Returnes the RMSD of two molecules calculated using the kabsch method for rotation.

    Parameters
    ----------

    molecule1: Atoms-object (or derivatives)
        First molecule to analyse.
    molecule2: Atoms-object (or derivative)
        Second molecule to analyse.
    """
    molecule_1 = deepcopy(molecule1.positions)
    molecule_2 = deepcopy(molecule2.positions)
    molecule_2 = rmsd.kabsch_rotate(molecule_2, molecule_1)
    return rmsd.kabsch_rmsd(molecule_1, molecule_2)


def hungarianRMSD(molecule1, molecule2):
    """
    Returnes the RMSD of two molecules calculated using the hungarian method for permutation.

    Parameters
    ----------

    molecule1: Atoms-object (or derivatives)
        First molecule to analyse.
    molecule2: Atoms-object (or derivative)
        Second molecule to analyse.
    """
    molecule_1 = deepcopy(molecule1.positions)
    molecule_2 = deepcopy(molecule2.positions)
    moleculeOrder = rmsd.reorder_hungarian(
        np.array(molecule1.symbols), np.array(molecule2.symbols), molecule_1, molecule_2
    )
    molecule_2 = sortList(molecule_2, moleculeOrder)
    return rmsd.rmsd(molecule_1, molecule_2)


def kabschHungRMSD(molecule1, molecule2):
    """
    Returnes the RMSD of two molecules calculated using a combination of the hungarian algorithm (permutation) and the kabsch-method
    (rotation). If the combination-approach is too time-intensive, THIS method is recommended, as it yields the best results in
    the tests.

    Parameters
    ----------

    molecule1: Atoms-object (or derivatives)
        First molecule to analyse.
    molecule2: Atoms-object (or derivative)
        Second molecule to analyse.
    """
    molecule_1 = deepcopy(molecule1.positions)
    molecule_2 = deepcopy(molecule2.positions)
    molecule_2 = rmsd.kabsch_rotate(molecule_2, molecule_1)
    moleculeOrder = rmsd.reorder_hungarian(
        np.array(molecule1.symbols), np.array(molecule2.symbols), molecule_1, molecule_2
    )
    molecule_2 = sortList(molecule_2, moleculeOrder)
    return rmsd.kabsch_rmsd(molecule_1, molecule_2)


def quarternionRMSD(molecule1, molecule2):
    """
    Returnes the RMSD of two molecules calculated using the quarternion method.

    Parameters
    ----------

    molecule1: Atoms-object (or derivatives)
        First molecule to analyse.
    molecule2: Atoms-object (or derivative)
        Second molecule to analyse.
    """
    molecule_1 = deepcopy(molecule1.positions)
    molecule_2 = deepcopy(molecule2.positions)
    return rmsd.quaternion_rmsd(molecule_1, molecule_2)


def reorderRMSD(molecule1, molecule2):
    """
    Returnes the RMSD of two molecules calculated using the reorder_distance method.

    Parameters
    ----------

    molecule1: Atoms-object (or derivatives)
        First molecule to analyse.
    molecule2: Atoms-object (or derivative)
        Second molecule to analyse.
    """
    molecule_1 = deepcopy(molecule1.positions)
    molecule_2 = deepcopy(molecule2.positions)
    mol_1POS = deepcopy(molecule1.positions)
    mol_2POS = deepcopy(molecule2.positions)
    mol_1sym = np.array(deepcopy(molecule1.symbols))
    mol_2sym = np.array(deepcopy(molecule2.symbols))
    moleculeOrder = rmsd.reorder_distance(mol_1sym, mol_2sym, mol_1POS, mol_2POS)
    molecule_2 = sortList(molecule_2, moleculeOrder)
    return rmsd.rmsd(molecule_1, molecule_2)


def reorderKabschRMSD(molecule1, molecule2):
    """
    Returnes the RMSD of two molecules calculated using a combination of the reorder_distance and kabsch-method.

    Parameters
    ----------

    molecule1: Atoms-object (or derivatives)
        First molecule to analyse.
    molecule2: Atoms-object (or derivative)
        Second molecule to analyse.
    """
    molecule_1 = deepcopy(molecule1.positions)
    molecule_2 = deepcopy(molecule2.positions)
    mol_1POS = deepcopy(molecule1.positions)
    mol_2POS = deepcopy(molecule2.positions)
    mol_1sym = np.array(deepcopy(molecule1.symbols))
    mol_2sym = np.array(deepcopy(molecule2.symbols))
    moleculeOrder = rmsd.reorder_distance(mol_1sym, mol_2sym, mol_1POS, mol_2POS)
    molecule_2 = sortList(molecule_2, moleculeOrder)
    return rmsd.kabsch_rmsd(molecule_1, molecule_2)
