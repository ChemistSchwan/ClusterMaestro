import numpy as np
from copy import deepcopy
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

def rmsd(V, W):
    D = len(V[0])
    N = len(V)
    result = 0.0
    for v, w in zip(V, W):
        result += sum([(v[i] - w[i])**2.0 for i in range(D)])
    return np.sqrt(result/N)

def kabsch_rmsd(P, Q, translate=False):
    if translate:
        Q = Q - centroid(Q)
        P = P - centroid(P)
    P = kabsch_rotate(P, Q)
    return rmsd(P, Q)

def kabsch_rotate(P, Q):
    U = kabsch(P, Q)
    P = np.dot(P, U)
    return P

def kabsch(P, Q):
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = np.dot(V, W)
    return U

def quaternion_rmsd(P, Q):
    rot = quaternion_rotate(P, Q)
    P = np.dot(P, rot)
    return rmsd(P, Q)

def quaternion_transform(r):
    Wt_r = makeW(*r).T
    Q_r = makeQ(*r)
    rot = Wt_r.dot(Q_r)[:3, :3]
    return rot

def quaternion_rotate(X, Y):
    N = X.shape[0]
    W = np.asarray([makeW(*Y[k]) for k in range(N)])
    Q = np.asarray([makeQ(*X[k]) for k in range(N)])
    Qt_dot_W = np.asarray([np.dot(Q[k].T, W[k]) for k in range(N)])
    W_minus_Q = np.asarray([W[k] - Q[k] for k in range(N)])
    A = np.sum(Qt_dot_W, axis=0)
    eigen = np.linalg.eigh(A)
    r = eigen[1][:, eigen[0].argmax()]
    rot = quaternion_transform(r)
    return rot

def makeW(r1, r2, r3, r4=0):
    W = np.asarray([
        [r4, r3, -r2, r1],
        [-r3, r4, r1, r2],
        [r2, -r1, r4, r3],
        [-r1, -r2, -r3, r4]])
    return W

def makeQ(r1, r2, r3, r4=0):
    Q = np.asarray([
        [r4, -r3, r2, r1],
        [r3, r4, -r1, r2],
        [-r2, r1, r4, r3],
        [-r1, -r2, -r3, r4]])
    return Q

def hungarian(A, B):
    distances = cdist(A, B, 'euclidean')
    indices_a, indices_b = linear_sum_assignment(distances)
    return indices_b

def reorder_hungarian(p_atoms, q_atoms, p_coord, q_coord):
    unique_atoms = np.unique(p_atoms)
    view_reorder = np.zeros(q_atoms.shape, dtype=int)
    view_reorder -= 1
    for atom in unique_atoms:
        p_atom_idx, = np.where(p_atoms == atom)
        q_atom_idx, = np.where(q_atoms == atom)
        A_coord = p_coord[p_atom_idx]
        B_coord = q_coord[q_atom_idx]
        view = hungarian(A_coord, B_coord)
        view_reorder[p_atom_idx] = q_atom_idx[view]
    return view_reorder

def reorder_distance(p_atoms, q_atoms, p_coord, q_coord):
    unique_atoms = np.unique(p_atoms)
    view_reorder = np.zeros(q_atoms.shape, dtype=int)
    for atom in unique_atoms:
        p_atom_idx, = np.where(p_atoms == atom)
        q_atom_idx, = np.where(q_atoms == atom)
        A_coord = p_coord[p_atom_idx]
        B_coord = q_coord[q_atom_idx]
        A_norms = np.linalg.norm(A_coord, axis=1)
        B_norms = np.linalg.norm(B_coord, axis=1)
        reorder_indices_A = np.argsort(A_norms)
        reorder_indices_B = np.argsort(B_norms)
        translator = np.argsort(reorder_indices_A)
        view = reorder_indices_B[translator]
        view_reorder[p_atom_idx] = q_atom_idx[view]
    return view_reorder

def centroid(X):
    C = X.mean(axis=0)
    return C

def sortList(list_1, sorting_list):
    sorted_list = []
    for i in sorting_list:
        sorted_list.append(list_1[i])
    return np.array(sorted_list)

def combineOLDRMSD(molecule1, molecule2):
    rmsd_List = []
    rmsd_List.append(kabsch_rmsd(molecule1.positions, molecule2.positions))            #1
    moleculeOrder = reorder_hungarian(np.array(molecule1.symbols), np.array(molecule2.symbols), molecule1.positions, molecule2.positions)
    molecule2_Positions = deepcopy(molecule2.positions)
    molecule2_Positions = sortList(molecule2_Positions, moleculeOrder)
    rmsd_List.append(kabsch_rmsd(molecule1.positions, molecule2_Positions))            #2
    molecule_1 = deepcopy(molecule1.positions)
    molecule_2 = deepcopy(molecule2.positions)
    molecule_2 = kabsch_rotate(molecule_2, molecule_1)
    moleculeOrder = reorder_hungarian(np.array(molecule1.symbols), np.array(molecule2.symbols), molecule_1, molecule_2)
    molecule_2 = sortList(molecule_2, moleculeOrder)
    rmsd_List.append(kabsch_rmsd(molecule_1, molecule_2))                              #3
    return min(rmsd_List)



def combineAllRMSD(molecule1, molecule2, ):
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
    return rmsd(molecule1.positions, molecule2.positions)

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
    molecule_2 = kabsch_rotate(molecule_2, molecule_1)
    return kabsch_rmsd(molecule_1, molecule_2)

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
    moleculeOrder = reorder_hungarian(np.array(molecule1.symbols), np.array(molecule2.symbols), molecule_1, molecule_2)
    molecule_2 = sortList(molecule_2, moleculeOrder)
    return rmsd(molecule_1, molecule_2)

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
    molecule_2 = kabsch_rotate(molecule_2, molecule_1)
    moleculeOrder = reorder_hungarian(np.array(molecule1.symbols), np.array(molecule2.symbols), molecule_1, molecule_2)
    molecule_2 = sortList(molecule_2, moleculeOrder)
    return kabsch_rmsd(molecule_1, molecule_2)

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
    return quaternion_rmsd(molecule_1, molecule_2)

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
    moleculeOrder = reorder_distance(mol_1sym, mol_2sym, mol_1POS, mol_2POS)
    molecule_2 = sortList(molecule_2, moleculeOrder)
    return rmsd(molecule_1, molecule_2)

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
    moleculeOrder = reorder_distance(mol_1sym, mol_2sym, mol_1POS, mol_2POS)
    molecule_2 = sortList(molecule_2, moleculeOrder)
    return kabsch_rmsd(molecule_1, molecule_2)

