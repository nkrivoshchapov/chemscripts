import numpy as np
from copy import deepcopy


def align_molecules(main_xyzs, ref_xyzs):
    main_np = np.array(main_xyzs)
    ref_np = np.array(ref_xyzs)
    assert main_np.shape[1] == 3 and ref_np.shape[1] == 3
    assert main_np.shape[0] == ref_np.shape[0]
    _, rot, trans = calc_rmsd(ref_np, main_np, np.array(['A'] * ref_np.shape[0]))
    return rot, trans

def calc_rmsd(p_all, q_all, atom_syms, ignore_atom=None, symmetry=False):  # All three args are numpy arrays
    """
    Here we use implementation of Kabsch algorithm with Hungarian reordering from https://github.com/charnley/rmsd
    """
    if ignore_atom is not None:
        myview = np.where(atom_syms != ignore_atom)
        p_coord = deepcopy(p_all[myview])
        q_coord = deepcopy(q_all[myview])
        atoms_view = deepcopy(atom_syms[myview])
    else:
        p_coord = deepcopy(p_all)
        q_coord = deepcopy(q_all)
        atoms_view = deepcopy(atom_syms)

    p_cent = p_coord.mean(axis=0)
    print("Start: " + repr(p_cent))
    q_cent = q_coord.mean(axis=0)
    print("End: " + repr(q_cent))
    p_coord -= p_cent
    q_coord -= q_cent

    if symmetry:
        from scipy.optimize import linear_sum_assignment
        from scipy.spatial.distance import cdist
        unique_atoms = np.unique(atoms_view)
        view_reorder = np.zeros(atoms_view.shape, dtype=int)
        view_reorder -= 1
        for atom in unique_atoms:
            (atom_idx,) = np.where(atoms_view == atom)
            A_coord = p_coord[atom_idx]
            B_coord = q_coord[atom_idx]

            distances = cdist(A_coord, B_coord, "euclidean")
            indices_a, indices_b = linear_sum_assignment(distances)
            view_reorder[atom_idx] = atom_idx[indices_b]
        q_coord = q_coord[view_reorder]

    C = np.dot(np.transpose(p_coord), q_coord)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = np.dot(V, W)
    p_coord = np.dot(p_coord, U)
    diff = np.array(p_coord) - np.array(q_coord)
    N = len(p_coord)
    rmsd = np.sqrt((diff * diff).sum() / N)
    q_new_cent = U @ q_cent
    return rmsd, U, p_cent - q_new_cent
