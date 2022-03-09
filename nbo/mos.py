import numpy as np

def is_occupied(dm_elem):
    if dm_elem > 1.9:
        return True
    elif dm_elem < 0.1:
        return False
    else:
        raise Exception("Unusial MO occupancy!")

def get_boundary_indices(dmmo):
    occ = np.diag(dmmo)
    for i in range(len(occ)):
        if is_occupied(occ[i]) and not is_occupied(occ[i + 1]):
            return i, i + 1
