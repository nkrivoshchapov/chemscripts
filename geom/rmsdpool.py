from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from copy import copy, deepcopy
import numpy as np
import multiprocessing
import glob
import ntpath
import sys

from ..utils import ELEMENT_NAMES, NAMES_ELEMENT, H2KC
from .. import utils

class RmsdPool: # TODO Implement sorting of this pool
    def __init__(self):
        self.structures = []
        self.energies = []
        self.comments = []
        self.atom_ints = None
        self.atom_symbols = None
    
    @staticmethod
    def struct_from_lines(lines):
        struct = []
        for line in lines:
            parts = line.split()
            struct.append([float(parts[1]), float(parts[2]), float(parts[3])])
        return np.array(struct)
    
    @staticmethod
    def symbols_from_lines(lines):
        res = []
        for line in lines:
            res.append(line.split()[0])
        return res
    
    @staticmethod
    def atom_to_int(atom_sym):
        return NAMES_ELEMENT[atom_sym.capitalize()]
    
    def _add_structure(self, xyz, ener, comment):
        self.structures.append(xyz)
        self.energies.append(ener)
        self.comments.append(comment)
    
    def include_structure(self, xyz, sym, energy=None, comment=None):
        if self.atom_symbols is not None:
            assert self.atom_symbols == sym
        else:
            self.atom_symbols = sym
            atom_ints = []
            for atom in self.atom_symbols:
                atom_ints.append(RmsdPool.atom_to_int(atom))
            self.atom_ints = np.array(atom_ints)
        self._add_structure(xyz, energy, comment)
    
    def include(self, xyzfile, energy_func=None, save_comment=False):
        lines = open(xyzfile, 'r').readlines()
        startline = lines[0]
        
        delimiters = []
        for i, line in enumerate(lines):
            if line == startline:
                delimiters.append(i)
        delimiters.append(len(lines))
        
        symbols = RmsdPool.symbols_from_lines(lines[2:delimiters[1]])
        if self.atom_symbols is not None:
            assert self.atom_symbols == symbols
        
        # minener = float(lines[1].replace('\n', ''))
        for i, idx in enumerate(delimiters[:len(delimiters) - 1]):
            comment = lines[idx + 1].replace('\n', '')
            if energy_func is not None:
                energy = energy_func(comment)
            else:
                energy = None
            
            if save_comment:
                mycomment = comment
            else:
                mycomment = None
            
            structure = RmsdPool.struct_from_lines(lines[idx + 2:delimiters[i + 1]])
            self._add_structure(structure, energy, comment=mycomment)
            assert symbols == RmsdPool.symbols_from_lines(lines[idx + 2:delimiters[i + 1]])
        
        if self.atom_symbols is None:
            self.atom_symbols = symbols
            atom_ints = []
            for atom in self.atom_symbols:
                atom_ints.append(RmsdPool.atom_to_int(atom))
            self.atom_ints = np.array(atom_ints)
    
    def energy_filter(self, maxener, etype='kcal'): # 'au' or 'kcal'
        minener = None
        for item in self.energies:
            if minener is None or minener > item:
                minener = item
        
        if etype == 'au':
            keep = lambda de: de < maxener
        elif etype == 'kcal':
            keep = lambda de: de * H2KC < maxener
        
        i = len(self.energies) - 1
        while i >= 0:
            # print("i = " + str(i))
            curener = self.energies[i]
            if not keep(curener - minener):
                # print("Removed")
                del self.structures[i]
                del self.energies[i]
                del self.comments[i]
            i -= 1
    
    def distance_filter(self, first_atom, second_atom, keep_if):
        i = len(self.energies) - 1
        while i >= 0:
            # print("i = " + str(i))
            dist = utils.get_length((first_atom, second_atom), self.structures[i])
            if not keep_if(dist):
                # print("Removed")
                del self.structures[i]
                del self.energies[i]
                del self.comments[i]
            i -= 1
            
    def dihedral_filter(self, indices, keep_if):
        i = len(self.energies) - 1
        while i >= 0:
            dihedral = utils.get_dihedral(atoms=indices, xyz=self.structures[i])
            if not keep_if(dihedral):
                del self.structures[i]
                del self.energies[i]
                del self.comments[i]
            i -= 1
    
    def save(self, filename):
        parts = []
        for i, item in enumerate(self.structures):
            if self.comments[i] is not None:
                comment = self.comments[i]
            elif self.energies[i] is not None:
                comment = "%14.6f"%self.energies[i]
            else:
                comment = ""
            parts.append(utils.to_xyz(item, self.atom_symbols, description=comment))
        
        with open(filename, 'w') as f:
            f.write('\n'.join(parts))

    def filter_rmsd(self, maxrmsd):
        i = len(self.energies) - 1
        while i > 0:
            print("i = " + str(i))
            curgeom = self.structures[i]
            j = i - 1
            while j >= 0:
                testgeom = self.structures[j]
                rmsd = RmsdPool.calc_rmsd(curgeom, testgeom, self.atom_ints)
                if rmsd < maxrmsd:
                    break
                j -= 1
            
            if j != -1:
                print("Match detected. RMSD = " + str(rmsd))
                del self.structures[i]
                del self.energies[i]
                del self.comments[i]
            i -= 1

    @staticmethod
    def calc_rmsd(p_all, q_all, atom_syms, ignore_atom=None):  # All three args are numpy arrays
        """
        Here we use implementation of Kabsch algorithm with Hungarian reordering from https://github.com/charnley/rmsd
        """
        if ignore_atom is not None: # TODO ignore_atoms from config
            myview = np.where(atom_syms != ignore_atom)
            p_coord = deepcopy(p_all[myview])
            q_coord = deepcopy(q_all[myview])
            atoms_view = deepcopy(atom_syms[myview])
        else:
            p_coord = deepcopy(p_all)
            q_coord = deepcopy(q_all)
            atoms_view = deepcopy(atom_syms)

        p_cent = p_coord.mean(axis=0)
        q_cent = q_coord.mean(axis=0)
        p_coord -= p_cent
        q_coord -= q_cent

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
        res = np.sqrt((diff * diff).sum() / N)
        return res
    