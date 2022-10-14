import numpy as np
from numpy.linalg import norm
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from copy import copy, deepcopy
import networkx as nx
from networkx.algorithms import isomorphism
import numpy as np
import multiprocessing
import glob
import ntpath
import sys

from ..utils import ELEMENT_NAMES, NAMES_ELEMENT, H2KC
from .. import utils

RADII = {'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.69, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.7, 'Ti': 1.6, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 'Co': 1.5, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.2, 'As': 1.19, 'Se': 1.2, 'Br': 1.2, 'Kr': 1.16, 'Rb': 2.2, 'Sr': 1.95, 'Y': 1.9, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.4, 'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.9, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.7, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.4, 'At': 1.5, 'Rn': 1.5, 'Fr': 2.6, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 'Pa': 2.0, 'U': 1.96, 'Np': 1.9, 'Pu': 1.87, 'Am': 1.8, 'Cm': 1.69}

def same_element(n1_attrib, n2_attrib):
    return n1_attrib['symbol'] == n2_attrib['symbol']

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
    
    def distance_count(self, first_atom, second_atom, crit):
        count = 0
        for struct in self.structures:
            dist = utils.get_length((first_atom, second_atom), struct)
            if crit(dist):
                count += 1
        return count
            
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
                comment = "%14.9f" % self.energies[i]
            else:
                comment = ""
            parts.append(utils.to_xyz(item, self.atom_symbols, description=comment))
        
        with open(filename, 'w') as f:
            f.write('\n'.join(parts))

    def write_pair(self, rmsd, p_coord, q_coord, i, j):
        parts = []
        for idx, coord in [[i, p_coord], [j, q_coord]]:
            comment = "idx = {}; rmsd = {}".format(idx, rmsd)
            parts.append(utils.to_xyz(coord, self.atom_symbols, description=comment))
        
        for idx in [i, j]:
            coord = self.structures[idx]
            comment = "idx = {}; rmsd = {}".format(idx, rmsd)
            parts.append(utils.to_xyz(coord, self.atom_symbols, description=comment))
        
        with open(f"pair_{i}_{j}.xyz", 'w') as f:
            f.write('\n'.join(parts))

    def deduce_topology(self, index, filename='check_topology', mult=1.0):
        self.G = nx.Graph()
        for i in range(len(self.atom_symbols)):
            self.G.add_node(i)
            self.G.nodes[i]['symbol'] = self.atom_symbols[i]
        for nodeA in range(len(self.atom_symbols)):
            for nodeB in range(nodeA):
                max_dist = mult * (RADII[self.atom_symbols[nodeA]] + RADII[self.atom_symbols[nodeB]])
                if norm(self.structures[index][nodeA] - self.structures[index][nodeB]) < max_dist:
                    self.G.add_edge(nodeA, nodeB)
        
        lines = ["", "", ""]
        lines.append("%3d%3d  0  0  0  0  0  0  0  0999 V2000" % (self.G.number_of_nodes(), self.G.number_of_edges()))
        #   -5.5250    1.6470    1.0014 C   0  0  0  0  0  0  0  0  0  0  0  0
        for atom in range(self.G.number_of_nodes()):
            lines.append("%10.4f%10.4f%10.4f%3s  0  0  0  0  0  0  0  0  0  0  0  0" % (
                self.structures[index][atom][0],
                self.structures[index][atom][1],
                self.structures[index][atom][2],
                self.atom_symbols[atom]))

        for edge in self.G.edges:
            lines.append("%3s%3s%3s  0" % (edge[0] + 1,
                                           edge[1] + 1,
                                           1))
        lines.append("M  END\n")

        with open(f'{filename}.sdf', "w") as f:
            f.write("\n".join(lines))

    def filter_rmsd(self, maxrmsd, energy_thr=None, ignore_elements=[]):
        save_atoms = []
        ignores_atoms = []
        for node in self.G.nodes:
            if self.G.nodes[node]['symbol'] not in ignore_elements:
                save_atoms.append(node)
            else:
                ignores_atoms.append(node)
        mysubgraph = self.G.subgraph(save_atoms)
        self.GM = isomorphism.GraphMatcher(mysubgraph, mysubgraph, node_match=same_element)

        isomorphisms = []
        for isom in self.GM.isomorphisms_iter():
            isomorphisms.append(np.array([isom[i] for i in range(mysubgraph.number_of_nodes())]))
        simple_reorder = list(mysubgraph.nodes)
        print(repr(isomorphisms))
        i = len(self.energies) - 1
        while i > 0:
            print("i = " + str(i))
            curgeom = self.structures[i]
            j = i - 1
            while j >= 0:
                if energy_thr is None or abs(self.energies[i] - self.energies[j]) < energy_thr:
                    testgeom = self.structures[j]
                    minrmsd = None
                    best_isom = None
                    for n_iso, isom in enumerate(isomorphisms):
                        rmsd, p_coord, q_coord = RmsdPool.calc_rmsd_new(curgeom, testgeom, isom, simple_reorder, self.atom_ints)
                        # print(f"RMSD {i} - {j} ({n_iso})= {rmsd}")
                        if minrmsd is None or rmsd < minrmsd:
                            minrmsd = rmsd
                            best_isom = n_iso

                    # rmsd, p_coord, q_coord = RmsdPool.calc_rmsd_new(curgeom, testgeom, isomorphisms[best_isom], simple_reorder, self.atom_ints)
                    # self.write_pair(rmsd, p_coord, q_coord, i, j)
                    if rmsd < maxrmsd:
                        break
                j -= 1
            
            if j != -1:
                print("Match detected. RMSD = " + str(rmsd))
                del self.structures[i]
                del self.energies[i]
                del self.comments[i]
            i -= 1

    def energy_sort(self):
        self.energies, self.structures, self.comments = map(list, zip(*sorted(zip(self.energies, self.structures, self.comments), key=lambda attrs: attrs[0])))

    @staticmethod
    def calc_rmsd_new(p_all, q_all, view_reorder, simple_reorder, atom_syms):
        p_coord = deepcopy(p_all[simple_reorder])
        q_coord = deepcopy(q_all[view_reorder])

        p_cent = p_coord.mean(axis=0)
        q_cent = q_coord.mean(axis=0)
        p_coord -= p_cent
        q_coord -= q_cent

        C = np.dot(np.transpose(p_coord), q_coord)
        V, S, W = np.linalg.svd(C)
        # d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
        # if d:
        #     S[-1] = -S[-1]
        #     V[:, -1] = -V[:, -1]
        U = np.dot(V, W)
        p_coord = np.dot(p_coord, U)
        diff = np.array(p_coord) - np.array(q_coord)
        N = len(p_coord)
        res = np.sqrt((diff * diff).sum() / N)
        return res, p_coord, q_coord

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
        # if d:
        #     S[-1] = -S[-1]
        #     V[:, -1] = -V[:, -1]
        U = np.dot(V, W)
        p_coord = np.dot(p_coord, U)
        diff = np.array(p_coord) - np.array(q_coord)
        N = len(p_coord)
        res = np.sqrt((diff * diff).sum() / N)
        return res, p_coord, q_coord
