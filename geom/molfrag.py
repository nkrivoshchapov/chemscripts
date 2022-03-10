import networkx as nx
import numpy as np
from numpy.linalg import inv, norm
import copy

from ..nbo import NBO3LogParser, NBO6LogParser

class Molecule:
    def __init__(self, sdf=None, nbolog=None):
        if sdf is not None:
            self.readsdf(sdf)
        elif nbolog is not None:
            self.readnbotopology(nbolog) # Must read symbols and xyzs by calling 'from_xyz'
        else:
            raise Exception("No args were given")

    def readsdf(self, file):
        lines = open(file, "r").readlines()
        natoms = int(lines[3][0:3])
        nbonds = int(lines[3][3:6])
        self.G = nx.Graph()
        for i in range(4, 4 + natoms):
            self.G.add_node(i-4)
            parts = lines[i].replace("\n", "").split()
            self.G.nodes[i-4]['xyz'] = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
            self.G.nodes[i-4]['symbol'] = parts[3]
        for i in range(4 + natoms, 4 + natoms + nbonds):
            at1 = int(lines[i][0:3])
            at2 = int(lines[i][3:7])
            bondtype = int(lines[i][7:10])
            self.G.add_edge(at1 - 1, at2 - 1)
            self.G[at1 - 1][at2 - 1]['type'] = bondtype

    def readnbotopology(self, nbolog):
        if nbolog.endswith('.log'):
            parser = NBO3LogParser(nbolog)
        elif nbolog.endswith('.out'):
            parser = NBO6LogParser(nbolog)
        else:
            raise Exception(NotImplementedError)
        bonds = parser.get_bondlist()
        self.G = nx.Graph() # TODO What to do with unbonded atoms??
        for bond in bonds:
            self.G.add_edge(bond[0] - 1, bond[1] - 1)
            self.G[bond[0] - 1][bond[1] - 1]['type'] = 1

    def from_xyz(self, xyzs, syms):
        for i, xyz in enumerate(xyzs):
            self.G.nodes[i]['xyz'] = copy.deepcopy(xyz)
            if 'symbol' in self.G.nodes[i]:
                assert self.G.nodes[i]['symbol'] == syms[i]
            else:
                self.G.nodes[i]['symbol'] = syms[i]

    def as_xyz(self):
        xyzs = []
        syms = []
        for atom in range(self.G.number_of_nodes()):
            syms.append(self.G.nodes[atom]['symbol'])
            xyzs.append(self.G.nodes[atom]['xyz'])
        return xyzs, syms

    def remove_bond(self, a, b):
        self.G.remove_edge(a - 1, b - 1)

    def has_bond(self, a, b):
        return self.G.has_edge(a - 1, b - 1)

    def xyz(self, i):
        return self.G.nodes[i - 1]['xyz']

    def save_sdf(self, sdfname):
        lines = ["", "", ""]
        lines.append("%3d%3d  0  0  0  0  0  0  0  0999 V2000" % (self.G.number_of_nodes(), self.G.number_of_nodes()))
        #   -5.5250    1.6470    1.0014 C   0  0  0  0  0  0  0  0  0  0  0  0
        for atom in range(self.G.number_of_nodes()):
            lines.append("%10.4f%10.4f%10.4f%3s  0  0  0  0  0  0  0  0  0  0  0  0" % (
                self.G.nodes[atom]['xyz'][0],
                self.G.nodes[atom]['xyz'][1],
                self.G.nodes[atom]['xyz'][2],
                self.G.nodes[atom]['symbol']))

        for edge in self.G.edges:
            lines.append("%3s%3s%3s  0" % (edge[0] + 1,
                                           edge[1] + 1,
                                           self.G[edge[0]][edge[1]]['type']))
        lines.append("M  END\n")

        with open(sdfname, "w") as f:
            f.write("\n".join(lines))


class Fragment:
    def __init__(self, mol, startatom):
        self.carried_atoms = nx.node_connected_component(mol.G, startatom - 1) # Numbering of atoms starts from 0
        self.G = mol.G.subgraph(self.carried_atoms)
        self.mol = mol

    def build_frame(self, central_atom, dir_atom, plane_atom):
        central_atom -= 1
        dir_atom -= 1
        plane_atom -= 1
        center_xyz = self.mol.G.nodes[central_atom]['xyz']
        dir_xyz = self.mol.G.nodes[dir_atom]['xyz']
        plane_xyz = self.mol.G.nodes[plane_atom]['xyz']
        
        xv = dir_xyz - center_xyz
        xv /= norm(xv)
        dir2 = plane_xyz - center_xyz
        zv = np.cross(xv, dir2)
        zv /= norm(zv)
        yv = np.cross(zv, xv)
        
        newframe = np.zeros((4, 4))
        newframe[:3, 0] = xv
        newframe[:3, 1] = yv
        newframe[:3, 2] = zv
        newframe[:3, 3] = center_xyz
        newframe[3, 3] = 1
        return newframe

    def translate_x(self, frame, length):
        translation_matrix = np.array([[1, 0, 0, length],
                                       [0, 1, 0, 0],
                                       [0, 0, 1, 0],
                                       [0, 0, 0, 1]])
        for atom in self.G.nodes:
            # print("%d) Before " % (atom+1) + repr(self.G.nodes[atom]['xyz']))
            self.G.nodes[atom]['xyz'] = (frame @ translation_matrix @ inv(frame) @ np.concatenate(
                (self.G.nodes[atom]['xyz'], [1]), axis=0))[:3]
            # print("%d) After " % (atom+1) + repr(self.G.nodes[atom]['xyz']))
        # return translation_matrix @ frame

    def rotate_x(self, frame, ang):
        rotmat = np.array([[1, 0, 0, 0],
                           [0, np.cos(ang), -np.sin(ang), 0],
                           [0, np.sin(ang), np.cos(ang), 0],
                           [0, 0, 0, 1]])
        for atom in self.G.nodes:
            # print("%d) Before " % (atom+1) + repr(self.G.nodes[atom]['xyz']))
            self.G.nodes[atom]['xyz'] = (frame @ rotmat @ inv(frame) @ np.concatenate((self.G.nodes[atom]['xyz'], [1]),
                                                                                      axis=0))[:3]
            # print("%d) After " % (atom+1) + repr(self.G.nodes[atom]['xyz']))
        # return rotmat @ frame

    def rotate_z(self, frame, ang):
        rotmat = np.array([[np.cos(ang), -np.sin(ang), 0, 0],
                           [np.sin(ang), np.cos(ang), 0, 0],
                           [0, 0, 1, 0],
                           [0, 0, 0, 1]])
        for atom in self.G.nodes:
            # print("%d) Before " % (atom+1) + repr(self.G.nodes[atom]['xyz']))
            self.G.nodes[atom]['xyz'] = (frame @ rotmat @ inv(frame) @ np.concatenate((self.G.nodes[atom]['xyz'], [1]),
                                                                                      axis=0))[:3]
            # print("%d) After " % (atom+1) + repr(self.G.nodes[atom]['xyz']))
        # return frame @ rotmat
