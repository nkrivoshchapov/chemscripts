import sys, os, glob, ntpath, re
import numpy as np

from .logparsers import NBO3LogParser, NBO6LogParser

E2_EDIFF_THRESHOLD = 0.00001


def e2_energy(fnbo, dmnbo, i, j):
    occup_sum = dmnbo.arr[i][i] + dmnbo.arr[j][j]
    if occup_sum <= 2:
        if abs(fnbo.arr[i][i] - fnbo.arr[j][j]) < E2_EDIFF_THRESHOLD:
            res = 0
        else:
            res = occup_sum * fnbo.arr[i][j] ** 2 / abs(fnbo.arr[i][i] - fnbo.arr[j][j])
    else:
        if abs(fnbo.arr[i][i] - fnbo.arr[j][j]) < E2_EDIFF_THRESHOLD:
            res = 0
        else:
            res = (4 - occup_sum) * fnbo.arr[i][j] ** 2 / abs(fnbo.arr[i][i] - fnbo.arr[j][j])
    return res


def df_energy(fnbo, dmnbo, i, j):
    return -dmnbo.arr[i][j] * fnbo.arr[i][j]


INT_FORMULAS = {
                  'E2': e2_energy,
                  'DF': df_energy,
               }
TOTAL_KEYS = [(key + "_total") for key in INT_FORMULAS.keys()]
SUM_KEYS = [(key + "_sum") for key in INT_FORMULAS.keys()]


class NboNonSymmMatrix:
    def __init__(self, nicefile, nbasis, section=None):
        self.n = nbasis
        self.arr = np.empty([self.n, self.n])
        self.count = 0
        self.cur_row = 0
        self.cur_col = 0
        nicelines = open(nicefile, "r").readlines()
        for line in nicelines[3:]:
            parts = line.split()
            for part in parts:
                self.append_elem(float(part))
            if self.cur_row == self.n and self.cur_col == self.n - 1:
                break

    def append_elem(self, newnum):
        if self.cur_row == self.n:
            self.cur_row = 0
            self.cur_col += 1
        self.arr[self.cur_row][self.cur_col] = newnum
        self.cur_row += 1


class NboSymmMatrix:
    def __init__(self, nicefile, nbasis, section=None):
        self.n = nbasis
        self.arr = np.empty([self.n, self.n])
        self.count = 0
        self.cur_row = 0
        self.cur_col = 0

        if nicefile.endswith('.log'):
            # No longer nice - It's Gaussian logfile
            lines = open(nicefile, "r").readlines()
            coreH_start = None
            for i, line in enumerate(lines):
                if "****** Core Hamiltonian ******" in line:
                    coreH_start = i
                    break
            assert coreH_start is not None
            coreH_end = None
            for i, line in enumerate(lines):
                if i < coreH_start:
                    continue
                if "SVDSVc" in line or "GSVD" in line or "Symmetry" in line:
                    coreH_end = i
                    break
            assert coreH_end is not None
            columns = None
            for line in lines[coreH_start+1:coreH_end]:
                parts = line.split()
                if "D" in line:
                    assert columns is not None
                    row_idx = int(parts[0])
                    for i, part in enumerate(parts[1:]):
                        col_idx = columns[i]
                        value = float(part.replace('D', 'E'))
                        self.arr[row_idx - 1][col_idx - 1] = value
                        self.arr[col_idx - 1][row_idx - 1] = value
                else:
                    columns = [int(p) for p in parts]
        elif nicefile.endswith('.47'):
            lines = open(nicefile, "r").readlines()
            mysection_start = None
            mysection_end = None
            for i, line in enumerate(lines):
                if section in line:
                    mysection_start = i
                elif "$END" in line and mysection_start is not None:
                    mysection_end = i
                    break
            assert mysection_start is not None and mysection_end is not None
            for line in lines[mysection_start + 1:mysection_end]:
                parts = line.split()
                for part in parts:
                    self.append_elem(float(part))
        else:
            nicelines = open(nicefile, "r").readlines()
            for line in nicelines[3:]:
                parts = line.split()
                for part in parts:
                    self.append_elem(float(part))

    def append_elem(self, newnum):
        if self.cur_col > self.cur_row:
            self.cur_row += 1
            self.cur_col = 0
        self.arr[self.cur_row][self.cur_col] = newnum
        self.arr[self.cur_col][self.cur_row] = newnum
        self.cur_col += 1


class NboCalculation:
    def __init__(self, nboname):
        self.nboname = nboname # Either log of g16/nbo3 calculation or nbo6 out-file

        if self.nboname.endswith('.log'):
            self.parser = NBO3LogParser(self.nboname)
            self.suffix = '.log'
        elif self.nboname.endswith('.out'):
            self.parser = NBO6LogParser(self.nboname)
            self.suffix = '.out'
        # print(repr(self.parser.NBOs))
        # raise Exception("kik")
        self.nbasis = self.parser.nbasis

        if self.suffix == '.log' and os.path.isfile(nboname.replace('.log', '.out')):
            self.outname = nboname.replace('.log', '.out')

        if self.suffix == '.out' and os.path.isfile(nboname.replace('.out', '.log')):
            self.logname = nboname.replace('.out', '.log')

        if self.nboname.endswith('.log') or hasattr(self, "logname"):
            self.scfener = self.obtain_scf_energy()
            # self.chao = NboSymmMatrix(self.nboname, self.nbasis) if self.nboname.endswith('.log') \
            #             else NboSymmMatrix(self.logname, self.nbasis)

        nbo47name = self.nboname.replace(self.suffix, '.47')
        if os.path.isfile(nbo47name):
            self.nbo47name = nbo47name
            self.fao = NboSymmMatrix(self.nbo47name, self.nbasis, section='$FOCK')
            self.dmao = NboSymmMatrix(self.nbo47name, self.nbasis, section='$DENSITY')

        fmoname = self.nboname.replace(self.suffix, '.fmo')
        if os.path.isfile(fmoname):
            self.fmoname = fmoname
            self.fmo = NboSymmMatrix(self.fmoname, self.nbasis)

        dmmoname = self.nboname.replace(self.suffix, '.dmmo')
        if os.path.isfile(dmmoname):
            self.dmmoname = dmmoname
            self.dmmo = NboSymmMatrix(self.dmmoname, self.nbasis)

        fnboname = self.nboname.replace(self.suffix, '.fnbo')
        if os.path.isfile(fnboname):
            self.fnboname = fnboname
            self.fnbo = NboSymmMatrix(self.fnboname, self.nbasis)

        dmnboname = self.nboname.replace(self.suffix, '.dmnbo')
        if os.path.isfile(dmnboname):
            self.dmnboname = dmnboname
            self.dmnbo = NboSymmMatrix(self.dmnboname, self.nbasis)

        aonboname = self.nboname.replace(self.suffix, '.aonbo')
        if os.path.isfile(aonboname):
            self.aonboname = aonboname
            self.aonbo = NboNonSymmMatrix(self.aonboname, self.nbasis)

        delname = self.nboname.replace('_nbo' + self.suffix, '_del.log')
        if os.path.isfile(delname):
            self.delname = delname
            self.delener = self.obtain_deletion_energy()

    def obtain_scf_energy(self):
        if self.nboname.endswith('.log'):
            lines = open(self.nboname, 'r').readlines()
        elif hasattr(self, "logname"):
            lines = open(self.logname, 'r').readlines()
        else:
            raise Exception(RuntimeError)
        for line in reversed(lines):
            if "SCF Done" in line:
                return float(line.split('=')[1].split('A.U.')[0])

    def obtain_deletion_energy(self):
        lines = open(self.delname, "r").readlines()
        for line in reversed(lines):
            if "Energy change :" in line:
                return float(line.split(':')[1].split('a.u.')[0])

    def get_data(self, keys=("E2_sum",), donor_patterns=(), acceptor_patterns=(), donors=(), acceptors=()):
        # print("My donors pat = " + repr(donor_patterns))
        # print("My acceptors pat = " + repr(acceptor_patterns))
        if len(donor_patterns) > 0 and len(acceptor_patterns) > 0:
            assert len(donors) == 0 and len(acceptors) == 0, "NBO indices and patterns were given simultaniously!"
            donors = []
            for donor_pat in donor_patterns:
                donors += [item['index'] for item in self.parser.find_by_regex(donor_pat)]
            acceptors = []
            for acceptor_pat in acceptor_patterns:
                acceptors += [item['index'] for item in self.parser.find_by_regex(acceptor_pat)]

        # print("My donors = " + repr(donors))
        # print("My acceptors = " + repr(acceptors))
        # raise Exception("kek")
        res = {}
        if 'ScfEner' in keys:
            res['ScfEner'] = self.scfener
        if 'Del_total' in keys:
            res['Del_total'] = self.delener

        for totalkey in TOTAL_KEYS:
            if totalkey in keys:
                res[totalkey] = 0 # Initialize before summing up

        for i in range(self.nbasis):
            for j in range(i):
                for totalkey in TOTAL_KEYS:
                    if totalkey in keys:
                        res[totalkey] += INT_FORMULAS[totalkey.split('_')[0]](self.fnbo, self.dmnbo, i, j)

        for sumkey in SUM_KEYS:
            if sumkey in keys:
                res[sumkey] = 0 # Initialize before summing up

        for donor in donors:
            for acceptor in acceptors:
                for sumkey in SUM_KEYS:
                    if sumkey in keys:
                        res[sumkey] += INT_FORMULAS[sumkey.split('_')[0]](self.fnbo, self.dmnbo,
                                                                          donor - 1, acceptor - 1)
        return res

