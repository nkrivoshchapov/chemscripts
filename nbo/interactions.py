import sys, os, glob, ntpath, re
import numpy as np


START_RE = re.compile(r"\(occupancy\).*bond.*orbital.*\/.*coefficients.*\/.*hybrids", re.IGNORECASE)
END_RE = re.compile("nho.*directionality.*and.*bond.*bending", re.IGNORECASE)
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


def get_scf_energy(lines):
    for line in reversed(lines):
        if "SCF Done" in line:
            return float(line.split('=')[1].split('A.U.')[0])


def get_deletion_energy(logfile):
    lines = open(logfile, "r").readlines()
    for line in reversed(lines):
        if "Energy of deletion :" in line:
            return float(line.split(':')[1].replace('\n', ''))


class NboCalculation:
    def __init__(self, logname, keys, donor_patterns, acceptor_patterns):
        self.donor_patterns = donor_patterns
        self.acceptor_patterns = acceptor_patterns
        self.keys = keys
        self.logname = logname
        self.molname = ntpath.basename(logname).split('.')[0].split('_')[0]
        self.nbo47name = logname.replace('.log', '.47')
        self.fnboname = logname.replace('.log', '.fnbo')
        self.dmnboname = logname.replace('.log', '.dmnbo')
        self.aonboname = logname.replace('.log', '.aonbo')
        self.nbooutname = logname.replace('.log', '.out')
        self.delname = logname.replace('_nbo.log', '_del.log')

    def parse_nbo_info(self, loglines):
        # Get NBasis
        self.nbasis = None
        for line in loglines:
            if 'NBasis' in line and 'RedAO' in line:
                self.nbasis = int(line.split('=')[1].split('RedAO')[0])
                break
        assert self.nbasis is not None

        # Get NBO symbols
        start_idx = None
        end_idx = None
        nbolines = open(self.nbooutname, 'r').readlines()
        for i, line in enumerate(nbolines):
            if re.search(START_RE, line):
                start_idx = i
            elif re.search(END_RE, line):
                end_idx = i
                break
        assert start_idx is not None and end_idx is not None

        self.donors = []
        self.acceptors = []
        for i, line in enumerate(nbolines[start_idx:end_idx], start=start_idx):
            for donor_pattern in self.donor_patterns[self.molname]:
                if re.search(donor_pattern, line):
                    self.donors.append(int(line.split('.')[0]))
            for acc_pattern in self.acceptor_patterns[self.molname]:
                if re.search(acc_pattern, line):
                    self.acceptors.append(int(line.split('.')[0]))

    def get_df_item(self):
        res = {}
        loglines = open(self.logname, "r").readlines()
        self.parse_nbo_info(loglines)

        fao = NboSymmMatrix(self.nbo47name, self.nbasis, section='$FOCK')
        fnbo = NboSymmMatrix(self.fnboname, self.nbasis)
        dmao = NboSymmMatrix(self.nbo47name, self.nbasis, section='$DENSITY')
        dmnbo = NboSymmMatrix(self.dmnboname, self.nbasis)
        aonbo = NboNonSymmMatrix(self.aonboname, self.nbasis)
        chao = NboSymmMatrix(self.logname, self.nbasis)

        if 'ScfEner' in self.keys:
            res['ScfEner'] = get_scf_energy(loglines)

        if 'Del_total' in self.keys:
            res['Del_total'] = get_deletion_energy(self.delname)

        for totalkey in TOTAL_KEYS:
            if totalkey in self.keys:
                res[totalkey] = 0 # Initialize before summing up

        for i in range(self.nbasis):
            for j in range(i):
                for totalkey in TOTAL_KEYS:
                    if totalkey in self.keys:
                        res[totalkey] += INT_FORMULAS[totalkey.split('_')[0]](fnbo, dmnbo, i, j)

        for sumkey in SUM_KEYS:
            if sumkey in self.keys:
                res[sumkey] = 0 # Initialize before summing up

        for donor in self.donors:
            for acceptor in self.acceptors:
                for sumkey in SUM_KEYS:
                    if sumkey in self.keys:
                        res[sumkey] += INT_FORMULAS[sumkey.split('_')[0]](fnbo, dmnbo, donor - 1, acceptor - 1)
        return res
