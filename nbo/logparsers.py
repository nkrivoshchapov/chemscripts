import re

"""
Some ideas for these classes:
    1) Put get_choose.py here; for Gaussian (1) and for NBO (2)
    2) Conversion from nbo data to molecular graph (interface with geom/molfrag)
    3) Hybridization data
"""
BOND_PATTERN = re.compile(r"BD\(.*\).*-.*", re.IGNORECASE)


class NBOLogParser:
    def __init__(self, logfile):
        self.NBOs = [] # Elements are dicts
                       # Keys: 'index' (starts with 1), 'name'
        self.parse_log(logfile)
        self.nbasis = len(self.NBOs)

    def parse_log(self, file):
        raise Exception(NotImplementedError)

    def find_by_regex(self, pattern):
        if isinstance(pattern, str):
            pattern = re.compile(pattern)
        assert isinstance(pattern, re.Pattern)

        result = []
        for nbo in self.NBOs:
            if re.search(pattern, nbo['name']):
                result.append(nbo)
        return result

    def get_bondlist(self):
        nbos = self.find_by_regex(BOND_PATTERN)
        bonds = []
        for nbo in nbos:
            atoms = nbo['name'].split(')')[1].split('-')
            idxs = [re.sub('[^0-9]','', a) for a in atoms]
            bonds.append(idxs)
        return bonds


class NBO6LogParser(NBOLogParser):
    def __init__(self, logfile):
        self.NBO_BEGIN = re.compile(r"\(occupancy\).*bond.*orbital.*\/.*coefficients.*\/.*hybrids", re.IGNORECASE)
        self.NBO_END = re.compile("nho.*directionality.*and.*bond.*bending", re.IGNORECASE)
        super().__init__(logfile)

    def parse_log(self, file):
        lines = open(file, 'r').readlines()
        start_idx = None
        end_idx = None
        for i, line in enumerate(lines):
            if re.search(self.NBO_BEGIN, line):
                start_idx = i
            elif re.search(self.NBO_END, line):
                end_idx = i
                break
        assert start_idx is not None and end_idx is not None

        for line in lines[start_idx + 2:end_idx - 2]:
            indexstr = line[:4].replace(' ', '')
            if len(indexstr) > 0 and indexstr != '---':
                nbo_index = int(indexstr)
                nbo_name = line[16:34].replace(' ', '')
                newitem = {
                    'index': nbo_index,
                    'name': nbo_name,
                }
                if "BD" not in nbo_name:
                    p_character = line[57:63]
                    if len(p_character) == 0:
                        p_character = 0.0
                    else:
                        p_character = float(p_character)
                    newitem['p_character'] = p_character / 100 # Ratios (not percents)
                self.NBOs.append(newitem)


class NBO3LogParser(NBOLogParser):
    def __init__(self, logfile):
        self.NBO_BEGIN = re.compile(r"\(occupancy\).*bond.*orbital.*\/.*coefficients.*\/.*hybrids", re.IGNORECASE)
        self.NBO_END = re.compile("nho.*directionality.*and.*bond.*bending", re.IGNORECASE)
        super().__init__(logfile)

    def parse_log(self, file):
        lines = open(file, 'r').readlines()
        start_idx = None
        end_idx = None
        for i, line in enumerate(lines):
            if re.search(self.NBO_BEGIN, line):
                start_idx = i
            elif re.search(self.NBO_END, line):
                end_idx = i
                break
        assert start_idx is not None and end_idx is not None

        for line in lines[start_idx + 2:end_idx - 2]:
            indexstr = line[:6].replace(' ', '')
            if len(indexstr) > 0:
                nbo_index = int(indexstr)
                nbo_name = line[18:39].replace(' ', '')
                newitem = {
                              'index': nbo_index,
                              'name': nbo_name,
                          }
                if "BD" not in nbo_name:
                    p_character = line[59:65]
                    if len(p_character) == 0:
                        p_character = 0.0
                    else:
                        p_character = float(p_character)
                    newitem['p_character'] = p_character / 100 # Ratios (not percents)
                self.NBOs.append(newitem)


def highest_p_character(nbolist):
    return sorted(nbolist, key=lambda x: x['p_character'], reverse=True)[0]
