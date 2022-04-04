import sys
from numpy import exp
from math import log
from openpyxl import load_workbook, Workbook
from .equations import _parse_equation
from .names import Names


h = 6.62607E-34
k = 1.38065E-23
T = 273.15
R = 0.001987204
def rateconst(dG):
    return k*T/h*exp(-dG/T/R)


def _get_c0(datalist, molname):
    res = None
    for item in datalist:
        if item[Names.MOLNAME_COL] == molname:
            res = item[Names.CO_COL]
    assert res is not None
    return res

def create_script(sheet):
    eq_block_idx = None
    for i, db in enumerate(sheet.datablocks):
        if db['name'] == Names.EQ_BLOCK:
            eq_block_idx = i
    assert eq_block_idx is not None

    mol_block_idx = None
    for i, db in enumerate(sheet.datablocks):
        if db['name'] == Names.MOL_BLOCK:
            mol_block_idx = i
    assert mol_block_idx is not None

    eq_list = []
    for item in sheet.datablocks[eq_block_idx]['data']:
        fixed_line, _, _ = _parse_equation(item[Names.EQ_COL])
        eq_list.append(fixed_line.replace("+", " + ").replace("<->", " <-> "))

    chem_dict = {}
    reac_dict = {}
    i = -1
    j = -1

    reac_list = []

    for eq_idx, line in enumerate(eq_list):
        if line == '':
            continue
        elif line[0] == '#':
            continue
        else:
            reac = []
            mol_reac = []
            mol_prod = []
            prod = []
            rate_str = ""
            v_str = ""
            R = True
            N = False
            for term in line.split():
                if term == '+':
                    continue
                elif term == "<->":
                    R = False
                    continue
          
                mol_list = term.split("*")
                if len(mol_list) == 1:
                    mol_list = [1] + mol_list
                else:
                    mol_list = [int(mol_list[0])] + mol_list[1:]
                mol = mol_list[1]
                stoi = int(mol_list[0])
                
                if not mol in chem_dict:
                    i += 1
                    chem_dict.update({mol : i})
                    reac_dict.update({chem_dict[mol] : ""})

                if R:
                    reac.append(mol_list)
                    mol_reac.append(mol)
                else:
                    prod.append(mol_list)
                    mol_prod.append(mol)

            j += 1
            reac_args = list(set(mol_reac + mol_prod)) #Remove repeating elements
            rate_str = "v_" + str(j) + "("
            for term in reac_args:
                rate_str += "y[" + str(chem_dict[term]) + "], "        
            rate_str = rate_str[0:-2] + ")"

            for term in reac:
                reac_dict[chem_dict[term[1]]] += " -" + str(term[0]) + "*" + rate_str
                
            for term in prod:
                reac_dict[chem_dict[term[1]]] += " +" + str(term[0]) + "*" + rate_str

            v_str = "#" + line + "\n"
            v_str += "v_" + str(j) + " = lambda "
            for term in reac_args:
                v_str += term + ", "
            v_str = v_str[0:-2] + " : "
            
            v_str += "k" + str(j) + " * "
        
            for term in reac:
                v_str += term[1] + "**" + str(term[0]) +  " * "
            v_str = v_str[0:-3]
                
            v_str += " - k" + str(j) + "r * "
            for term in prod:
                v_str += term[1] + "**" + str(term[0]) + " * "
            v_str = v_str[0:-3]
            v_str += "\nk%d = %.15E" % (j, rateconst(sheet.datablocks[eq_block_idx]['data'][eq_idx][Names.FORW_COL]))
            v_str += "\nk%dr = %.15E" % (j, rateconst(sheet.datablocks[eq_block_idx]['data'][eq_idx][Names.BACKW_COL]))
            reac_list.append(v_str)

    ofstr = "from numpy import *\nimport scipy.integrate as itg\n\n"
    chem_dict_r = {}
    for term in chem_dict:
        chem_dict_r.update({chem_dict[term] : term})
    # ofstr += "'''\n## Reaction ##\n\n"
    # ofstr += st + "\n\n"
    # ofstr += "## Mapping ##\n\n"
    # for term in reac_dict:
        # ofstr += chem_dict_r[term] + "\t" + str(term) + "\t" + reac_dict[term] + "\n"
    # ofstr += "'''\n\n"
        
    ost = "dy = lambda y, t: array([\\\n"

    for term in reac_dict:
        ost += reac_dict[term] + ",\\\n"
    ost = ost[0:-3] + "\\\n])"
    ofstr += ost
    ofstr += "\n\n#Initial concentrations:\ny0 = array([\\\n"
    for n in range(0, i + 1):
        c0 = _get_c0(sheet.datablocks[mol_block_idx]['data'], chem_dict_r[n])
        if n != i:
            ofstr += "#" + chem_dict_r[n] + "\n"
            ofstr += "%f,\\\n" % c0
        else:
            ofstr += "#" + chem_dict_r[n] + "\n"
            ofstr += "%f,\\\n])\n\n" % c0

    for term in reac_list:
        ofstr += term + "\n\n"

    ofstr += """mainT = 0
csvfile=open("res.csv","w")
csvfile.write(";".join(["{mols}"])+"\\n")
csvfile.close()
csvfile=open("res.csv","a")
count = 0
count2 = 0
while mainT < 3600000:
    count += 1
    count2 += 1
    t = arange(0, 0.2, 0.1)
    Y = itg.odeint(dy, y0, t)
    y0 = Y[len(Y)-1]
    
    y00= []
    for i in y0:
        y00.append(str(i))
    if count2 == 20:
        csvfile.write(";".join(y00)+"\\n")
        count2 = 0
    if count == 1000:
        print(str(mainT/360000)+"%    " + repr(y0))
        count = 0
    mainT += 1
csvfile.close()""".format(mols='","'.join([chem_dict_r[n] for n in range(0, i + 1)]))
    return ofstr
