from numpy import exp
from .equations import _parse_equation
from .names import Names


h = 6.62607E-34
k = 1.38065E-23
R = 0.001987204
def rateconst(dG, T):
    return k*T/h*exp(-dG/T/R)


def _get_c0(datalist, molname):
    res = None
    for item in datalist:
        if item[Names.MOLNAME_COL] == molname:
            res = item[Names.CO_COL]
            break
    assert res is not None
    return res


def create_script(sheet, resfile="res.csv"):
    eq_block = sheet.block(Names.EQ_BLOCK)
    mol_block = sheet.block(Names.MOL_BLOCK)
    setup_block = sheet.block(Names.SETUP_BLOCK)
    assert len(setup_block['data']) == 1, "Expected exactly one row in '%s' table" % Names.SETUP_BLOCK
    setup_data = setup_block['data'][0]
    temperature = setup_data[Names.TEMP_COL]

    eq_list = []
    for item in eq_block['data']:
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
            v_str += "\nk%d = %.15E" % (j, rateconst(eq_block['data'][eq_idx][Names.FORW_COL], temperature))
            v_str += "\nk%dr = %.15E" % (j, rateconst(eq_block['data'][eq_idx][Names.BACKW_COL], temperature))
            reac_list.append(v_str)

    ofstr = "from numpy import *\nimport scipy.integrate as itg\nimport os\n\n"
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
        c0 = _get_c0(mol_block['data'], chem_dict_r[n])
        if n != i:
            ofstr += "#" + chem_dict_r[n] + "\n"
            ofstr += "%f,\\\n" % c0
        else:
            ofstr += "#" + chem_dict_r[n] + "\n"
            ofstr += "%f,\\\n])\n\n" % c0

    for term in reac_list:
        ofstr += term + "\n\n"

    ofstr += """TIMESTEP = {timestep} # sec
MAXTIME = {maxtime} # sec
RESFILE = "{resfile}"
NINT = 1
if os.path.exists(RESFILE):
    os.remove(RESFILE)
csvfile=open(RESFILE, 'a')
csvfile.write(";".join(["Time, s", "{mols}"])+"\\n")

mainT = 0
print_count = 0
write_count = 0
while mainT < MAXTIME:
    print_count += 1
    write_count += 1
    mainT += TIMESTEP
    t = arange(0, (1 + NINT) * TIMESTEP, TIMESTEP)
    Y = itg.odeint(dy, y0, t)
    y0 = Y[Y.shape[0] - 1]
    y0_str = [str(mainT)] + [str(i) for i in y0]

    if write_count == 1:
        csvfile.write(";".join(y0_str)+"\\n")
        write_count = 0
    if print_count == 50:
        print("%f seconds have passed (%f %%)\\nConcentrations = %s" % (mainT, mainT/MAXTIME*100, ",".join(y0_str)))
        print_count = 0
csvfile.close()
print("Kinetic modelling of {maxtime} seconds is finished")
""".format(mols='", "'.join([chem_dict_r[n] for n in range(0, i + 1)]),
                          timestep=setup_data[Names.STEP_COL],
                          maxtime=setup_data[Names.MAXTIME_COL],
                          resfile=resfile
                          )
    return ofstr
