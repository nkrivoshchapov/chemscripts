import sys, os, ntpath, string
import chemscripts.excelutils as xl
import chemscripts.energyutils as ener
from chemscripts import utils

LOG_PATHS = "./*/*.log"
XYZ_FILENAME = "geometries.xyz"

class Names(object):
    DIR_COL = "Directory"
    NAME_COL = "Logname"
    IDX_COL = "Index"
    TYPE_COL = "Type"

    # The values below are optional, they depend on the data you want the xyz-file to contain
    INTERM_BLOCK = "Intermediates"
    REL_ENERGY_COL = "dG relative, kcal/mol"
    TS_BLOCK = "Transition states"
    EA_ENERGY_COL = "dG of activation, a.u."
    DESCR_COL = "Description"
    CHR_COL = "Charge"
    MULT_COL = "Multiplicity"
    FORMATS = {
        INTERM_BLOCK: "{DESCR_COL}; Molecule {IDX_COL}; Chrg={CHR_COL}; Mult={MULT_COL}; Electronic energy={SCF_ENERGY_COL} a.u.; Free energy={QH_ENERGY_COL} a.u.; Relative energy={REL_ENERGY_COL} kcal/mol",
        TS_BLOCK: "{DESCR_COL}; Transition state {IDX_COL}; Chrg={CHR_COL}; Mult={MULT_COL}; Electronic energy={SCF_ENERGY_COL} a.u.; Free energy={QH_ENERGY_COL} a.u.; Activation free energy={EA_ENERGY_COL} kcal/mol",
    }
    
    # These columns will be added automatically (in the same order from left to right)
    ADD_COLUMNS = [DESCR_COL, CHR_COL, MULT_COL, TYPE_COL, IDX_COL]


if __name__ == "__main__":
    """
    1) No args -- obtain a table of energies in LOG_PATHS
    2) si_build.py smth.xlsx enum -- enumerate structures of different types
    3) si_build.py smth.xlsx build -- build xyz-file from given sheet
    """
    if len(sys.argv) == 1: # No args
        excelsheet = ener.get_energy_sheet(LOG_PATHS, get_scf=True)
        excelsheet.remove_key(ener._Names.ENERGY_BLOCK, ener._Names.C0_COL) # Delete concentration column

        # Add requested keys
        for key in reversed(Names.ADD_COLUMNS):
            excelsheet.add_key(ener._Names.ENERGY_BLOCK, key, position=0)

        # Split filenames into two columns
        excelsheet.add_key(ener._Names.ENERGY_BLOCK, Names.NAME_COL, position=0)
        excelsheet.add_key(ener._Names.ENERGY_BLOCK, Names.DIR_COL, position=0)
        block = excelsheet.block(ener._Names.ENERGY_BLOCK)
        for item in block['data']:
            item[Names.DIR_COL] = os.path.dirname(item[ener._Names.LOGNAME_COL])
            item[Names.NAME_COL] = ntpath.basename(item[ener._Names.LOGNAME_COL])
        excelsheet.remove_key(ener._Names.ENERGY_BLOCK, ener._Names.LOGNAME_COL)
        excelsheet.save_xlsx("energies.xlsx")
    elif sys.argv[1].endswith(".xlsx") and "enum" in sys.argv:
        tablename = sys.argv[1]
        excelsheet = xl.ExcelSheet()
        excelsheet.read_xlsx(tablename)

        counters = {}
        for block in excelsheet.datablocks:
            if Names.IDX_COL not in block['keys'] or Names.TYPE_COL not in block['keys']:
                print("Skipping the section '%s'" % block['name'])
                continue

            for item in block['data']:
                cur_type = item[Names.TYPE_COL]
                if cur_type not in counters:
                    counters[cur_type] = 1
                else:
                    counters[cur_type] += 1
                item[Names.IDX_COL] = counters[cur_type]

        print("Total number of structures of each type: " + repr(counters))
        excelsheet.save_xlsx(tablename.replace('.xlsx', '_enum.xlsx'), oldfile=tablename)
    elif sys.argv[1].endswith(".xlsx") and "build" in sys.argv:
        tablename = sys.argv[1]
        newname = xl.fix_xlsx(tablename)
        excelsheet = xl.ExcelSheet()
        excelsheet.read_xlsx(newname, get_values=True)
        print(repr(excelsheet.datablocks))
        xyzdata = [] # each element (string) corresponds to different structure
        for block in excelsheet.datablocks:
            assert block['name'] in Names.FORMATS.keys()
            for item in block['data']:
                keys = [t[1] for t in string.Formatter().parse(Names.FORMATS[block['name']]) if t[1] is not None]
                # print(repr(item))
                description_values = {}
                for key in keys:
                    sheet_key = None
                    if hasattr(ener._Names, key):
                        sheet_key = getattr(ener._Names, key) # SCF_ENERGY_COL and QH_ENERGY_COL are assigned here
                    elif hasattr(Names, key):
                        sheet_key = getattr(Names, key)
                    
                    if isinstance(item[sheet_key], float):
                        description_values[key] = "%12.8f" % item[sheet_key]
                    else:
                        description_values[key] = item[sheet_key]
                description = Names.FORMATS[block['name']].format(**description_values)
                
                print("Reading geometry from " + os.path.join(item[Names.DIR_COL], item[Names.NAME_COL]))
                xyzs, syms = utils.parse_log(os.path.join(item[Names.DIR_COL], item[Names.NAME_COL]))
                xyzdata.append(utils.to_xyz(xyzs, syms, description=description))

        with open(XYZ_FILENAME, 'w') as f:
            f.write('\n'.join(xyzdata))
