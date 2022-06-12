from chemscripts import kin
import chemscripts.excelutils as xl
from chemscripts.energyutils import energies_to_excel
import sys, os

"""
Can be executed in 4 modes:
1) Parse equations from *.k and initialize the excel worksheet RAWTABLENAME
2) Calculate reaction barriers (saved to TABLENAME) from *.xlsx and create sctipt for kinetic integration SCRIPTNAME.  
3) Create sctipt for kinetic integration SCRIPTNAME from *.xlsx. Additional 'noprep' keyword is required.
4) Can also run GoodVibes calculations and save results to ENERGYFILE. 'energy' keyword is required 
"""

ENERGYFILE = "energies.xlsx"
RAWTABLENAME = "mytable.xlsx"
TABLENAME = "mytable_new.xlsx"
SCRIPTNAME = "calc_kinetics.py"

if __name__ == "__main__":
    firstarg = sys.argv[1]
    if firstarg == "energy":
        energies_to_excel("*.log", ENERGYFILE)
    elif firstarg.endswith(".k"):
        equations = kin.equations_from_k(firstarg)
        mols = kin.molecules_from_equations(equations)
        excelsheet = kin.initialize_sheet(mols, equations)
        excelsheet.save_xlsx(RAWTABLENAME)
    elif firstarg.endswith(".xlsx") and "noprep" in sys.argv:
        excelsheet = xl.ExcelSheet()
        excelsheet.read_xlsx(firstarg, get_values=True)
        with open(SCRIPTNAME, "w") as f:
            f.write(kin.create_script(excelsheet))
    elif firstarg.endswith(".xlsx"):        
        excelsheet = xl.ExcelSheet()
        excelsheet.read_xlsx(firstarg)
        kin.obtain_barriers(excelsheet)
        excelsheet.save_xlsx(TABLENAME, oldfile=firstarg)

        excelsheet = xl.ExcelSheet()
        excelsheet.read_xlsx(TABLENAME, get_values=True)
        with open(SCRIPTNAME, "w") as f:
            f.write(kin.create_script(excelsheet))
