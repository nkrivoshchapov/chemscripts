import sys
import myscripts.excelutils as msexcel
from myscripts.energyutils import get_energy_sheet, prepare_sheet

"""
Runs GoodVibes free energy calculations from Gaussian logs
Saves results as Excel table
"""

if __name__ == "__main__":
    if len(sys.argv) == 1: # No args
        get_energy_sheet("./*.log").save_xlsx("energies.xlsx")
    elif sys.argv[1].endswith(".xlsx"):
        excelsheet = msexcel.ExcelSheet()
        excelsheet.read_xlsx(sys.argv[1], get_values=True)
        get_energy_sheet(excelsheet=excelsheet).save_xlsx("energies.xlsx")
    elif sys.argv[1] == "prep":
        prepare_sheet("./*.log").save_xlsx("concentrations.xlsx")
