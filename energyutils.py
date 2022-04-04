from .excelutils import ExcelSheet
from .utils import get_goodvibes_g
import glob

class _Names(object):
    ENERGY_BLOCK = "Energies"
    LOGNAME_COL = "Filename"
    QH_ENERGY_COL = "QH free energy, a.u."

def get_energy_sheet(filemask):
    excelsheet = ExcelSheet()
    excelsheet.add_block(blockname=_Names.ENERGY_BLOCK,
                         cols=[_Names.LOGNAME_COL,
                               _Names.QH_ENERGY_COL])
    for file in glob.glob(filemask):
        excelsheet.add_row(blockname=_Names.ENERGY_BLOCK, data={
                                _Names.LOGNAME_COL: file,
                                _Names.QH_ENERGY_COL: get_goodvibes_g(file)})
    return excelsheet

def energies_to_excel(filemask, excelname):
    mysheet = get_energy_sheet(filemask)
    mysheet.save_xlsx(excelname)
