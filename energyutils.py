from .excelutils import ExcelSheet
from .utils import get_goodvibes_g
import glob


class _Names(object):
    ENERGY_BLOCK = "Energies"
    LOGNAME_COL = "Filename"
    C0_COL = "C0, M"
    QH_ENERGY_COL = "QH free energy, a.u."


def get_energy_sheet(filemask=None, excelsheet=None):
    assert filemask is not None or excelsheet is not None
    assert not (filemask is not None and excelsheet is not None)
    if filemask is not None:
        excelsheet = ExcelSheet()
        excelsheet.add_block(blockname=_Names.ENERGY_BLOCK,
                             cols=[_Names.LOGNAME_COL,
                                   _Names.C0_COL,
                                   _Names.QH_ENERGY_COL])

        for file in glob.glob(filemask):
            excelsheet.add_row(blockname=_Names.ENERGY_BLOCK, data={
                                    _Names.LOGNAME_COL: file,
                                    _Names.C0_COL: 1.0})

    energy_block = excelsheet.block(_Names.ENERGY_BLOCK)
    for item in energy_block['data']:
        item[_Names.QH_ENERGY_COL] = get_goodvibes_g(item[_Names.LOGNAME_COL], conc=item[_Names.C0_COL])
    return excelsheet


def energies_to_excel(filemask, excelname): # This is useless
    mysheet = get_energy_sheet(filemask)
    mysheet.save_xlsx(excelname)


def prepare_sheet(filemask):
    excelsheet = ExcelSheet()
    excelsheet.add_block(blockname=_Names.ENERGY_BLOCK,
                         cols=[_Names.LOGNAME_COL,
                               _Names.C0_COL,
                               _Names.QH_ENERGY_COL])
    for file in glob.glob(filemask):
        excelsheet.add_row(blockname=_Names.ENERGY_BLOCK, data={
                                _Names.LOGNAME_COL: file,
                                _Names.C0_COL: 1.0})
    return excelsheet
