from .excelutils import ExcelSheet
from .utils import get_goodvibes_g, get_gaussian_scfener
import glob


class _Names(object):
    ENERGY_BLOCK = "Energies"
    LOGNAME_COL = "Filename"
    C0_COL = "C0, M"
    SCF_ENERGY_COL = "E, a.u."
    QH_ENERGY_COL = "QH free energy, a.u."


def get_columns(get_scf: bool, get_free: bool):
    res = [_Names.LOGNAME_COL]
    if get_scf:
        res += [_Names.SCF_ENERGY_COL]
    if get_free:
        res += [_Names.C0_COL,
                _Names.QH_ENERGY_COL]
    return res

def get_energy_sheet(filemask=None, excelsheet=None, get_scf=False, get_free=True):
    assert filemask is not None or excelsheet is not None
    assert not (filemask is not None and excelsheet is not None)
    if filemask is not None:


        excelsheet = ExcelSheet()
        excelsheet.add_block(blockname=_Names.ENERGY_BLOCK,
                             cols=get_columns(get_scf, get_free))

        for file in glob.glob(filemask):
            newitem = {_Names.LOGNAME_COL: file}
            if get_free:
                newitem[_Names.C0_COL] = 1.0
            excelsheet.add_row(blockname=_Names.ENERGY_BLOCK, data=newitem)

    energy_block = excelsheet.block(_Names.ENERGY_BLOCK)
    for item in energy_block['data']:
        if _Names.QH_ENERGY_COL in energy_block['keys']:
            item[_Names.QH_ENERGY_COL] = get_goodvibes_g(item[_Names.LOGNAME_COL], conc=item[_Names.C0_COL])
        if _Names.SCF_ENERGY_COL in energy_block['keys']:
            item[_Names.SCF_ENERGY_COL] = get_gaussian_scfener(item[_Names.LOGNAME_COL])
    return excelsheet


def energies_to_excel(filemask, excelname): # This is useless
    mysheet = get_energy_sheet(filemask)
    mysheet.save_xlsx(excelname)


def prepare_sheet(filemask, get_scf=False, get_free=True):
    excelsheet = ExcelSheet()
    excelsheet.add_block(blockname=_Names.ENERGY_BLOCK,
                         cols=get_columns(get_scf))
    for file in glob.glob(filemask):
        newitem = {_Names.LOGNAME_COL: file}
        if get_free:
            newitem[_Names.C0_COL] = 1.0
        excelsheet.add_row(blockname=_Names.ENERGY_BLOCK, data=newitem)
    return excelsheet
