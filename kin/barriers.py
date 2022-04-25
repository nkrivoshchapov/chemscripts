from .names import Names
from .equations import equations_from_sheet
from ..utils import H2KC


def _sumup(datalist, indices, key):
    res = []
    for idx in indices:
        res.append(datalist[idx][key])
    return "+".join(res)


def obtain_barriers(sheet):
    equations = equations_from_sheet(sheet)
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

    for eq_idx, eq in enumerate(equations):
        assert sheet.datablocks[eq_block_idx]['data'][eq_idx][Names.EQ_COL] == eq['original_line']
        reagent_idxs = []
        for mol in eq['reagents']:
            for cmol_idx, cmol in enumerate(sheet.datablocks[mol_block_idx]['data']):
                if cmol[Names.MOLNAME_COL] == mol['name']:
                    reagent_idxs.append(cmol_idx)
                    break
        product_idxs = []
        for mol in eq['products']:
            for cmol_idx, cmol in enumerate(sheet.datablocks[mol_block_idx]['data']):
                if cmol[Names.MOLNAME_COL] == mol['name']:
                    product_idxs.append(cmol_idx)
                    break
        have_tsener = sheet.datablocks[eq_block_idx]['data'][eq_idx][Names.TSENER_COL]
        have_fakebarrier = sheet.datablocks[eq_block_idx]['data'][eq_idx][Names.FAKEB_COL]
        assert have_fakebarrier or have_tsener, "Not enough data to compute rate constants"
        assert not (have_fakebarrier and have_tsener), "Have TS energy AND fake barrier. Which one to use?!"
        if have_tsener is not None:
            ts_ener_cell = sheet.datablocks[eq_block_idx]['cells'][eq_idx][Names.TSENER_COL]
            sheet.datablocks[eq_block_idx]['data'][eq_idx][Names.FORW_COL] = "=%f*(%s-(%s))" % (H2KC, ts_ener_cell,
                                          _sumup(sheet.datablocks[mol_block_idx]['cells'], reagent_idxs, Names.E_COL))
            sheet.datablocks[eq_block_idx]['data'][eq_idx][Names.BACKW_COL] = "=%f*(%s-(%s))" % (H2KC, ts_ener_cell,
                                          _sumup(sheet.datablocks[mol_block_idx]['cells'], product_idxs, Names.E_COL))
        elif have_fakebarrier is not None:
            barrier_cell = sheet.datablocks[eq_block_idx]['cells'][eq_idx][Names.FAKEB_COL]
            dg = "(%s)-(%s)" % (_sumup(sheet.datablocks[mol_block_idx]['cells'], product_idxs, Names.E_COL),
                                _sumup(sheet.datablocks[mol_block_idx]['cells'], reagent_idxs, Names.E_COL))
            sheet.datablocks[eq_block_idx]['data'][eq_idx][Names.FORW_COL] = "={fb}+{H2KC}*({dg}-ABS({dg}))/2".format(
                                                                            H2KC=H2KC, fb=barrier_cell, dg=dg)
            sheet.datablocks[eq_block_idx]['data'][eq_idx][Names.BACKW_COL] = "={fb}-{H2KC}*({dg}+ABS({dg}))/2".format(
                                                                            H2KC=H2KC, fb=barrier_cell, dg=dg)
