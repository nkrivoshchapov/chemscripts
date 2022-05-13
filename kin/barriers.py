from .names import Names
from .equations import equations_from_sheet
from ..utils import H2KC


def _sumup(datalist, moldata, key):
    res = []
    for idx, coeff in moldata:
        if coeff == 1:
            res.append(datalist[idx][key])
        else:
            res.append("%d*%s" % (coeff, datalist[idx][key]))
    return "+".join(res)


def obtain_barriers(sheet):
    eq_block = sheet.block(Names.EQ_BLOCK)
    mol_block = sheet.block(Names.MOL_BLOCK)
    equations = equations_from_sheet(sheet)

    for eq_idx, eq in enumerate(equations):
        assert eq_block['data'][eq_idx][Names.EQ_COL] == eq['original_line']
        reagent_idxs = []
        for mol in eq['reagents']:
            for cmol_idx, cmol in enumerate(mol_block['data']):
                if cmol[Names.MOLNAME_COL] == mol['name']:
                    reagent_idxs.append((cmol_idx, mol['n']))
                    break
        product_idxs = []
        for mol in eq['products']:
            for cmol_idx, cmol in enumerate(mol_block['data']):
                if cmol[Names.MOLNAME_COL] == mol['name']:
                    product_idxs.append((cmol_idx, mol['n']))
                    break
        # reagent_idxs and product_idxs contain pairs (Index of molecule in mol_block, stoichiometric coeffs)

        have_tsener = eq_block['data'][eq_idx][Names.TSENER_COL]
        have_fakebarrier = eq_block['data'][eq_idx][Names.FAKEB_COL]
        assert have_fakebarrier is not None or have_tsener is not None, "Not enough data to compute rate constants"
        assert not (have_fakebarrier is not None and have_tsener is not None), "Have TS energy AND fake barrier. Which one to use?!"

        if have_tsener is not None:
            ts_ener_cell = eq_block['cells'][eq_idx][Names.TSENER_COL]
            eq_block['data'][eq_idx][Names.FORW_COL] = "=%f*(%s-(%s))" % (H2KC, ts_ener_cell,
                                                       _sumup(mol_block['cells'], reagent_idxs, Names.E_COL))
            eq_block['data'][eq_idx][Names.BACKW_COL] = "=%f*(%s-(%s))" % (H2KC, ts_ener_cell,
                                                        _sumup(mol_block['cells'], product_idxs, Names.E_COL))
        elif have_fakebarrier is not None:
            barrier_cell = eq_block['cells'][eq_idx][Names.FAKEB_COL]
            dg = "(%s)-(%s)" % (_sumup(mol_block['cells'], product_idxs, Names.E_COL),
                                _sumup(mol_block['cells'], reagent_idxs, Names.E_COL))
            eq_block['data'][eq_idx][Names.FORW_COL] = "={fb}+{H2KC}*({dg}-ABS({dg}))/2".format(
                                                                            H2KC=H2KC, fb=barrier_cell, dg=dg)
            eq_block['data'][eq_idx][Names.BACKW_COL] = "={fb}-{H2KC}*({dg}+ABS({dg}))/2".format(
                                                                            H2KC=H2KC, fb=barrier_cell, dg=dg)
