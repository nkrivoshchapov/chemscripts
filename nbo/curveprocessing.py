from math import sqrt
from scipy.optimize import minimize
import pandas as pd
import scipy.stats

def get_coeff(df, main, ref):
    def myrmsd(args):
        a = args[0]
        return sqrt(sum(
                        ((main_i * a - ref_i) ** 2 for main_i, ref_i in zip(main_l, ref_l))
                       ) / len(main_l)
                   )
    main_l = df.to_dict(orient='list')[main] # SORT
    ref_l = df.to_dict(orient='list')[ref] # SORT
    res = minimize(myrmsd, [1.0])
    return res['x'][0]


def get_coeff_derivative(df, main, ref):
    NPOINTS = 5
    main_l = df.to_dict(orient='list')[main]
    ref_l = df.to_dict(orient='list')[ref]
    main_l = main_l[len(main_l) - NPOINTS:len(main_l)]
    ref_l = ref_l[len(ref_l) - NPOINTS:len(ref_l)]
    
    x = list(range(NPOINTS))
    main_slope = scipy.stats.linregress(x, main_l).slope
    ref_slope = scipy.stats.linregress(x, ref_l).slope
    if main_slope == 0.0:
        return None
    else:
        return ref_slope/main_slope

def prepare_scale_df(csvname):
    df = pd.read_csv(csvname, sep=';')
    bases = list(df['Basis'].unique())
    molecules = list(df['Molecule'].unique())
    ENERGY_TYPES = ['ScfEner', 'DF_total', 'E2_total', 'DF_sum', 'E2_sum', 'Del_total']

    df_parts = []
    scaledf = {
                  'Molecule': [],
                  'Basis': [],
                  'NBasis': [],
                  'EnergyType': [],
                  'Factor': [],
              }
    for basis in bases:
        for molecule in molecules:
            subdf = pd.DataFrame(df[(df['Basis'] == basis) & (df['Molecule'] == molecule)])
            nbasis = subdf.iloc[0]['NBasis']
            
            if molecule == 'ethane':
                for ener_type in ENERGY_TYPES:
                    if ener_type == 'ScfEner':
                        curmin = subdf[ener_type].min()
                        subdf[ener_type] = subdf[ener_type].transform(lambda e: e - curmin)
                    else:
                        curmax = subdf[ener_type].max()
                        subdf[ener_type] = subdf[ener_type].transform(lambda e: curmax - e)
            elif molecule == 'nh3hf' or molecule == 'oxygen':
                for ener_type in ENERGY_TYPES:
                    if ener_type == 'ScfEner':
                        curmax = subdf[ener_type].max()
                        subdf[ener_type] = subdf[ener_type].transform(lambda e: e - curmax)
                    else:
                        curmin = subdf[ener_type].min()
                        subdf[ener_type] = subdf[ener_type].transform(lambda e: curmin - e)
            
            df_parts.append(subdf)
            
            subdf.sort_values(by=['GeomVariable'], inplace=True)
            for ener_type in ENERGY_TYPES:
                if ener_type != 'ScfEner':
                    scaledf['Molecule'].append(molecule)
                    scaledf['Basis'].append(basis)
                    scaledf['NBasis'].append(nbasis)
                    scaledf['EnergyType'].append(ener_type)
                    if molecule == 'ethane':
                        scaledf['Factor'].append(get_coeff(subdf, ener_type, 'ScfEner'))
                    elif molecule == 'nh3hf' or molecule == 'oxygen':
                        scaledf['Factor'].append(get_coeff_derivative(subdf, ener_type, 'ScfEner'))

    # normdf = pd.concat(df_parts)
    # get_coeff(normdf, "DF_total", "ScfEner")
    return pd.DataFrame(scaledf)
