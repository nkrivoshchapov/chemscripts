
def parse_k(filename):
    eqs = []
    lines = open(filename, 'r').readlines()
    for line in lines:
        fixed_line = line.replace('\n', '').replace('\r', '').replace(' ', '')
        parts = fixed_line.split('<->')
        reagents = parts[0].split('+')
        for i, item in enumerate(reagents):
            if '*' in item:
                parts = item.split('*')
                n = int(parts[0])
                molname = '*'.join(parts[1:])
            else:
                n = 1
                molname = item
            reagents[i] = {'name': molname, 'n': n}
        products = parts[1].split('+')
        for i, item in enumerate(products):
            if '*' in item:
                parts = item.split('*')
                n = int(parts[0])
                molname = '*'.join(parts[1:])
            else:
                n = 1
                molname = item
            products[i] = {'name': molname, 'n': n}
        eqs.append({
            'original_line': fixed_line,
            'reagents': reagents,
            'products': products,
        })
    return eqs


def molecules_from_equations(eqs):
    mols = []
    for eq in eqs:
        for mol in eq['reagents'] + eq['products']:
            if mol['name'] not in mols:
                mols.append(mol['name'])
    return mols
