from .names import Names


def _parse_equation(eq):
    fixed_line = eq.replace('\n', '').replace('\r', '').replace(' ', '')
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
    return fixed_line, reagents, products


def equations_from_k(filename):
    eqs = []
    lines = open(filename, 'r').readlines()
    for line in lines:
        fixed_line, reagents, products = _parse_equation(line)
        eqs.append({
            'original_line': fixed_line,
            'reagents': reagents,
            'products': products,
        })
    return eqs


def equations_from_sheet(sheet):
    eqs = []
    eq_block = sheet.block(Names.EQ_BLOCK)
    for item in eq_block['data']:
        fixed_line, reagents, products = _parse_equation(item[Names.EQ_COL])
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
