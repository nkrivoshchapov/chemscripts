from chemscripts.geom import Molecule, Fragment
from chemscripts import utils
import copy
import numpy as np

FILENAME = "ozoneattack_start"

if __name__ == "__main__":
    mymol = Molecule(FILENAME + '.sdf')
    mymol.from_xyz(*utils.parse_gjf(FILENAME + '.gjf'))
    dummyposition = (mymol.atom_xyz(4) + mymol.atom_xyz(5)) / 2
    dummy_idx = mymol.add_dummy(dummyposition)
    
    mymol.remove_bond(5, 27)
    mymol.remove_bond(4, 25)
    ozonepart = Fragment(mymol, 25)
    
    idx = 0
    for length in list(np.arange(0.0, 0.7, 0.05)):
        idx += 1
        newfrag = copy.deepcopy(ozonepart)
        newframe = newfrag.build_frame(26, dummy_idx, 25)
        newfrag.translate_x(newframe, length)
        newfrag.mol.remove_dummies()
        xyzs, syms = newfrag.mol.as_xyz()
        
        gjfname = "inpfiles/point_%d.gjf" % (idx)
        utils.write_gjf(xyzs, syms, 'preopt', gjfname,
                        subs={
                                'nproc': 12
                             })
    
    # mymol.save_sdf("check.sdf")