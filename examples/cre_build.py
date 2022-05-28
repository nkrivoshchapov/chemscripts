import sys, glob, ntpath, os
import chemscripts.excelutils as xl
from chemscripts.energyutils import prepare_sheet, get_energy_sheet
from chemscripts.geom import RmsdPool
from chemscripts import utils

class _Names(object):
    MAIN_BLOCK = "Main"
    XYZNAME_COL = "Filename"
    FIRSTATOM_COL = "First atom"
    SECATOM_COL = "Second atom"
    MAXDIST_COL = "Maximal distance"

# RAWXYZ_DIR = "./raw_cis"
# PREPROC_DIR = "./preproc_cis"
RAWXYZ_DIR = "./raw_trans"
PREPROC_DIR = "./preproc_trans"
MAXENER = 10 # kcal/mol
MAXRMSD = 0.2

if __name__ == "__main__":
    if len(sys.argv) == 1: # No args
        excelsheet = xl.ExcelSheet()
        excelsheet.add_block(blockname=_Names.MAIN_BLOCK,
                             cols=(_Names.XYZNAME_COL, _Names.FIRSTATOM_COL, _Names.SECATOM_COL, _Names.MAXDIST_COL))
        
        for xyzfile in glob.glob(os.path.join(RAWXYZ_DIR, "*.xyz")):
            justname = ntpath.basename(xyzfile)
            newitem = {_Names.XYZNAME_COL: justname}
            excelsheet.add_row(blockname=_Names.MAIN_BLOCK, data=newitem)
        excelsheet.save_xlsx("template.xlsx")
    elif sys.argv[1].endswith(".xlsx"):
        utils.checkout_directory(PREPROC_DIR)
        nstructs = []
        
        excelsheet = xl.ExcelSheet()
        excelsheet.read_xlsx(sys.argv[1], get_values=True)
        myblock = excelsheet.block(_Names.MAIN_BLOCK)
        for item in myblock['data']:
            xyzfile = item[_Names.XYZNAME_COL]
            print("Processing " + xyzfile)
            first_atom = item[_Names.FIRSTATOM_COL]
            second_atom = item[_Names.SECATOM_COL]
            maxdist = item[_Names.MAXDIST_COL]
            
            xyzpath = os.path.join(RAWXYZ_DIR, xyzfile)
            mypool = RmsdPool()
            mypool.include(xyzpath, energy_func=lambda comment_line: float(comment_line.strip()))
            nstructs.append({'filename': xyzfile, 
                             'before': len(mypool.structures)
                            })
            
            mypool.energy_filter(MAXENER)
            if first_atom is not None:
                mypool.distance_filter(first_atom, second_atom, keep_if=lambda dist: dist < maxdist)
            else:
                print("CAUTION! Constraints weren't given for " + xyzfile)
            newname = xyzfile.replace('done', 'preproc')
            mypool.save(os.path.join(PREPROC_DIR, newname))
            nstructs[len(nstructs) - 1]['after'] = len(mypool.structures)
        
        print("-- Summary --")
        before_total = 0
        after_total = 0
        for item in nstructs:
            print("{xyzfile} contains {after}({before}) geometries".format(xyzfile=item['filename'],
                        after=item['after'],
                        before=item['before'],))
            before_total += item['before']
            after_total += item['after']
        print("Total: {after}({before})".format(after=after_total, before=before_total))
    elif len(sys.argv) == 2 and "rmsd" in sys.argv:
        mainpool = RmsdPool()
        for xyzpath in glob.glob(os.path.join(PREPROC_DIR, "*.xyz")):
            xyzfile = ntpath.basename(xyzpath)
            mainpool.include(xyzpath, energy_func=lambda comment_line: float(comment_line.strip()))
            print("{file} was successfully parsed. Current pool size = {n}".format(file=xyzfile,
                    n=len(mainpool.structures)))
        mainpool.energy_filter(MAXENER)
        print("Pool size after energy filtering = " + str(len(mainpool.structures)))
        mainpool.filter_rmsd(MAXRMSD)
        print("Pool size after RMSD filtering = " + str(len(mainpool.structures)))
        mainpool.save("final_confs.xyz")
        