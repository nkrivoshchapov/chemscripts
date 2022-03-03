import subprocess, time, os, glob
from myscripts.mylogging import createLogger


def nbo_to_idx(nbo, reorder_nbo):
    return reorder_nbo.index(nbo) + 1

CUBENPROC = 6
def generate_nbo_cube(logname, fchkname, cubename, nbo_indices, gdriver, logger): 
    # cubename is expected to have {nbo} at least in cases when len(nbo_indices) > 1
    logger.info("Processing logfile " + logname)
    reorder_nbo = []
    loglines = open("%s.log" % logname, "r").readlines()
    for line in loglines:
        if "Reordering of NBOs for storage:" in line:
            parts = line.replace("Reordering of NBOs for storage:", "").replace("\n", "").split()
            for part in parts:
                reorder_nbo.append(int(part))

    cubefiles = []
    for nbo in nbo_idx:
        if "{nbo}" in cubename:
            newcube = cubename.format(nbo=nbo)
        else:
            assert len(nbo_indices) == 1
            newcube = cubename
        cubefiles.append(newcube)
        logger.info("Generating cubefile %s" % newcube)
        with gdriver['todo_lock']:
            gdriver['todo_files'].append({
                                            'command': 'cubegen {nproc} MO={corr_idx} {fchk} {cube} 150'.format(
                                            nproc=CUBENPROC, corr_idx=nbo_to_idx(nbo, reorder_nbo), fchk=fchkname, cube=newcube),
                                            'wd': os.path.dirname(newcube),
                                            'nproc': CUBENPROC,
                                            'resfile': newcube,
                                         })
    wait_for_termination([gjfname], gdriver)


if __name__ == "__main__":
    logger = createLogger("Main")
    data = []
    csvlines = open("nbo_indices.csv", 'r').readlines()
    for line in csvlines:
        if line.startswith('#'):
            continue
        parts = line.replace('\n', '').split(';')
        name = parts[0]
        nbos = [int(c) for c in parts[1].replace('[', '').replace(']', '').split(',')]
        data.append({
            'Name': name,
            'Nbos':nbos,
        })
    for item in data:
        logger.info("Running calc for " + item['Name'])
        prepfiles(item['Name'], item['Nbos'], logger)
