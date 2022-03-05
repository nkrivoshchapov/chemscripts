import numpy as np
import mcubes, sklearn.decomposition
import glob

from myscripts.utils import wait_for_termination


def nbo_to_idx(nbo, reorder_nbo):
    return reorder_nbo.index(nbo) + 1


CUBENPROC = 6
def generate_nbo_cube(logname, fchkname, cubename, nbo_indices, gdriver, logger, wait=False):
    # cubename is expected to have {nbo} at least in cases when len(nbo_indices) > 1
    logger.info("Processing logfile " + logname)
    reorder_nbo = []
    loglines = open(logname, "r").readlines()
    for line in loglines:
        if "Reordering of NBOs for storage:" in line:
            parts = line.replace("Reordering of NBOs for storage:", "").replace("\n", "").split()
            for part in parts:
                reorder_nbo.append(int(part))

    waittasks = []
    cubefiles = []
    for nbo in nbo_indices:
        if "{nbo}" in cubename:
            newcube = cubename.format(nbo=nbo)
        else:
            assert len(nbo_indices) == 1
            newcube = cubename
        cubefiles.append(newcube)
        logger.info("Generating cubefile %s" % newcube)
        command = 'cubegen {nproc} MO={corr_idx} {fchk} {cube} 150'.format(
            nproc=CUBENPROC, corr_idx=nbo_to_idx(nbo, reorder_nbo), fchk=fchkname, cube=newcube)
        waittasks.append(command)

        with gdriver['todo_lock']:
            gdriver['todo_files'].append({
                                            'command': command,
                                            'wd': '.', # os.path.dirname(newcube),
                                            'nproc': CUBENPROC,
                                            'resfile': newcube,
                                         })
    if wait:
        wait_for_termination(waittasks, gdriver)
    return waittasks


BOHR2A = 0.529177
def generate_isosurface(cubename, meshfile_template, ival):
    cubelines = open(cubename, 'r').readlines()[2:] # First two lines are comments
    fline_parts = cubelines[0].split()
    natoms = int(fline_parts[0])
    if natoms < 0:
        natoms = abs(natoms) # Idk what the negative sign means here
    assert int(fline_parts[4]) == 1
    origin = np.array([
                                float(fline_parts[1]),
                                float(fline_parts[2]),
                                float(fline_parts[3]),
                           ])
    
    x_parts = cubelines[1].split()
    y_parts = cubelines[2].split()
    z_parts = cubelines[3].split()

    frame_data = {
                   'x': {'N': int(x_parts[0]), 'step': float(x_parts[1]), 'current': 0},
                   'y': {'N': int(y_parts[0]), 'step': float(y_parts[2]), 'current': 0},
                   'z': {'N': int(z_parts[0]), 'step': float(z_parts[3]), 'current': 0},
                 }
    assert float(x_parts[2]) == 0.0 and float(x_parts[3]) == 0.0
    assert float(y_parts[1]) == 0.0 and float(y_parts[3]) == 0.0
    assert float(z_parts[1]) == 0.0 and float(z_parts[2]) == 0.0

    syms = []
    xyzs = []
    for line in cubelines[4:4 + natoms]:
        parts = line.split()
        syms.append(int(parts[0]))
        xyzs.append(BOHR2A * np.array([
                                           float(parts[2]),
                                           float(parts[3]),
                                           float(parts[4])
                                      ]))
    
    voldata = np.empty([
                           frame_data['x']['N'],
                           frame_data['y']['N'],
                           frame_data['z']['N'],
                       ])
    for line in cubelines[5 + natoms:]:
        parts = line.split()
        for part in parts:
            voldata[frame_data['x']['current']] \
                   [frame_data['y']['current']] \
                   [frame_data['z']['current']] = float(part)
            frame_data['z']['current'] += 1
            if frame_data['z']['current'] == frame_data['z']['N']:
                frame_data['z']['current'] = 0
                frame_data['y']['current'] += 1
                if frame_data['y']['current'] == frame_data['y']['N']:
                    frame_data['y']['current'] = 0
                    frame_data['x']['current'] += 1

    for mode in [{'name': 'plus', 'ival': ival}, {'name': 'minus', 'ival': -ival}]:
        vertices, triangles = mcubes.marching_cubes(voldata, mode['ival'])
        for v in vertices:
            v[0] = (origin[0] + v[0] * frame_data['x']['step']) * BOHR2A
            v[1] = (origin[1] + v[1] * frame_data['y']['step']) * BOHR2A
            v[2] = (origin[2] + v[2] * frame_data['z']['step']) * BOHR2A
        edges = []
        for t in triangles:
            if {t[0], t[1]} not in edges:
                edges.append({t[0], t[1]})
            if {t[0], t[2]} not in edges:
                edges.append({t[0], t[2]})
            if {t[1], t[2]} not in edges:
                edges.append({t[1], t[2]})
        for i in range(len(edges)):
            edges[i] = tuple(edges[i])
        edges = np.array(edges).astype('int32')
        triangles = triangles.astype('int32')
        with open(meshfile_template.format(sign=mode['name'], type='vertices'), 'wb') as f:
            np.save(f, vertices)
        with open(meshfile_template.format(sign=mode['name'], type='edges'), 'wb') as f:
            np.save(f, edges)
        with open(meshfile_template.format(sign=mode['name'], type='triangles'), 'wb') as f:
            np.save(f, triangles)


def get_nbo_direction(surfnames):
    vertices = []
    for surfname in glob.glob(surfnames.format(type='vertices')):
        curvertices = list(np.load(surfname.format(type='vertices')))
        vertices += curvertices

    center = np.zeros(3)
    for v in vertices:
        center += v
    center /= len(center)

    for item in vertices:
        item -= center
    pca = sklearn.decomposition.PCA(3)
    pca.fit(vertices)
    return pca.components_[0]
