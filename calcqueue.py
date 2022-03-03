import multiprocessing, threading, subprocess, time, os, ntpath
from . import utils

MAXPROC = 22
MAXTRIES = 3
SLEEP_DURATION = 0.1

def get_next_task(todo_files):
    if len(todo_files) == 0:
        return None
    elif todo_files[0] == "Done":
        return ("Done", 0)
    else:
        return (todo_files[0], utils.getnproc(todo_files[0]))


def gauss_driver(todo_files, todo_lock, done_files, done_lock):
    procs = []
    nfiles = 0
    termination_mode = False
    ntries = {}
    occupied_proc = 0
    
    with todo_lock:
        next_task = get_next_task(todo_files) # Tuple (filename, nproc)
    while nfiles > 0 or len(procs) > 0 or not termination_mode:
        for i in reversed(range(len(procs))):
            if not procs[i]['proc'].poll() == None:
                if utils.is_normal_termination(procs[i]['logfile']):
                    # print("Normal termination of " + procs[i]['logfile'])
                    with done_lock:
                        done_files.append(procs[i]['inpfile'])
                else:
                    print("Not normal termination of " + procs[i]['logfile'])
                    todo_files.append(procs[i]['inpfile'])
                    if procs[i]['inpfile'] not in ntries:
                        ntries[procs[i]['inpfile']] = 1
                    else:
                        ntries[procs[i]['inpfile']] += 1
                occupied_proc -= procs[i]['nproc']
                del procs[i]

        if next_task is None:
            with todo_lock:
                if len(todo_files) > 0:
                    next_task = get_next_task(todo_files)

        while nfiles > 0 and next_task[1] <= MAXPROC - occupied_proc:
            with todo_lock:
                calc_file = todo_files.pop(0)
                curnproc = next_task[1]
                next_task = get_next_task(todo_files)

            if calc_file in ntries and ntries[calc_file] > MAXTRIES:
                continue

            if calc_file == "Done":
                print("Enabled termination mode")
                termination_mode = True
            else:
                tempwd = os.path.dirname(calc_file)
                mainwd = os.getcwd()
                os.chdir(tempwd)
                if calc_file.endswith('.gjf'):
                    procs.append({
                                    'proc': subprocess.Popen("rung " + ntpath.basename(calc_file), shell = True),
                                    'inpfile': calc_file,
                                    'logfile': calc_file.replace('.gjf', '.log'),
                                    'wd': tempwd,
                                    'nproc': curnproc
                                 })
                elif calc_file.endswith('.47'):
                    procs.append({
                                    'proc': subprocess.Popen("NBO6 " + ntpath.basename(calc_file), shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL),
                                    'inpfile': calc_file,
                                    'logfile': calc_file.replace('.47', '_NBO.out'),
                                    'wd': tempwd,
                                    'nproc': curnproc
                                 })
                os.chdir(mainwd)
                occupied_proc += curnproc
            with todo_lock:
                nfiles = len(todo_files)
        with todo_lock:
                nfiles = len(todo_files)
        time.sleep(SLEEP_DURATION)


def init_thread():
    todo_files = []
    todo_lock = multiprocessing.Lock()
    done_files = []
    done_lock = multiprocessing.Lock()
    thread = threading.Thread(target=gauss_driver, args=(todo_files, todo_lock, done_files, done_lock))
    thread.start()
    gdriver = {
                  'todo_lock': todo_lock,
                  'done_lock': done_lock,
                  'todo_files': todo_files,
                  'done_files': done_files
              }
    return gdriver
