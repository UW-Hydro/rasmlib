"""
share.py
"""

HOURSPERDAY = 24.
SECSPERHOUR = 3600.
MINSPERHOUR = 60.
MINSPERDAY = HOURSPERDAY * MINSPERHOUR
SECSPERDAY = HOURSPERDAY * SECSPERHOUR
MONTHSPERYEAR = 12

NCOFORMAT = "%%Y-%%m-%%d %%H:%%M:%%S"

import subprocess


def argsort(seq):
    """ Equivalent to numpy argsort """
    #http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
    #by ubuntu
    return sorted(range(len(seq)), key=seq.__getitem__)

def call_nco(fargs):
    """Thin subprocess wrapper to call nco and handle exceptions"""
    # flatten nexted lists if present
    fargs = flatten(fargs)

    returnval = subprocess.Popen(fargs, stdout=subprocess.PIPE)
    stdout, stderr = returnval.communicate()
    if returnval.returncode:
        print(stderr)
        raise Exception("%s returned %s", fargs, returnval.returncode)
    elif stderr:
        print('WARNING, %s wrote to stderr: %s', fargs[0], stderr)
    return stdout, stderr

def flatten(foo):
    for x in foo:
        if hasattr(x, '__iter__'):
            for y in flatten(x):
                yield y
        else:
            yield x