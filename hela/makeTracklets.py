from hela import *

import numpy as np

REAL = '<f16'
SSHORT = 'S20'
SMIN = 'S5'

def load_detections(fn, **kwargs):
    """Load detections file

    File format is assumed to be: ...?
    """
    intype = np.dtype(list(dict(
        MJD=REAL,
        RA=REAL,
        Dec=REAL,
        x=REAL,
        y=REAL,
        z=REAL,
        idstring=SSHORT,
        mag=float,
        band=SMIN,
        obscode=SMIN,
        index=int,
    ).items()))

    # x,y,z,index are not passed to the c++ layer, we just need dummy values
    # to make pybind11 play nice with the array as a input.
    #cols = [0,2,5,7,25,31,37,9,10,11,36]
    cols = [2,5,7,9,10,11,0,31,25,37,36]
    data = np.loadtxt(fn, usecols=cols, dtype=intype, delimiter=',', skiprows=1, converters={36: lambda x: int(float(x)), },**kwargs)
    return data

def load_observatory(fn, **kwargs):
    """Load-in observatory file
    """
    intype = np.dtype(list(dict(
        obscode=SMIN,
        obslon=REAL,
        plxcos=REAL,
        plxsin=REAL,
    ).items()))

    cols = [0,1,2,3]
    data = np.loadtxt(fn, usecols=cols, dtype=intype, skiprows=1, **kwargs)
    return data

def load_earth_ephemerides(fn, **kwargs):
    """ Load Earth ephemerides downloaded from horizons """

    def jpl_file(fn):
        with open(fn) as fp:
            # skip the header
            while True:
                line = fp.readline()
                if line == "$$SOE\n": break
                if line == '': raise EOFError

            # yield lines until we hit $$EOE
            while True:
                line = fp.readline()
                if line == "$$EOE\n": break
                if line == '': raise EOFError
                yield line

    MJDOFF=2400000.5
    data = np.genfromtxt(
        jpl_file(fn),
        dtype=[
                ('MJD', REAL),
                ('X', REAL), ('Y', REAL), ('Z', REAL),
                ('VX', REAL), ('VY', REAL), ('VZ', REAL)
        ],
        usecols=[0, 2, 3, 4, 5, 6, 7],
        autostrip=True,
        delimiter=',',
        **kwargs
    )
    data["MJD"] -= MJDOFF
    return data

def load_images(fn, **kwarg):
    intype = np.dtype(list(dict(
        MJD=float,
        RA=float,
        Dec=float,
        obscode=SMIN
    ).items()))

    if fn != "":
        data = np.loadtxt(fn, dtype=intype, delimiter=' ', **kwargs)
        return data

    print("No image file. One will be made.")
    return data

if __name__ == "__main__":

    config = MakeTrackletsConfig()
    config.indetfile = "../tests/sol_month04a.csv"
    config.inimfile = ""
    config.earthfile = "../tests/Earth1day2020s_02a.txt"
    config.obscodefile = "../tests/ObsCodes.txt"

    obsv = load_observatory(config.obscodefile)
    dets = load_detections(config.indetfile)

    makeTracklets(config, obsv, dets)