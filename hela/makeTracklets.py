from hela import *

import numpy as np
import numpy.lib.recfunctions as rfn
import pandas as pd

REAL = '<f16'
SSHORT = 'S20'
SMIN = 'S5'

def load_detections_dataframe(fn, **kwargs):
    """Load detections from a csv through pandas to numpy

    *Testing method to convert from dataframe to structured array*
    """
    intype = {
        "MJD":REAL,
        "RA":float,
        "Dec":float,
        "idstring":SSHORT,
        "mag":float,
        "band":SMIN,
        "obscode":SMIN,
        "x":REAL,
        "y":REAL,
        "z":REAL,
        "index":int
    }

    df = pd.read_csv(fn, index_col=0, **kwargs)
    data = df.to_records(index=False, column_dtypes=intype)
    return data

def load_detections(fn, **kwargs):
    """Load detections file

    File format is assumed to be: ...?
    """
    intype = np.dtype(list(dict(
        MJD=REAL,
        RA=float,
        Dec=float,
        idstring=SSHORT,
        mag=float,
        band=SMIN,
        obscode=SMIN,
    ).items()))

    defaults_dtype = np.dtype(list(dict(
        x=REAL,
        y=REAL,
        z=REAL,
        index=int
    ).items()))

    cols = [2,5,7,0,31,25,37]
    data = np.loadtxt(fn, usecols=cols, dtype=intype, delimiter=',', skiprows=1, **kwargs)
    defaults = np.zeros(len(data), dtype=defaults_dtype)
    return rfn.merge_arrays((data, defaults), flatten=True)

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
        obscode=SMIN,
        startind=int,
        endind=int
    ).items()))

    if fn != "":
        data = np.loadtxt(fn, dtype=intype, delimiter=' ', **kwargs)
        return data

    print("No image file. One will be made.")
    return np.array([], dtype=intype)

if __name__ == "__main__":

    config = MakeTrackletsConfig()
    config.indetfile = "../tests/sol_month04a.csv"
    config.inimfile = ""
    config.earthfile = "../tests/Earth1day2020s_02a.txt"
    config.obscodefile = "../tests/ObsCodes.txt"

    obsv = load_observatory(config.obscodefile)
    dets = load_detections(config.indetfile)
    imgs = load_images(config.inimfile)

    # Need proper input files
    # earthsv = load_earth_ephemerides(config.earthfile)

    outs = makeTracklets(config, obsv, dets, imgs)

    df_imgs = pd.DataFrame(outs[0])
    df_pairs = pd.DataFrame(outs[1])

    df_imgs.to_csv("imgs_test.csv")
    df_pairs.to_csv("pairdets_test.csv")

