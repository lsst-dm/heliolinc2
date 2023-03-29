import numpy as np
import numpy.lib.recfunctions as rfn
import sys

sys.path.insert(1, '/home/aheinze/CppCode')
import heliohypy

SSHORT = 'S20'
SMIN = 'S5'

intype = np.dtype(list(dict(
    MJD=float,
    RA=float,
    Dec=float,
    mag=float,
    sigmag=float,
    idstring=SSHORT,
    band=SMIN,
    obscode=SMIN,
    traillen=float,
    trailpa=float,
    sigalong=float,
    sigacross=float,
).items()))

hldet = np.dtype(list(dict(
    MJD=float,
    RA=float,
    Dec=float,
    mag=float,
    trail_len=float,
    trail_PA=float,
    sigmag=float,
    sig_across=float,
    sig_along=float,
    image=int,
    idstring=SSHORT,
    band=SMIN,
    obscode=SMIN,
    index=int,
).items()))

shortimage = np.dtype(list(dict(
    MJD=float,
    RA=float,
    Dec=float,
    obscode=SMIN,
).items()))

hlimage = np.dtype(list(dict(
    MJD=float,
    RA=float,
    Dec=float,
    obscode=SMIN,
    X = float,
    Y = float,
    Z = float,
    VX = float,
    VY = float,
    VZ = float,
    startind = int,
    endind = int,
).items()))

hlradhyp = np.dtype(list(dict(
    HelioRad=float,
    R_dot=float,
    R_dubdot=float,
).items()))

hlclust = np.dtype(list(dict(
    clusternum=int,
    posRMS = float,
    velRMS = float,
    totRMS = float,
    astromRMS = float,
    pairnum = int,
    timespan = float,
    uniquepoints = int,
    obsnights = int,
    metric = float,
    rating = SSHORT,
    heliohyp0 = float, 
    heliohyp1 = float,
    heliohyp2 = float,
    posX=float,
    posY=float,
    posZ=float,
    velX=float,
    velY=float,
    velZ=float,
    orbit_a=float,
    orbit_e=float,
    orbit_MJD=float,
    orbitX=float,
    orbitY=float,
    orbitZ=float,
    orbitVX=float,
    orbitVY=float,
    orbitVZ=float,
    orbit_eval_count=int;
).items()))


def celeproj(ra, dec, **kwargs):
    """ Project an RA, Dec position onto the Cartesian unit sphere """
    x=np.cos(np.pi*dec/180.0)*np.cos(np.pi*ra/180.0)
    y=np.cos(np.pi*dec/180.0)*np.sin(np.pi*ra/180.0)
    z=np.sin(np.pi*dec/180.0)
    return(np.array([x,y,z]))

def celedeproj(cartvec, **kwargs):
    """  Deproject a (potentially un-normalized) Cartesian x,y,z position onto RA, Dec """
    norm = np.sqrt(cartvec[0]**2 + cartvec[1]**2 + cartvec[2]**2)
    x = cartvec[0]/norm
    y = cartvec[1]/norm
    z = cartvec[2]/norm
    if(abs(z)<1.0): Dec = np.arcsin(z)*180.0/np.pi
    if(y==0.0 and x<0.0): RA = 180.0
    elif(y==0.0): RA=0.0
    elif(y>0.0): RA = 90.0 - np.arctan(x/y)*180.0/np.pi
    elif(y<0.0): RA = 270.0 - np.arctan(x/y)*180.0/np.pi
    else: print("ERROR: illogical case in celedeproj")
    return(np.array([RA,Dec]))
    #This function needs some error catching, which I don't know how to do in python.

def make_image_table(detection_table, **kwargs):
    """  Construct an image table from an input detection catalog with MJD, RA, Dec, and obscode """
    imct=0
    timetol = 2.0
    imarray = np.array([detection_table['MJD'][0],detection_table['RA'][0],detection_table['Dec'][0]],ndmin=2)
    print ("on image", imct)
    image = np.array([detection_table['MJD'][0],detection_table['RA'][0],detection_table['Dec'][0],detection_table['obscode'][0]],ndmin=2)
    for i in range(detection_table.size-1):
        a = (detection_table['MJD'][i+1] - detection_table['MJD'][i])*86400.0
        if(a>timetol or detection_table['obscode'][i+1]!=detection_table['obscode'][i]):
            imct += 1
            print("on image", imct, "at entry", i, "size is", int(imarray.size/3))
            x = np.zeros((int(imarray.size/3),3))
            for j in range(int(imarray.size/3)):
                x[j] = celeproj(imarray[j][1],imarray[j][2])
            cartcen = np.array([np.mean((x.T[0].min(),x.T[0].max())),np.mean((x.T[1].min(),x.T[1].max())),np.mean((x.T[2].min(),x.T[2].max()))])
            boresight = celedeproj(cartcen)
            print(boresight)
            a = np.array([detection_table['MJD'][i],boresight[0],boresight[1],detection_table['obscode'][i]],ndmin=2)
            if(imct == 1):
                image[0][1]=boresight[0]
                image[0][2]=boresight[1]
            else:
                b=np.append(image,a,axis=0)
                image = b
            del(imarray)
            imarray = np.array([detection_table['MJD'][i+1],detection_table['RA'][i+1],detection_table['Dec'][i+1]],ndmin=2)
        else:
            q = np.array([detection_table['MJD'][i+1],detection_table['RA'][i+1],detection_table['Dec'][i+1]],ndmin=2)
            r = np.append(imarray,q,axis=0)
            imarray = r
    # handle final image
    imct += 1
    print("on image", imct, "at entry", i, "size is", int(imarray.size/3))
    x = np.zeros((int(imarray.size/3),3))
    for j in range(int(imarray.size/3)):
        x[j] = celeproj(imarray[j][1],imarray[j][2])
    cartcen = np.array([np.mean((x.T[0].min(),x.T[0].max())),np.mean((x.T[1].min(),x.T[1].max())),np.mean((x.T[2].min(),x.T[2].max()))])
    boresight = celedeproj(cartcen)
    print(boresight)
    a = np.array([detection_table['MJD'][i],boresight[0],boresight[1],detection_table['obscode'][i]],ndmin=2)
    b=np.append(image,a,axis=0)
    image = b
    return(image)


def read_ObsCodes(obscode_file, **kwargs):
    """ Read an observatory code file from MPC, and return it as an array with obscode, Longitude, cos(lat), sin(lat) """
    file1 = open(obscode_file, 'r')
    lines = file1.readlines()
    c = []
    for i in range(len(lines)-1) :
        a = list((lines[i+1][0:3], lines[i+1][4:13], lines[i+1][13:21], lines[i+1][21:30]))
        c.append(a)

    obsarr = np.array(c)
    return(obsarr)


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
                ('MJD', float),
                ('X', float), ('Y', float), ('Z', float),
                ('VX', float), ('VY', float), ('VZ', float)
        ],
        usecols=[0, 2, 3, 4, 5, 6, 7],
        autostrip=True,
        delimiter=',',
        **kwargs
    )
    data["MJD"] -= MJDOFF
    return data


def image_add_observerpos(image, obsarr, earthpos, **kwargs):
    """ Given an input image file, calculate the observer's heliocentric position and velocity at the time of every image """
    
    b=np.array([],dtype=hlimage)
    for i in range(len(image)) :
        observatory=0
        for j in range(len(obsarr)) :
            if(image[i][3].decode('utf-8') == obsarr[j][0]) : observatory = obsarr[j]
        mjd = float(image[i][0])
        Long = float(observatory[1])
        pcos = float(observatory[2])
        psin = float(observatory[3])
        obsx = heliohypy.observer_vel(mjd,Long,pcos,psin,earthpos)
        a1 = 1
        a = np.array([a1],dtype=hlimage)
        a['MJD'][0] = image[i][0]
        a['RA'][0] = image[i][1]
        a['Dec'][0] = image[i][2]
        a['obscode'][0] = image[i][3]
        a['X'][0] = obsx[0]
        a['Y'][0] = obsx[1]
        a['Z'][0] = obsx[2]
        a['VX'][0] = obsx[3]
        a['VY'][0] = obsx[4]
        a['VZ'][0] = obsx[5]
        a['startind'][0]=0;
        a['endind'][0]=0;
        b=np.append(b,a)

    return(b)
