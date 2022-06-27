# This file is part of the LSST Solar System Processing lsst_dm_ssp.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


"""
nbody

LSST Solar System Processing

N-body propagation of asteroid orbits based on JPL DExxx planetary ephemeris

Implementation: Python 3.6, S. Eggl 20191115
"""


import numpy as np
import spiceypy as sp
from scipy.integrate import odeint

############################################
# Global Variables
###########################################
# list of spice kernels. DE4xx planetary ephemeris must come first.
spkFiles=['de430_1850-2150.bsp',
          'ast343de430.bsp']


############################################
# MODULE SPECIFIC EXCEPTION
###########################################

class Error(Exception):
    """Module specific exception."""
    pass


############################################
# Functions
###########################################

def getDEpos(target,t,frame='ECLIPJ2000'):
    return sp.spkezp(target,t,frame,'NONE',0)[0]

def getDEstate(target,t,frame='ECLIPJ2000'):
    return sp.spkez(target,t,frame,'NONE',0)[0]

try:
    getDEstateVec=np.vectorize(getDEstate)
    getDEposVec=np.vectorize(getDEpos)
except:
    raise Exception('Error: Could not vectorize getDEstate and getDEpos functions.')

def loadSpiceKernels(spkfiles=spkFiles):
    """ Load spice kernels. DE4xx planetary ephemeris must come first.
    
    Parameters:
    -----------
    spkfiles ... path to spice kernel files
    
    Returns:
    --------
    none
    """
    for f in spkfiles:
        try:
            sp.furnsh(f)
        except:
            raise Exception('Error in loadSpiceKernels: could not load spice file:', f)
    return

try:
    loadSpiceKernels(spkFiles)
except:
    raise Exception('Error: Could not load SPICE Kernels. Please check file path.')
                    
def readEphemerisHeader(filename):
    """Read JPL Digital Ephemeris DExxx header. 
    
    Parameters:
    -----------
    filename ... path to DExxx.bsp file
    
    Returns:
    --------
    header ... list of strings containing the header
    """
    with open(filename, 'r', encoding = "ISO-8859-1") as file:
        all_data = file.read(55000)
    header=(all_data.replace('\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00','')).split('\x00')
    header2=list(filter(None, header))[8:]
    return header2

def getEphemerisParameters(filename):
    """Get parameters from JPL DExxx ephemeris file.
    Parameters:
    -----------
    filename ... path to JPL DExxx file
    
    Returns:
    --------
    au2km     ... astronomical unit in kilometers
    clight    ... light speed
    GMlist    ... list of gravitational constant * mass for ephemeris bodies
    GMdict    ... dictionary which GM corresponds to which object
    Ephemdict ... DExxx numbers for spiceypy state calls
    ppndict   ... Consider Parametrized Post Newtonian terms for bodies where ppndict == True
    """
    header=readEphemerisHeader(filename)

    # G*mass for the sun, planets and the moon [AU**3/DAY**2] 
    GMs = []
    for i in range(11):
        GMs.append(float(header[16+i].split('   ')[2].replace("D", "E")))

    # G*mass for perturbing asteroids
    for i in range(343):
        GMs.append(float(header[30+i].split('  ')[2].replace("D", "E")))
    
    GMlist=GMs  

    # AU in kilometers
    au2km = float(header[381].split()[1])
    # light speed in km/sec
    clight = float(header[378].split()[1])
    
    # Solar System planets + Big16 perturbing asteroids
    Ephemdict={'Sun':10, 'Mercury':1, 'Venus':2, 'Earth':399, 'Moon':301, 
            'Mars':4, 'Jupiter':5,  'Saturn':6, 'Uranus':7, 'Neptune':8,
            'Pluto':9, 'Ceres':2000001, 'Pallas':2000002, 'Juno':2000003,
            'Vesta':2000004, 'Hebe':2000006, 'Iris':2000007, 'Hygiea':2000010,
            'Eunomia':2000015, 'Psyche':2000016, 'Amphitrite':2000029,
            'Europa':2000052, 'Cybele': 2000065, 'Sylvia':2000087,
            'Thisbe':2000088, 'Eros':2000433, 'Davida':2000511,    
             'Interamnia':2000704 }
    
#     asteroids=header[30:343]
#     i=0
#     for s in asteroids:
#         if('MA0704' in s):
#                print(i+11)    
#         i=i+1
    
    GMdict={'Sun':9, 'Mercury':0, 'Venus':1, 'Earth':2, 'Moon':10, 
            'Mars':3, 'Jupiter':4,  'Saturn':5, 'Uranus':6, 'Neptune':7,
            'Pluto':8, 'Ceres':11, 'Pallas':12, 'Juno':13,
            'Vesta':14, 'Hebe':16, 'Iris':17, 'Hygiea':20,
            'Eunomia':25, 'Psyche':26, 'Amphitrite':39,
            'Europa':61, 'Cybele': 71, 'Sylvia':90,
            'Thisbe':91, 'Eros':251, 'Davida':276,    
            'Interamnia':316}
    
    
    # Order contributions by acceleration mangitude
    Ephemorder=['Eros','Interamnia','Davida','Thisbe','Sylvia','Cybele',
                'Europa','Amphitrite','Psyche','Eunomia','Hygiea',
                'Iris','Hebe','Vesta','Juno','Pallas','Ceres',
                'Pluto','Mercury','Uranus','Neptune',
                'Venus','Mars','Moon','Earth','Saturn','Jupiter','Sun']
    
    ppndict=dict(zip(Ephemorder,[False for i in range(len(Ephemorder))]))
    ppndict['Sun']=True
    
    return au2km, clight, GMlist, GMdict, Ephemdict, Ephemorder, ppndict

try: 
    [au, clight, GMlist, GMdict, Ephemdict, Ephemorder, ppndict] = getEphemerisParameters(spkFiles[0])
except:
    raise Exception('Error: Could not acquire Ephemeris parameters.' 
                    'Please check file path for JPL DE ephemeris files.')
   

def mjd2et(MJD):
    """Converts modified Julian Date to JPL NAIF SPICE Ephemeris time.
    Only valid for TDB timescales.
    
    Parameters:
    -----------
    MJD ... modified Julian Day
    
    Returns:
    --------
    ET ... ephemeris time (ephemeris seconds beyond epoch J2000)

    """
    ET=(MJD+2400000.5-2451545.0)*86400
    return ET

############################################
# Force / acceleration calculation
###########################################


def acc_newton(y,t,GMlist,GMdict,Ephemdict,Ephemorder,au):
    """Newtonian accelerations for massless body in the Solar System.
    
    Parameters:
    -----------
    y          ... state (positions and velocities)
    t          ... epochs
    GMlist     ... List of G*mass for all perturbing bodies
    GMdict     ... Dictionary linking the name of perturbers to their GM id
    Ephemdict  ... Dictionary linking the name of perturbers to their Ephemeris id
    Ephemorder ... Order in which the accelerations should be summed. 
    au         ... astronomical unit [km]
    
    Returns:
    -------
    dydt       ... derivative of the state vector (dx/dt, dv/dt)
    """
    
    # gather positions of massive bodies
    n = len(Ephemorder)
    acc = np.zeros((n,3))
    dydt = np.zeros(6) 
    i=0
    
    for obj in Ephemorder:
        #pos = getDEposVec(Ephemdict[obj],mjd2et(t))/au
        pos = getDEpos(Ephemdict[obj],mjd2et(t))/au
        rel = pos[0:3]-y[0:3]
        d = np.linalg.norm(rel)
        GM=GMlist[GMdict[obj]]
        acc[i,:] = GM*rel/d**3
        i = i + 1
        
    dydt[0:3] = y[3:6]
    dydt[3:6] = np.sum(acc,axis=0)
    
#     print('y')
#     print(y)
#     print('t')
#     print(t)
#     print('acc')
#     print(acc)
#     print('dydt')
#     print(dydt)
    return dydt


def acc_ppn(y,t,GMlist,GMdict,Ephemdict,Ephemorder,ppndict,au,clight):
    """Parameterized Post Newtonian accelerations for massless body in the Solar System.
    
    Parameters:
    -----------
    y          ... state (positions and velocities)
    t          ... epochs
    GMlist     ... List of G*mass for all perturbing bodies
    GMdict     ... Dictionary linking the name of perturbers to their GM id
    Ephemdict  ... Dictionary linking the name of perturbers to their Ephemeris id
    Ephemorder ... Order in which the accelerations should be summed. 
    ppndict    ... Dictionary, which perturbers should have Post Newtonian terms added to their accelerations?
    au         ... astronomical unit [km]
    clight     ... speed of light [km/s]
    
    Returns:
    -------
    dydt       ... derivative of the state vector (dx/dt, dv/dt)
    
    
    """
    
    n = len(Ephemorder)
    acc = np.zeros((n,3))
    dydt = np.zeros(6) 
    kmps2aupd=86400./au
    c2 = (clight*kmps2aupd)**2
    
    i=0
    
    for obj in Ephemorder:
        #rvkm = sp.spkez(Ephemdict[obj],mjd2et(t),'ECLIPJ2000','NONE',0)[0]
        #rvkm = getDEstateVec(Ephemdict[obj],mjd2et(t))
        rvkm = getDEstate(Ephemdict[obj],mjd2et(t))
        r = rvkm[0:3]/au
        v = rvkm[3:6]*kmps2aupd
        rrel = y[0:3]-r[0:3]
        vrel = y[3:6]-v[0:3]
        d = np.linalg.norm(rrel)
        d3 = d**3
        GM = GMlist[GMdict[obj]]
        if(ppndict[obj]):
            acc[i,:] = GM/d3*(-rrel[:] + 
                          ((4.*GM/d-np.dot(v,v))*rrel[:] +
                            4.*np.dot(rrel,vrel)*vrel)/c2)
        else:
            acc[i,:] = GM/d3*(-rrel[:])
        #print(acc[i,:]-GM/d3*rrel[:])    
#         print(GM/d3*(((4.*GM/d-np.dot(v,v))*rrel[:] +
#         4.*np.dot(rrel,vrel)*vrel)/c2))
        i = i + 1
        
    dydt[0:3] = y[3:6]
    dydt[3:6] = np.sum(acc,axis=0)
    
#     print('y')
#     print(y)
#     print('t')
#     print(t)
#     print('acc')
#     print(acc)
#     print('dydt')
#     print(dydt)
    return dydt


def propagate_nbody(state0, epoch, tp, ppn=True):
    """ Propagate state through nbody-code using JPL DE430 ephemeris including BIG16 perturbing asteroids.
    
    Parameters:
    -----------
    state0       ... initial state (barycentric ecliptic frame) [au, au/day]
    epoch   ... epoch (TDB)of initial conditions [MJD]
    tp      ... list of target epochs (TDB) [MJD]
    
    Keyword Arguments:
    ------------------
    ppn     ... Parametrized Post Newtonian equations for the sun? (True, False)

    Returns:
    --------
        t       ... epochs of output (TDB) [MJD]
    statep      ... propagated state (barycentric ecliptic frame) [au, au/day]
    """

    #loadSpiceKernels(spkFiles)

    t = np.insert(tp,0,epoch)
    
    if(ppn):
        statep = odeint(acc_ppn, state0, t, args=(GMlist,GMdict,Ephemdict,
                                                  Ephemorder,ppndict,au,clight),rtol=4.E-10, atol=4.E-10)

    else:
        statep = odeint(acc_newton, state0, t, args=(GMlist,GMdict,
                                                 Ephemdict,Ephemorder,au), rtol=4.E-10, atol=4.E-10)

    return t, statep
    
