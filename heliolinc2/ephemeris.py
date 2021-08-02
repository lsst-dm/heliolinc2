# This file is part of the LSST Solar System Processing lsstssp.
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
ephemeris

LSST Solar System Processing routines for
calculating ephemerides of moving objects


Implementation: S. Eggl 20120308
"""

# Accelerators
import numpy as np
import numba

# External API's
from astroquery.jplhorizons import Horizons

#Interpolation
import scipy.interpolate as spi

# time scale transforms
from astropy.time import Time

# Accelerated vector operations
from . import vector as vec

# Constants such as the speed of light and GM
from . import constants as cn

# Coordinate transfroms
from . import transforms as tr

# State Propagation routines
from . import propagate as pr


__all__ = ['getObserverStates', 'getObserverInterpolant', 
           'getTargetStates', 'targetStatesFromHorizons',
           'icrf2ephemeris', 
           'topocentric2ephemeris', 'state2ephemeris', 'radecResiduals']





############################################
# MODULE SPECIFIC EXCEPTION
###########################################
class Error(Exception):
    """Vector module specific exception."""
    
    pass

###########################################
# OBSERVER STATES 
###########################################

def getObserverStates(observation_epochs,origin='SSB',observer_location='I11',ephemeris_dt='1h',frame='ecliptic'):
    """Produce sun-observer state vectors at observation epochs.
    
    Parameters:
    -----------
    observation_epochs         ... Numpy array of observation epochs [JD] 
    observer_location          ... Horizons identifyer of observer location, e.g. 'I11'
    ephemeris_dt               ... Time step for ephemeris query. 
                                   Typically 1h since the actual times will be interpolated later.
    
    Returns:
    --------
    observer_positions         ... Heliocentric observer positions at observation epochs in [au].
    
    
    External Function Requirements:
    -------------------------------

    # Interpolation
    import scipy.interpolate as spi
    
    # time transform
    mjd2jd                         ... change modified Julian date to Julian date, timescale TDB)
    
    # NASA JPL HORIZONS API call wrapper
    observerStatesFromHorizons  ... Wrapper function for JPL Horizons state query via astropy
    
    """

    tmin = np.min(observation_epochs)
    tmax = np.max(observation_epochs)
    
    #Start and stop times of the survey
    tstart = tmin-1.
    tstop = tmax+1.

    epochs = np.unique(observation_epochs)

    ir = getObserverInterpolant(tstart,tstop,origin=origin,
                                observer_location=observer_location,
                                ephemeris_dt=ephemeris_dt, frame=frame)
    
    # Interpolate heliocentric observer positions to the actual observation epochs
    observer_positions = ir(observation_epochs)
    # Interpolate heliocentric observer velocities to the actual observation epochs
    observer_velocities = ir.derivative(nu=1)(observation_epochs)
    
    return observer_positions.T, observer_velocities.T


def getObserverInterpolant(tmin,tmax,origin='SSB',observer_location='I11',ephemeris_dt='1h', frame='ecliptic'):
    """Produce sun-observer state vectors at observation epochs.
    
    Parameters:
    -----------
    tmin,tmax                  ... float, float, start and stop time for interpolant [MJD] 
    origin                     ... str, origin of the coordinate system (e.g. 'Sun' or 'SSB' =Solar System Barycenter)
    observer_location          ... str, Horizons identifyer of observer location, e.g. 'I11'
    ephemeris_dt               ... float, Time step for ephemeris query. 
                                   Typically 1h since the actual times will be interpolated later.
    frame                      ... str, Coordinate system reference frame: 'ecliptic' or 'ICRF'
    
    Returns:
    --------
    ipos                       ...  scipy interpolnt [au] Heliocentric observer position interpolant 
                                    as function of time [MJD].
                                    velocities can be interpolated via ipos.derivative(nu=1)
   
    
    External Function Requirements:
    -------------------------------

    # Interpolation
    import scipy.interpolate as spi
    
    # time transform
    mjd2jd                         ... change modified Julian date to Julian date, timescale TDB)
    
    # NASA JPL HORIZONS API call wrapper
    observerStatesFromHorizons  ... Wrapper function for JPL Horizons state query via astropy
    
    """
    
    tminjd = tr.mjd2jd(tmin)
    tmaxjd = tr.mjd2jd(tmax)
    
    #Start and stop times of the survey
    tstart = 'JD'+str(tminjd-1.)
    tstop = 'JD'+str(tmaxjd+1.)
    
    try:
        # Get observer locations (caution: choose the right plane of reference and direction of the vectors!)
        # check query by copy/pasting the output of print(observer_sun.uri) into a webbrowser if there are problems.
    
        
        if(origin == 'SSB' or origin == '@0'):
        
            observer_origin = Horizons(id='Sun', location=observer_location, id_type='majorbody', 
                          epochs = {'start':tstart, 'stop':tstop,
                          'step':ephemeris_dt})

            origin_barycenter = Horizons(id='Sun', location='@0', id_type='majorbody', 
                          epochs = {'start':tstart, 'stop':tstop,
                          'step':ephemeris_dt})


            if(frame == 'ecliptic'):
                oo = observer_origin.vectors(refplane='ecliptic')  
                ob = origin_barycenter.vectors(refplane='ecliptic')
                
            elif(frame == 'ICRF' or frame == 'J2000' or frame == 'earth' or frame == 'icrf'):
                oo = observer_origin.vectors(refplane='earth')  
                ob = origin_barycenter.vectors(refplane='earth')   
            else:
                raise Exception('Error: requested frame unknown.')
            
            observer_xyz = (-1)*(np.array([oo['x'],oo['y'],oo['z']]).astype('float') +
                                 np.array([ob['x'],ob['y'],ob['z']]).astype('float'))
            
            observer_vxyz =(-1)*(np.array([oo['vx'],oo['vy'],oo['vz']]).astype('float') +
                                 np.array([ob['vx'],ob['vy'],ob['vz']]).astype('float'))
            
            observer_jd = np.array(oo['datetime_jd']).astype('float')
        
        else:
            
            observer_origin = Horizons(id=origin, location=observer_location, id_type='majorbody', 
                          epochs = {'start':tstart, 'stop':tstop,
                          'step':ephemeris_dt})

            if(frame == 'ecliptic'):
                obs = observer_origin.vectors(refplane='ecliptic') 
            elif(frame == 'ICRF' or frame == 'J2000' or frame == 'earth' or frame == 'icrf'):
                obs = observer_orignin.vectors(refplane='earth') 
            else:
                raise Exception('Error: requested frame unknown.')
        
                    #We need the sun-observer vector not the observer-sun vector
            observer_xyz = (-1)*np.array([obs['x'],obs['y'],obs['z']]).astype('float')
            observer_vxyz =(-1)*np.array([obs['vx'],obs['vy'],obs['vz']]).astype('float')
            observer_jd = np.array(obs['datetime_jd']).astype('float')
        
    except:
        print("Error: potential online ephemeris query failure. Make sure internet connectivity is available.")
        raise
        
    
    observer_mjd = tr.jd2mjd(observer_jd)
        
    # Interpolate heliocentric observer positions to the actual observation epochs
    ipos = spi.CubicHermiteSpline(observer_mjd, observer_xyz,observer_vxyz, axis=1, extrapolate=None)
    
    # Interpolate heliocentric observer velocities to the actual observation epochs
    # ivel=ipos.derivative(nu=1)
    
    return ipos


def getTargetStates(observation_epochs,target_id='Duende', observer_id='I11', ephemeris_dt='1h', frame='ecliptic'):
    """Produce sun-observer state vectors at observation epochs.
    
    Parameters:
    -----------
    observation_epochs         ... Numpy array of observation epochs [JD] 
    target_id                  ... Horizons identifier target, e.g. 'Ceres'
    observer_id                ... Horizons identifier observer, e.g. 'I11'
    ephemeris_dt               ... Time step for ephemeris query. 
                                   Typically 1h since the actual times will be interpolated later.
    frame                      ... coordinate frame (default 'ecliptic' or 'icrf')
    
    Returns:
    --------
    target_positions            ... Heliocentric target states at observation epochs in [au, au/day].
    target_velocities
    
    External Function Requirements:
    -------------------------------

    # Interpolation
    import scipy.interpolate as spi
    
    # time transform
    mjd2jd                         ... change modified Julian date to Julian date, timescale TDB)
    
    # NASA JPL HORIZONS API call wrapper
    targetStatesFromHorizons  ... Wrapper function for JPL Horizons state query via astropy
    
    """

    tmin = tr.mjd2jd(np.min(observation_epochs))
    tmax = tr.mjd2jd(np.max(observation_epochs))
    
    #Start and stop times of the survey
    tstart = 'JD'+str(tmin-1.)
    tstop = 'JD'+str(tmax+1.)

    epochs = np.unique(observation_epochs)

    [target_mjd,target_xyz,target_vxyz] = targetStatesFromHorizons(target_id, observer_id,
                                                                  tstart,tstop, ephemeris_dt)
      
    # Interpolate heliocentric observer positions to the actual observation epochs
    ir = spi.CubicHermiteSpline(target_mjd, target_xyz, target_vxyz, axis=1, extrapolate=None)
    target_positions = ir(observation_epochs).T
    # Interpolate heliocentric observer velocities to the actual observation epochs
    dirdt=ir.derivative(nu=1)
    target_velocities = dirdt(observation_epochs).T
    
    return target_positions, target_velocities

def targetStatesFromHorizons(target_id, observer_id,
                             tstart, tstop, ephemeris_dt='1h', frame='ecliptic'):
    """Query JPL Horizons via astroquery to get sun-observer state vectors.
    
    Parameters:
    -----------
    target_id          ... Horizons identifier of target, e.g. 'Ceres'
    observer_id        ... Horizons identifier observer, e.g. 'I11'
    tstart             ... start time for ephemeris in Horizons format, e.g. 'JD2456789.5'
    tstop              ... end time for ephemeris in Horizons format, e.g. 'JD2456799.5'
    ephemeris_dt       ... Time step for ephemeris query. 
                           Typically 1h since the actual times will be interpolated later.
    frame              ... coordinate frame ('ecliptic' or 'icrf')
    
    Returns:
    --------
    observer_xyz       ... Heliocentric observer positions at gridded epochs in [au].
    observer_vxyz      ... Heliocentric observer velocities at gridded epochs in [au].
    observer_jd        ... Gridded ephemeris epochs (JD / TDB)
    
    External Function Requirements:
    -------------------------------
    # External API's
    from astroquery.jplhorizons import Horizons
    """
    try:
        # Get observer locations (caution: choose the right plane of reference and direction of the vectors!)
        # check query by copy/pasting the output of print(observer_sun.uri) into a webbrowser if there are problems.
        obs2target = Horizons(id=target_id, location=observer_id,
                      epochs = {'start':tstart, 'stop':tstop,
                      'step':ephemeris_dt})
        
        
        if(frame == 'ecliptic'):
            vec = obs2target.vectors(refplane='ecliptic') 
        elif(frame == 'ICRF' or frame == 'J2000' or frame == 'earth' or frame == 'icrf'):
            vec = obs2target.vectors(refplane='earth') 
        else:
            raise Exception('Error: requested frame unknown.')
        
        
        target_jd = np.array(vec['datetime_jd']).astype('float')
        
        #We need the sun-observer vector not the observer-sun vector
        target_xyz = np.array([vec['x'],vec['y'],vec['z']]).astype('float')
        target_vxyz = np.array([vec['vx'],vec['vy'],vec['vz']]).astype('float')
        
    except:
        print("Error: potential online ephemeris query failure. Make sure internet connectivity is available.")
        raise
          
    target_mjd = tr.jd2mjd(target_jd)
        
    return target_mjd, target_xyz, target_vxyz


def targetEphemerisFromHorizons(target_id, observer_location, 
                                tstart, tstop, ephemeris_dt='12h'):
    """Query JPL Horizons via astroquery to get sun-observer state vectors.
    
    Parameters:
    -----------
    target_id          ... Horizons identifyer target, e.g. 'Ceres'
    observer_location  ... Horizons identifyer of observer location, e.g. 'I11'
    tstart             ... start time for ephemeris in Horizons format, e.g. 'JD2456789.5'
    tstop              ... end time for ephemeris in Horizons format, e.g. 'JD2456799.5'
    ephemeris_dt       ... Time step for ephemeris query. 
                           Typically 1h since the actual times will be interpolated later.
    
    Returns:
    --------
    RA, DEC               ... Right Ascension and Declination of target wrt observer
    ephemeris_jd        ... Ephemeris epochs (JD / UTC)
    
    External Function Requirements:
    -------------------------------
    # External API's
    from astroquery.jplhorizons import Horizons
    """
    try:
        # Get observer locations (caution: choose the right plane of reference and direction of the vectors!)
        # check query by copy/pasting the output of print(observer_sun.uri) into a webbrowser if there are problems.
        
#         tmin = 'JD'+str(tstart)
#         tmax = 'JD'+str(tstop+1.)
    
        tmin = tstart
        tmax = tstop
        observer2target = Horizons(id=target_id, location=observer_location,
                      epochs = {'start':tmin, 'stop':tmax,
                      'step':ephemeris_dt})

        ra, dec = observer2target.ephemerides()['RA','DEC']
        jd = (observer2target.vectors())['datetime_jd']
       
        
    except:
        print("Error in observer_state_from_horizons: potential online ephemeris query failure.")
        raise
        
    return ephemeris_jd, ra, dec


###########################################
# EPHEMERIS CALCULATION
###########################################

def icrf2ephemeris(epoch, state, timescale_epoch='utc', 
                   timescale_state='tdb', time_format='mjd', deg=True, lttc=True):
    """Transform topocentric ICRF states to
    Right Ascension (RA) and Declination (DEC) observations.
    
    Parameters:
    -----------
    epoch             ... epoch of RADEC observation (astropy.time object)
    state             ... topocentric state vector of object (ICRF) [au] 
                      ... (given at epoch in TDB)
    epoch_timescale   ... time scale for epoch ('utc', 'tdb'), default utc
    state_timescale   ... time scale for state epoch ('utc', 'tdb'), default tdb
    deg               ... True (default): angles in degrees, False: angles in radians
    lttc              ... True (default): correct for light travel time \
                          (needs entire state including velocities)
    
    Returns:
    --------
    RA               ... Right Ascension [default deg]
    DEC              ... Declination [default deg]
    dRA/dt*cos(DEC)  ... sky plane motion in RA [default deg/day]
    dDEC/dt          ... sky plane motion in DEC [default deg/day]
    r                ... distance to object [au]
    dr/dt            ... line of sight velocity [au/day]

    """
    
    #numpy functions for math and vector operations
    multiply = np.multiply
    divide = np.divide
    add = np.add
    norm = np.linalg.norm
    dot = np.dot
    mod = np.mod
    arcsin = np.arcsin
    sqrt = np.sqrt
    cos = np.cos
    
    # other numpy functions and constants
    array = np.array
    hstack = np.hstack
    rad2deg = np.rad2deg
    
    #PI * 2
    PIX2=np.pi*2
    

    # Correct for time shift between UTC and TDB
    time=Time(epoch, scale=timescale_epoch, format=time_format)
    
    if(timescale_state != timescale_epoch):
        dt = ((Time(time, scale=timescale_epoch,format=time_format)).value.astype('float') - 
             (Time(time, scale=timescale_state,format=time_format)).value.astype('float'))
    else:
        dt = 0.
      
    if(state.ndim == 1):  
        # Check if we have the full state 
        if(len(state)!=6):
                   raise Exception('Error in icrf2radec: Full state including velocities \
                        needed for light travel time and timescale correcion.')           
           
        # Add Light travel time correction to TDB-UTC
        if(lttc):
            # determine state corrected for light travel time
            r0 = norm(state[0:3]) 
            dt = dt + r0/cn.CAUPD
        
        state_corrected = hstack([state[0:3] - dt*state[3:6] 
                                    # - dt**2*cnst.GM/r0**3*state[0:3], 
                                     , state[3:6]])
       
        # Calculate Right Ascension and Declination
        r1 = norm(state_corrected[0:3]) 
        
        rn = state_corrected[0:3]/r1 
        rdot = dot(rn, state_corrected[3:6])
        
        RA = mod(np.arctan2(rn[1], rn[0])+PIX2, PIX2)
        DEC = arcsin(rn[2])
        

        # Analytic Derivatives
        # dalpha/dt = (x dy/dt - y dx/dt)/(x^2+y^2)
        # ddelta/dt = (dz/dt r - z dr/dt)/(r^2 sqrt(1-z^2/r^2))
        dRAdt = (state_corrected[0]*state_corrected[4]-
                 state_corrected[1]*state_corrected[3]
                 )/(state_corrected[0]**2+state_corrected[1]**2)
        dDECdt = (state_corrected[5]*r1-state_corrected[2]*rdot)/(r1*
                  sqrt(r1**2-state_corrected[3]**2))
        
#         print('dRAdt,dDECdt, analytics') 
#         print([np.rad2deg(dRAdt)*np.cos(DEC),np.rad2deg(dDECdt)])
        
         # Finite Differences for derivatives
#         RA_later = np.mod(np.arctan2(state[1],state[0])+pix2,pix2)
#         DEC_later = np.arcsin(state[2]/r0)
#         dRAdt=(RA_later-RA)/dt
#         dDECdt=(DEC_later-DEC)/dt
       
#         print('dRAdt,dDECdt, finite diff')     
#         print([np.rad2deg(dRAdt)*np.cos(DEC),np.rad2deg(dDECdt)])
        
    else:
        if(state.shape[1]!=6):
                raise Exception('Error in icrf2radec: Full state including velocities \
                        needed for light travel time and timescale correcion.')  
        
         # Add Light travel time correction to TDB-UTC
        if(lttc):
            # determine state corrected for light travel time
            r0 = norm(state[:,0:3],axis=1)
            # print(dt)
            # print(np.divide(r0,cnst.c_aupd))
            dt = add(dt,divide(r0,cn.CAUPD))
            # print(dt)
            
        state_xyz = array([state[:,0] - multiply(dt,state[:,3]),
                           state[:,1] - multiply(dt,state[:,4]),
                           state[:,2] - multiply(dt,state[:,5])]).T 
        
        r1 = norm(state_xyz[:,0:3],axis=1)                     
        rn = array([divide(state_xyz[:,0],r1),
                    divide(state_xyz[:,1],r1),
                    divide(state_xyz[:,2],r1)]).T
        
        
        #rdot = np.tensordot(rn[:,0:3],state[:,3:6],axes=1)
#         rdot = vec.dot2D(rn,state,
#                           array([0,1,2],dtype='Int32'),array([3,4,5],dtype='Int32'))
        rdot = vec.dot2D(rn,state,
                         array([0,1,2]),array([3,4,5]))
        
        RA = mod(np.arctan2(rn[:,1],rn[:,0])+PIX2, PIX2)
        DEC = arcsin(rn[:,2])
        
        dRAdt = divide((state_xyz[:,0]*state[:,4]-state_xyz[:,1]*state[:,3]),
                        state_xyz[:,0]**2+state_xyz[:,1]**2)
        #dDECdt=divide(divide(state[:,5]),rdot)
        dDECdt = divide(multiply(state[:,5],r1)-multiply(state_xyz[:,2],rdot),
                        multiply(r1,sqrt(r1**2 - state[:,5]**2)))
        
    if(deg):
#         print('RA, DEC')
#         print(np.rad2deg(RA))
#         print(np.rad2deg(DEC))
#         print('dRAdt,dDecdt')
#         print(np.rad2deg(multiply(dRAdt,np.cos(DEC))))
#         print(np.rad2deg(dDECdt))
#         print('r, rdot')
#         print(r1)
#         print(rdot)
        
        radecr = array([rad2deg(RA), rad2deg(DEC),
                           rad2deg(multiply(dRAdt,cos(DEC))),rad2deg(dDECdt),
                           r1, rdot]).T
    else:
        radecr = array([RA, DEC,  multiply(dRAdt,cos(DEC)), dDECdt, r1, rdot]).T
    
    return radecr

def topocentric2ephemeris(epoch, state, frame='icrf', **kwargs):
    """Transform topocentric ICRF states to
    Right Ascension (RA) and Declination (DEC) observations.
    
    Parameters:
    -----------
    epoch             ... epoch of RADEC observation (astropy.time object)
    state             ... topocentric state vector of object [au] 
    frame             ... reference frame ('icrf' or 'ecliptic')
                      ... (given at epoch in TDB)
    Kwargs:
    -------
    epoch_timescale   ... time scale for epoch ('utc', 'tdb'), default 'utc'
    state_timescale   ... time scale for state epoch ('utc', 'tdb'), default 'tdb'
    deg               ... True (default): angles in degrees, False: angles in radians
    lttc              ... True (default): correct for light travel time (needs full state including velocities)
    
    Returns:
    --------
    RA               ... Right Ascension [default deg]
    DEC              ... Declination [default deg]
    dRA/dt*cos(DEC)  ... sky plane motion in RA [default deg/day]
    dDEC/dt          ... sky plane motion in DEC [default deg/day]
    r                ... distance to object [au]
    dr/dt            ... line of sight velocity [au/day]
    """

    if (frame == 'ecliptic'):    
        state_icrf = tr.ecliptic2icrf(state)
        ephemeris = icrf2ephemeris(epoch, state_icrf, **kwargs)
        
    elif(frame == 'icrf'):
        ephemeris = icrf2ephemeris(epoch, state, **kwargs)
         
    else:
        raise Exception('Error in coordinate_transform: \
                         3D positions or 6D state vector required') 
        
    return ephemeris   


def state2ephemeris(epoch, state_asteroid, state_observer, **kwargs):
    """Transform topocentric ICRF states to
    Right Ascension (RA) and Declination (DEC) observations.
    
    Parameters:
    -----------
    epoch             ... epoch of RADEC observation (astropy.time object)
    state_asteroid    ... state vector of asteroid [au, au/day] 
    state_observer    ... state vector of observer [au, au/day]    
    
    Kwargs:
    -------
    frame             ... reference frame for both states ('icrf' or 'ecliptic'), default 'icrf'
    epoch_timescale   ... time scale for epoch ('utc', 'tdb'), default 'utc'
    state_timescale   ... time scale for state epoch ('utc', 'tdb'), default 'tdb'
    deg               ... True: angles in degrees, False: angles in radians
    lttc              ... True: correct for light travel time (needs full state including velocities)
    
    Returns:
    --------
    RA               ... Right Ascension [default deg]
    DEC              ... Declination [default deg]
    dRA/dt*cos(DEC)  ... sky plane motion in RA [default deg/day]
    dDEC/dt          ... sky plane motion in DEC [default deg/day]
    r                ... distance to object [au]
    dr/dt            ... line of sight velocity [au/day]
    """

    # observer to asteroid vectors
    topocentric_state = state_asteroid - state_observer
    # calculate ephemeris
    ephemeris = topocentric2ephemeris(epoch, topocentric_state, **kwargs)
        
    return ephemeris 

def radecResiduals(df, epoch, state_asteroid, output_units='deg', 
                   frame_state_ast='ecliptic', frame_state_obs='ecliptic', propagator='2body',
                   **kwargs):
    """Calculate O-C values in Right Ascension (RA) and Declination (DEC) for a given asteroid state. 
    The state is propagated to all observation epochs and the corresponding RA and DEC values are compared
    to the corresponding observations. Heliocentric ecliptic states for the observer are required for
    every obserbation epoch. 
    
    Parameters:
    -----------
    df                ... Pandas DataFrame containing nightly RA and DEC [deg], time [JD, MJD] UTC,
                          heliocentric ecliptic observer positions and velocities [au]
                
    epoch             ... epoch of asteroid state [JD, MJD], timescale: TDB
    state_asteroid    ... heliocentric ecliptic state of asteroid at epoch (x,y,z,vx,vy,vz), [au, au/day]
    
    Keyword arguments:
    ------------------
    output_units      ... units for O-C results: 'deg', 'arcsec', 'rad'
    frame_state_ast   ... coordinate frame of asteroid state ('icrf' or 'ecliptic')
    frame_state_obs   ... coordinate frame of observer state ('icrf' or 'ecliptic')
    propagator        ... propagator type (e.g. '2body')
    epoch_timescale   ... time scale for observation epoch ('utc', 'tdb'), default utc
    state_timescale   ... time scale for asteroid and observer state epoch ('utc', 'tdb'), default tdb
    deg               ... True (default): angles in degrees, False: angles in radians
    lttc              ... True (default): correct for light travel time \
                          (needs entire state including velocities) 
    propagator        ... propagation algorithm for asteroid state propagation: 
                          '2body' (default), 'linear', 'nbody'                      
    Returns:
    --------
    rms               ... root mean square (RMS) of RA and DEC O-Cs [arcseconds]
    dra               ... Right Ascension O-C values for all observations [arcseconds]
    ddec              ... Declination O-C values for all observations [arcseconds]
    """
    obs_times=df['time'].to_numpy().astype('float')
    dt = epoch-obs_times

    if (frame_state_ast == 'icrf'):
        state_ast=tr.icrf2ecliptic(state_asteroid)
    elif(frame_state_ast == 'ecliptic'):
        state_ast=state_asteroid
    else:
        raise Exception('Error in radecResiduals: frame for asteroid state unkown.')    
        
    ephemeris = []
    ephemeris_app = ephemeris.append
    nobs = len(dt)
    for i in range(nobs):
        state_observer = df[['x_obs', 'y_obs', 'z_obs','vx_obs', 'vy_obs', 'vz_obs']].to_numpy().astype('float')[i]
    
        if (frame_state_obs == 'icrf'):
            state_obs=tr.icrf2ecliptic(state_observer)
        elif(frame_state_obs == 'ecliptic'):
            state_obs=state_observer
        else:
            raise Exception('Error in radecResiduals: frame for observer state unkown.')
    
    # propagate orbit to all observation time
        pstate = pr.propagateState(state_ast[0:3],state_ast[3:6], 
                                   epoch, obs_times[i], propagator=propagator)
        
        state_asteroid_prop = np.array(pstate[0:2]).flatten()
        
        ephemeris_app(state2ephemeris(obs_times[i], state_asteroid_prop, state_obs, 
                                      frame='ecliptic', lttc=True, timescale_state='tdb',
                                      timescale_epoch='utc')) 
    
    # O-C
    ephemeris_array = np.array(ephemeris)
    
#     print('observed', np.array([df['RA'].values, df['DEC'].values]).T)
#     print('calculated', np.array(ephemeris)[:,0:2])
    
    dra = df['RA'].to_numpy().astype('float') - ephemeris_array[:,0]
    ddec =  df['DEC'].to_numpy().astype('float') - ephemeris_array[:,1]
    
    dradec = np.array([dra, ddec]).flatten()
    
    rms = np.sqrt(np.dot(dradec.T, dradec)/nobs)

    if(output_units == 'deg'):
        return rms, dra, ddec
    elif(output_units == 'rad'):
        return np.deg2rad(rms), np.deg2rad(dra), np.deg2rad(ddec)
    elif(output_units == 'arcsec'):
        return rms*3600, dra*3600, ddec*3600
    else:
        raise Exception('Error in radecResiduals: unknown output unit.')