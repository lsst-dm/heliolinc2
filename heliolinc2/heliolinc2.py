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
lsstssp

LSST Solar System Processing routines for
linking of observations of moving objects.

Formerly known as MOPS. Algorithm: Based on HelioLinC (Holman et al. 2018)
we transform topocentric observations to heliocentric states assuming a
distance and radial velocity. The resulting 3D positions are collected into
tracklets. Tracklets contain at least two observations and can,
thus, be used to create velocity vectors. A tracklet + velocity vector is
called an "arrow". Arrows are propagated to a common epoch using spiceypy's
2body propagator, and then clustered using dbscan.

Implementation: S. Eggl 20191215
"""

# Temporary modules:
#from thor.orbits import iod

# Accelerators
import numpy as np
import numba

# Database
import pandas as pd

# Orbital Dynamics
import spiceypy as sp

# Clustering
import scipy.spatial as sc
import sklearn.cluster as cluster

# Trimmed mean
from scipy import stats


# Time transforms
from astropy.time import Time

# Accelerated vector operations
from . import vector as vc

# Constants such as the speed of light and GM
from . import constants as cn

# Coordinate transfroms
from . import transforms as tr

# Orbital Dynamics 
from . import propagate as pr

# Ephemeris generator
from . import ephemeris as ephem

# Utility functions
from . import utility as ut


__all__ = [ 'heliolinc2','obs2heliolinc',
           'selectTrackletsFromObsData', 'cullSameTimePairs',
           'makeHeliocentricArrows', 'observationsInCluster', 
           'meanArrowStatesInClusters','observationsInArrows',
           'collapseClusterSubsets', 'deduplicateClusters']


############################################
# MODULE SPECIFIC EXCEPTION
###########################################

class Error(Exception):
    """HelioLinC3D module specific exception."""
    pass


############################################
# MAIN 
###########################################

def heliolinc2(r, drdt, cr_obs, cr_arrows, ct_min, ct_max, dfobs=pd.DataFrame(), dftrails=pd.DataFrame(),
               observer_interpolant=[],observer_location='I11', origin='SSB', ephemeris_dt='1h',
               trail_dt=0.01, v_max=1, clustering_dimensions=6, light_time=True, 
               verbose=False, min_samples=3, n_jobs=1, mean_state_variance_limit=1):

    """HelioLinC2 (Heliocentric Linking in Cartesian Coordinates) algorithm.

    Parameters:
    -----------        
    r         ... assumed heliocentric distance [au]
    drdt      ... dr/dt assumed heliocentric radial velocity [au/day]
    cr_obs    ... clustering radius for observations mapped to heliocentric positions [au]
    cr_arrows ... clustering radius for propagated arrows [au]
    ct_min    ... minimum timespan between observations to allow for trackelt making [days]
    ct_max    ... maximum timespan between observations to allow for tracklet making [days]
    dfobs     ... Pandas DataFrame containing observation ID (obsId),
                  time [MJD], night, RA [deg], DEC [deg] 
    dftrails  ... Pandas DataFrame containing trail ID (obsId),
                  time, night, RA_center [deg], DEC_center [deg], RA_vel [deg/day], DEC_vel [deg/day] 


    
    Keyword Arguments:
    ------------------
    observer_interpolant     ... scipy.interpolate object, time interpolant for observer ecliptic xyz positions 
                                 (and their derivatives) wrt origin. If not provided, 
                                 heliolinc will generate one through online query to JPL Horizons.
    observer_location        ... Minor Planet Center observatory code for obsever location (e.g. I11)
    origin                   ... Coordinate origin for propagation, e.g. 'Sun' or 'SSB' : Solar System Barycenter
    ephemeris_dt             ... observer ephemeris query timestep
    trail_dt                 ... [days] time step for transforming RADEC +- RADEC_vel*dt into trail RADEC start and endpoint   
    v_max                    ... [au/day] maximum allowable arrow velocity
    clustering_dimensions    ... dimensions for clustering: 3 for position space, 6 for phase space
    light_time               ... Light travel time correction
    verbose                  ... Print progress statements
    min_samples              ... minimum number of samples for clustering with dbscan
    n_jobs                   ... number of processors used for 2body propagation, dbscan and KDTree query
    mean_state_variance_limit... mean_state filtering for variance (e.g. 1e-7) 
                                 increases purity but decreases completeness 

    Returns:
    --------
    obs_in_cluster_df    ... Pandas DataFrame containing linked observation clusters (no prereduction), 
                             r, dr/dt, mean cluster state (ICRF)
    """

    xar = []
    var = []
    tar = []
    obsids_night = []
    trailids_night = []
    
    xar_add = xar.append
    var_add = var.append
    tar_add = tar.append
    obsids_night_add = obsids_night.append
    trailids_night_add = trailids_night.append
    
    vnorm=vc.norm
    array = np.array                  
    unique = np.unique
    concat = np.concatenate 
    
    notrails = False
    noobs = False
    
    nobsarrows=0
    
    if(dftrails.dropna().empty):
        observation_epochs = unique(dfobs['time'].values)
        notrails=True
    elif(dfobs.dropna().empty):
        observation_epochs = unique(dftrails['time'].values)
        noobs=True
    else:    
        observation_epochs = unique(concat([dfobs['time'].values,dftrails['time'].values]))
    
    tmin = min(observation_epochs)
    tmax = max(observation_epochs)
    
    #tmin = min(dfobs['time'].min(),dftrails['time'].min())
    #tmax = max(dfobs['time'].max(),dftrails['time'].max())
                             
    # generate interpolating functions for observer heliocentric positions and velocities as function of time [MJD]  
    
    if(observer_interpolant):
            observerInterpolant = observer_interpolant
    else:
            observerInterpolant = ephem.getObserverInterpolant(tmin,tmax, origin='SSB',
                                                       observer_location=observer_location,
                                                       ephemeris_dt=ephemeris_dt,frame='ecliptic')

    # the following two arrays are for testing purposes only

    # objid_night = []
    # tobs_night = []
    
    # objid_night_add =  objid_night.append
    # tobs_night_add = tobs_night.append
            
    # reference time for propagation (center of linking window)
    tref=(tmin+tmax)*0.5 
    
    if (not noobs):
        obs_nights = unique(dfobs['night'].values)
        
        for n in obs_nights:      
        # GENERATE ARROWS FROM OBSERVATIONS FOR THIS NIGHT
            dfo = dfobs.loc[dfobs['night']==n].reset_index(drop=True)
            if (not dfo.dropna().empty):           
                [xarrow_night, 
                 varrow_night, 
                 tarrow_night, 
                 goodpairs_night]=makeHeliocentricArrows(dfo,r,
                                                     drdt,tref,cr_obs,ct_min, ct_max,
                                                     observerInterpolant,
                                                     v_max=v_max,lttc=light_time, eps=cr_obs,
                                                     filtering=True, verbose=False, 
                                                     leafsize=16, balanced_tree=False, 
                                                     n_jobs=n_jobs)

            # ADD TO PREVIOUS ARROWS
            if (len(xarrow_night)<1):
                if (verbose):
                    print('no arrows from observations in night ',n)
            else:
                xar_add(xarrow_night)
                var_add(varrow_night)
                tar_add(tarrow_night)
                obsids_night_add(dfo['obsId'].values[goodpairs_night])

        #objid_night_add(df['objId'].values[goodpairs_night])
        #tobs_night_add(df['time'].values[goodpairs_night])

            nobsarrows=len(tar)
                
    if (not notrails):
        trail_nights = unique(dftrails['night'].values)
        print('trail_nights',trail_nights)
        
        for n in trail_nights:                    
        # GENERATE HELIOCENTRIC ARROWS FROM TOPOCENTRIC TRAILS 
            dft = dftrails.loc[dftrails['night']==n].reset_index(drop=True)
            if (not dft.dropna().empty): 
                [xarrow_night, 
                 varrow_night, 
                 tarrow_night, 
                 trails_night]= makeArrowsFromTrails(dft, r, 
                                                 drdt, tref, observerInterpolant, 
                                                 dt = trail_dt, v_max=v_max,
                                                 GM=cn.GM, lttc=False)

            # ADD TO PREVIOUS ARROWS
            if (len(xarrow_night)<1):
                if (verbose):
                    print('no arrows from observations in night ',n)
            else:
                xar_add(xarrow_night)
                var_add(varrow_night)
                tar_add(tarrow_night)
                # trail IDs are negative
#                 print('trails',len(trails_night))
#                 print('arrows',len(tarrow_night))
#                 print(trails_night[0],trails_night[-1])
#                 print(n,dft['trailId'].values[trails_night])
                trailids_night_add(dft['trailId'].values[trails_night])                             
            ntrailarrows=len(tar)-nobsarrows
#   Collect all arrows for propagation                                   
    if (len(xar)<1):
        if (verbose):
            print('No arrows for the current r, dr/dt pair. ',r,drdt)
    else:    
        xarrow=np.vstack(xar)
        varrow=np.vstack(var)
        tarrow=np.hstack(tar)
        
#       Which observations are in each arrow? -> obsids   
        print('obsids_night',obsids_night)
        obsids=[]
        trailids=[]
        if(not noobs):
            obsids = np.vstack(obsids_night)
            print('obsids', obsids)
            print('nobsarrows',nobsarrows,len(obsids))
        if(not notrails):    
            trailids = np.hstack(trailids_night)
            print('trailids',trailids)
            print('ntrailarrows',ntrailarrows,len(trailids))
    
        if(verbose):
            print('tarrow', tarrow)
            print('xarrow', xarrow)
            print('varrow', varrow) 

#   # the following two arrays are for testing purposes only
        #objids=np.vstack(objid_night)
        #tobs=np.vstack(tobs_night)

#   # PROPAGATE ARROWS TO COMMON EPOCH
        if (verbose):
            print('Propagating arrows...')
            
#   # propagate arrows to the center of the observational arc    
        tprop=tref
        [xp, vp, dt] = pr.propagateState(xarrow, varrow, tarrow, tprop,
                                         propagator='2body', n_jobs=n_jobs)

        xpvp=np.hstack([xp,vp])
#       # nomalize v vector to have approximately the same magnitude as xp to help with 6D clustering
        if (clustering_dimensions==6):
            vnormr=(r/vnorm(vp))
            vpn=vp*np.array([vnormr,vnormr,vnormr]).T
            xpvpn=np.hstack([xp,vpn])
        
        if (verbose):
            print('Propagated arrows xyz', xpvp[:,0:3])
            print('Propagated arrows vxyz', xpvp[:,3:6])
        
#       # CLUSTER WITH DBSCAN
        if (verbose):
            print('Clustering arrows...')
            
#       # CLUSTER PROPAGATED STATES (COORDINATE SPACE (xp) OR PHASE SPACE (xpvp)               
#         if(clustering_algorithm=='dbscan'):
        if (clustering_dimensions==6):
            db=cluster.DBSCAN(eps=cr_arrows, min_samples=min_samples, n_jobs=n_jobs).fit(xpvpn)
        elif (clustering_dimensions==3):
            db=cluster.DBSCAN(eps=cr_arrows, min_samples=min_samples, n_jobs=n_jobs).fit(xp)
        else:
            raise Exception('Error: only clustering in 3 or 6 dimensions supported.')
            
        print(db.labels_)    
#       # CONVERT CLUSTER INDICES TO OBSERVATION INDICES IN EACH CLUSTER
        try:
            if (verbose):
                print('Calculating mean arrow states in clusters...')
            [mean_states, var_states, labels2] = meanArrowStatesInClusters(xpvp, db, garbage=False, trim=25) 

            if (verbose):
                print('Finding observations in clusters...')
#             if(notrails):
#                 [obs_in_cluster, labels] = observationsInCluster(obsids, db, garbage=False)
#                  obs_in_cluster_df=pd.DataFrame(zip(labels,obs_in_cluster), columns=['clusterId','obsId'])
#             elif:
            
            [obs_in_cluster, trails_in_cluster, labels] = observationsInCluster(obsids,trailids, db, garbage=False)
            obs_in_cluster_df = pd.DataFrame(zip(labels,obs_in_cluster,trails_in_cluster), 
                                             columns=['clusterId','obsId','trailId'])
#             [obs_in_cluster, labels] = observationsInCluster(obsids, db, garbage=False)
#             obs_in_cluster_df=pd.DataFrame(zip(labels,obs_in_cluster), columns=['clusterId','obsId'])

        except: 
            raise Exception('Error: Could not construct cluster dataframe.')

    
        # Add heliocentric r, dr/dt, epoch and clipped mean states (ICRF) to pandas DataFrame
        obs_in_cluster_df['r'] = r
        obs_in_cluster_df['drdt'] = drdt
        obs_in_cluster_df['cluster_epoch'] = tprop
        
        if (len(mean_states) < 1):
            obs_in_cluster_df[['x_ecl','y_ecl','z_ecl','vx_ecl',
                               'vy_ecl','vz_ecl','var_pos','var_vel']] = pd.DataFrame([np.zeros(8)],                                                                                                            index=obs_in_cluster_df.index)
        else:
            obs_in_cluster_df['x_ecl'] = mean_states[:,0]
            obs_in_cluster_df['y_ecl'] = mean_states[:,1]
            obs_in_cluster_df['z_ecl'] = mean_states[:,2]
            obs_in_cluster_df['vx_ecl'] = mean_states[:,3]
            obs_in_cluster_df['vy_ecl'] = mean_states[:,4]
            obs_in_cluster_df['vz_ecl'] = mean_states[:,5]
            
            obs_in_cluster_df['var_pos'] = vnorm(var_states[:,0:3])
            obs_in_cluster_df['var_vel'] = vnorm(var_states[:,3:6])
        
        if (mean_state_variance_limit>0):
            if (verbose):
                print('Filtering clusters for mean state variance...')
            obs_in_cluster_df = obs_in_cluster_df[
                                (obs_in_cluster_df['var_pos']<=mean_state_variance_limit)].reset_index(drop=True)  

        return obs_in_cluster_df


def observationsInCluster(pairs, trails, cluster, garbage=False):
    """List observations in each cluster of arrows.
    
    Parameters:
    -----------
    df      ... pandas dataframe with observations
    pairs   ... list of pairs of observations [obsid1,obsid2] linked into arrows
    cluster ... output of clustering algorithm (sklearn.cluster)
    
    Returns:
    --------
    obs_in_cluster ... list of observations in each cluster
    """
    #cluster names (beware: -1 is the cluster of all the leftovers)
    if(garbage):
        unique_arrow_labels = np.unique(cluster.labels_)
    else:
        unique_arrow_labels = np.unique(cluster.labels_)[0:]
    #number of clusters
    n_clusters = len(unique_arrow_labels)
    
    #which objects do observations in pairs (tracklets) belong to
    p = np.array(pairs)
    lenp = len(pairs)
    print('pairs',pairs)
    print('lenp',lenp)
    print('n_clusters',n_clusters)
    print('unique_arrow_labels',unique_arrow_labels)
    obs_in_cluster=[]
    obs_in_cluster_add=obs_in_cluster.append
    
    trails_in_cluster=[]
    trails_in_cluster_add=trails_in_cluster.append
    
    # cluster contains 
    for u in unique_arrow_labels:
        # which indices of arrows appear in a given cluster?
        idx = np.where(cluster.labels_ == u)[0]
#         print('idx',idx)
#         print('p[idx]',np.unique(p[idx].flatten()))
        # which observations are in this cluster
        obs_idx=np.where(idx < lenp)[0]
#         print('obs_idx',obs_idx)
#         print('p[obs_idx]',np.unique(p[idx[obs_idx]].flatten()))
        if obs_idx.size == 0:
            obs_in_cluster_add(-1)
        else:    
            obs_in_cluster_add(np.unique(p[idx[obs_idx]].flatten()))
        
        #which trails are in this cluster
        trail_idx = np.where(idx >= lenp)[0]
#         print('idx',idx)
#         print('trail_idx',trail_idx)
#         print('trails',trails)
#         print('trails_add',np.unique(trails[idx[trail_idx]-lenp]))
        if trail_idx.size == 0:
            trails_in_cluster_add(-1)
        else:
            trails_in_cluster_add(np.unique(trails[idx[trail_idx]-lenp]))
        
    return obs_in_cluster, trails_in_cluster, unique_arrow_labels      
    
#####
    
# def observationsInCluster(pairs, cluster, garbage=False):
#     """List observations in each cluster.
    
#     Parameters:
#     -----------
#     pairs   ... list of pairs of observations [obsid1,obsid2] linked into arrows
#     cluster ... output of clustering algorithm (sklearn.cluster)
    
#     Returns:
#     --------
#     obs_in_cluster ... list of observations in each cluster
#     """
#     #cluster names (beware: -1 is the cluster of all the leftovers)
#     if(garbage):
#         unique_labels = np.unique(cluster.labels_)
#     else:
#         unique_labels = np.unique(cluster.labels_)[1:]
#     #number of clusters
#     n_clusters = len(unique_labels)
    
#     #which objects do observations in pairs (tracklets) belong to
#     p = np.array(pairs)              

#     obs_in_cluster=[]
#     obs_in_cluster_add=obs_in_cluster.append
    
#     #cluster contains 
#     for u in unique_labels:
#         #which indices in pair array appear in a given cluster?
#         idx = np.where(cluster.labels_ == u)[0]      
# #         print('idx',idx)
# #         print('p[idx]',np.unique(p[idx].flatten()))
#         #which observations are in this cluster
#         obs_in_cluster_add(np.unique(p[idx].flatten()))
        
#     return obs_in_cluster, unique_labels     
    
############################################
# OBSERVATIONS, TRACKLETS AND ARROWS
###########################################
    
def obs2heliolinc(df, required_input_columns={'obsName':'obsName','time':'FieldMJD',
                                              'RA':'AstRA(deg)','DEC':'AstDec(deg)'},
                 required_output_columns={'obsName':'obsName','time':'time',
                                          'RA':'RA','DEC':'DEC','obsId':'obsId'}):
    
    """Converts observations in pandas DataFrame to HelioLinC3D ingestable format.
    
    Parameters:
    ------------
    df                     ... pandas DataFrame containing observations with time in UTC MJD 
    required_input_columns ... dictionary, input columns that must be present in input pandas DataFrame
    
                               obsName: str, name of observation
                               FieldMJD: float, UTC time of observation [MJD]
                               AstRA(deg): astrometric Right Ascension [deg] 
                               AstDec(deg): astrometric Declination [deg]                             
                              
    required_output_columns ... dictionary, output columns necessary for HelioLinC3D:
    
                                obsId: int, unique observation index
                                obsName: str, name of observation
                                time: float, TDB time of observation [MJD]
                                RA: astrometric Right Ascension [deg] 
                                DEC: astrometric Declination [deg]
                           
    Returns:
    --------
    df2                ... pandas DataFrame containing required output columns as well as time in TDB MJD
    
    
    Remarks:
    --------
    Input timescale is UTC and format is Modified Julian Date (MJD).
    Input astrometry is not corrected for aberration and light travel time. 
    This will be done later in the program.
    
    """

                      
#     required_output_columns=['obsId','obsName','time','RA','DEC',
#                       'x_obs','y_obs','z_obs',
#                      'vx_obs','vy_obs','vz_obs']

    for col in required_input_columns.values:
        if col not in df2.columns:
            raise Exception("Error in obs2heliolinc: required input column not present:",col)
                
    
    icols = [required_input_columns[r] for r in required_input_columns]
    ocols = [required_output_columns[r] for r in required_output_columns]
    
    # Create a copy with only required inputs
    df2 = df[icols]
    
    # Rename columns of observations DataFrame to work with HelioLinC3D
    column_mapping = dict(zip(icols,icols2))  
    df2.rename(columns = column_mapping, inplace=True) 
    
    # Transform time from UTC to TDB
    df2['time'] = Time(df2['time'].values,scale='utc',format='mjd').tdb.mjd
  
  
    # Sort by observation time
    # df2.sort_values(by=['time','RA'], inplace=True)
    
    df2.reset_index(inplace=True, drop=False)
    
    # Create observation ID         
    df2['obsId'] = df2.index
    
    # print(df2)
    # If a unique observation Identifier is present in the dataframe use that
    if (uniqueObsId):
        df2['obsName'] = df[uniqueObsIdName].values.astype(str)
    # Otherwise create your own
    else:    
        df2['obsName'] = df2['obsId'].values.astype(str)
      
    df2['night']=(ut.lsstNight(df2['time'],df2['time'].min())).astype(int)
    
    # Check if all required columns are present
    for col in required_output_columns:
           if col not in df2.columns:
                raise Exception("Error in obs2heliolinc: required columns not present.")
               
    
    return df2, observerInterpolant


def cullSameTimePairs(pairs, df, dt_min, dt_max, time_column_name):
    """Cull all observation pairs that occur at the same time.

    Parameters:
    -----------
    pairs             ... array of observation indices that form a tracklet
                          (only 2 entries are supported for now)
    df                ... pandas dataframe containing observation data
    dt_min            ... minimum timespan between two observations 
                          to be considered a pair, e.g. exposure duration (days)
    dt_max            ... maximum timespan between two observations 
                          to be considered a pair (days)
    time_column_name  ... string, the name of the pandas dataframe column
                          containing the epoch of each observation

    Returns:
    --------
    goodpairs         ... array of observation indices that form possible
                          tracklets (only 2 observation entries
                          per tracklet are supported for now)
    """

    tn = time_column_name
    # Tracklets from observation paris cannot be constructed from contemporal observations (delta_t==0)
    # nor observations that are too far apart (>= dt_max)
    delta_t = np.abs(df[tn].values[pairs[:,1]]-df[tn].values[pairs[:,0]])
    goodpairs = pairs[(delta_t>dt_min) & (delta_t<dt_max)]
    return np.array(goodpairs)
    
    
def selectTrackletsFromObsData(pairs, df, dt_min, dt_max, time_column_name):
    """Select data in trackelts from observations data frame.

    Parameters:
    -----------
    pairs             ... array of observation indices that form a tracklet
                          (only 2 entries are supported for now)
    df                ... pandas dataframe containing observation data
    dt_min            ... minimum timespan between two observations 
                          to be considered a pair, e.g. exposure duration (days)
    dt_max            ... maximum timespan between two observations 
                          to be considered a pair (days)
    time_column_name  ... string, the name of the pandas dataframe column
                          containing the epoch of each observation

    Returns:
    --------
    df2               ... pandas dataframe containing only observations
                          that occur in tracklets.
    goodpairs         ... array of observation indices that form possible
                          tracklets (only 2 observation entries per
                          tracklet are supported for now)
    """
    goodpairs = cullSameTimePairs(pairs, df, dt_min, dt_max, time_column_name)
    index_list = np.unique(goodpairs.flatten())
    #df2 = (df.iloc[index_list]).reset_index()

    return goodpairs

def makeHeliocentricArrows(df, r, drdt, tref, cr, ct_min, ct_max, observerInterpolant, 
                           v_max=1., eps=0,
                           lttc=False, filtering=True, verbose=False, 
                           leafsize=16, balanced_tree=True, n_jobs=1, 
                           GM=cn.GM):

    """Create tracklets/arrows from dataframe containing nightly RADEC observations
    and observer positions.

    Parameters:
    -----------
    df       ... Pandas DataFrame containing nightly RA and DEC [deg], time [MJD]
    r        ... assumed radius of heliocentric sphere used for arrow creation [au]
    drdt     ... assumed radial velocity [au/day]
    tref     ... reference time for arrow generation. Used to calculate how much the 
                 heliocentric distance changes between observations based on assumed dr/dt
    cr       ... maximum spacial clustering radius for arrow creation (au)
    ct_min   ... minimum temporal clusting radius for arrow creation (days)
    ct_max   ... maximum temporal clusting radius for arrow creation (days)
    observerInterpolant ... scipy interpolant for heliocentric observer position [au, ecliptic]


    Keyword arguments:
    ------------------
    v_max (optional)       ... velocity cutoff [au/day]
    lttc (optional)        ... light travel time correction
    filtering (optional)   ... filter created tracklets (exclude tracklets built 
                               from data with the same timestamp) 
    verbose (optional)     ... print verbose progress statements  
    eps (optional)         ... Branches of the Kdtree are not explored if their 
                               nearest points are further than r/(1+eps), 
                               and branches are added in bulk if their furthest points 
                               are nearer than r * (1+eps). eps has to be non-negative.
    leafsize               ... cKDTree leafsize 
    balanced_tree          ... cKDTree: if True median is used 
                               instead of mean to calculate box centers, more robust but 
                               leads to longer tree build times
    n_jobs                 ... number of available processors for simultaneous cKDTree query
    GM                     ... gravitational constant * mass for central body (e.g. gaussk**2*M_sun)                    


    Returns:
    --------
    x         ... tracklet/arrow position (3D) [au]
    y         ... tracklet/arrow velocity (3D) [au]
    t         ... tracklet/arrow reference epoch [JD/MJD]
    goodpairs ... index pairs of observations that go into each tracklet/arrow
    """
    
    sqrt = np.sqrt
    array = np.array
    divide = np.divide
    where = np.where
    take = np.take
    
    goodpairs=[]
    paris=[]

    # Transform RADEC observations into positions on the unit sphere (ICRF)
    xyz = tr.radec2icrfu(df['RA'], df['DEC'], deg=True)

    # Those are the line of sight (LOS) vectors
    los_icrf = array([xyz[0], xyz[1], xyz[2]]).T

    # Transform LOS vectors to frame of choice
    los = tr.frameChange(los_icrf, 'icrf', 'ecliptic') 
    
    # Use the position of the observer and the LOS to project the position of
    # the asteroid onto a heliocentric great circle with radius r
    observerXYZ = observerInterpolant(df['time'].values).T

    # Calculate how much the heliocentric distance changes
    # during the obsevations based on assumed dr/dt
    # we assume we have r at t=tref. if dr/dt>0  and t<tref we have
    # to move to smaller r
    
    dt = df['time'].values-tref
    
    # calculate angular momentum h = r x v, 
    # split velocity in parallel and orthogonal components to r: v = v_o + v_r 
    # then h= r x v_o + r x v_r 
    # since r x v_r = 0 (parallel) h = r x v_o 
    # v_o = v - v_r = v - rdot, and |h| = r*(v-rdot)
    # for circular orbits rdot = 0 and v = sqrt(GM/r)  
    # with F(r) = m d^2r/dt^2 - m r w^2 = m d^2r/dt^2 - m h^2/r^3
    #     dt = tref-df['time'].values
        
    h = r*(sqrt(GM/r)-drdt)    
    dr = drdt*dt + (-GM/(r*r)+h**2/r**3)*dt**2/2 
    r_plus_dr = r+dr
    print('tref,t0',tref,df['time'].values)
    print('r,r+dr',r,r_plus_dr)
    
#     r_plus_dr,rdot_new = vleapfrog(dt,r,drdt,0.5,GM,h)
    
    
    # Heliocentric postions of the observed asteroids
    posu = vc.sphereLineIntercept(los, observerXYZ, r_plus_dr)

    if(verbose):
        print('Heliocentric positions generated.')
        print('Building spacial KDTree...')
        
    # To generate tracklets we build our KDTree based on the positions
    # in heliocentric space
    kdtree_s = sc.cKDTree(posu, leafsize=leafsize, compact_nodes=True,
                          copy_data=False, balanced_tree=balanced_tree, 
                          boxsize=None)
    # rule out un-physical combinations of observations with kdtree

    # Query KDTree for good pairs of observations that lie within
    # the clustering radius cr
    if(verbose):
        print('KDTree generated. Creating tracklets...')
        
    pairs = kdtree_s.query_pairs(cr, p=2., eps=eps, 
                                 output_type='ndarray')

    if(verbose):
        print('Tracklet candidates found:',len(pairs))

    if (filtering):
        if(verbose):
            print('Filtering arrows by time between observations...')
        
        # Discard impossible pairs (same timestamp)
        goodpairs = selectTrackletsFromObsData(pairs, df, ct_min, ct_max, 'time')
        
        if(verbose):
            print('Tracklets filtered. New number of tracklets:',len(goodpairs))
    
    else:
        goodpairs=pairs
    
    # tracklet position for filtered pairs
    x = posu[goodpairs[:,0]]
    # tracklet time
    t = df['time'].values[goodpairs[:,0]]
    # tracklet velocity through forward differencing
    va = []
    vapp = va.append
    dt = df['time'].values[goodpairs[:,1]]-df['time'].values[goodpairs[:,0]]
    dx = posu[goodpairs[:,1]]-posu[goodpairs[:,0]]
    for d in range(0,3):
        vapp(divide(dx[:,d],dt))
    v = array(va).T
    
    if (filtering):
        if(verbose):
            print('Filtering arrows by max velocity...')
        vnorm = vc.norm(v)
        v_idx = where(vnorm<=v_max)[0]
    
        goodpairs=np.take(goodpairs,v_idx,axis=0)
        x = take(x,v_idx,axis=0)
        v = take(v,v_idx,axis=0)
        t = take(t,v_idx,axis=0)
    
    if(verbose):
        print('Tracklets created:',len(goodpairs))
    
    # correct arrows for down-leg light travel time. At the time of observation, the source has
    # moved from its previous position where it emitted light to a new position advancing on its orbit
    if(lttc):
        if(verbose):
            print('(Linear correction for light travel time aberration...')
        xo = observerInterpolant(t).T
        dist = vc.norm(x-xo)
        xl = x.T+dist/cn.CAUPD*v.T
        x = xl.T 
 
    return x, v, t, goodpairs


def makeArrowsFromTrails(df, r, drdt, tref, observerInterpolant, dt = 0.01, v_max=1.,
                         GM=cn.GM, verbose=False, filtering=True, lttc=False):

    """Create tracklets/arrows from dataframe containing nightly RADEC observations
    and observer positions.

    Parameters:
    -----------
    df       ... Pandas DataFrame containing trail data RA_center, DEC_center [deg], RA_vel, DEC_vel [deg/day]
                 time [JD, MJD]
    r        ... assumed radius of heliocentric sphere used for arrow creation [au]
    drdt     ... assumed radial velocity [au/day]
    tref     ... reference time for arrow generation. Used to calculate how much the 
                 heliocentric distance changes between observations based on assumed dr/dt
    observerInterpolant ... scipy interpolant for heliocentric observer position 
    dt       ... [days] time step for transforming RADEC +- RADEC_vel*dt into trail start and endpoint
                 dt should be based on exposure time. 


    Keyword arguments:
    ------------------
    v_max (optional)       ... velocity cutoff [au/day]
    frame                  ... frame of reference for heliocentric arrows: 'icrf' or 'ecliptic'
    GM                     ... gravitational constant * mass for central body (e.g. gaussk**2*M_sun)                    


    Returns:
    --------
    x         ... tracklet/arrow position (3D) [au]
    y         ... tracklet/arrow velocity (3D) [au]
    t         ... tracklet/arrow reference epoch [JD/MJD]
    goodpairs ... index pairs of observations that go into each tracklet/arrow
    """

    sqrt = np.sqrt
    trailidx = df.reset_index().index.values
    
    df['RA0'] = df['RA_center']-df['RA_vel']*dt
    df['RA1'] = df['RA_center']+df['RA_vel']*dt
    df['DEC0']= df['DEC_center']- df['DEC_vel']*dt
    df['DEC1']= df['DEC_center']+ df['DEC_vel']*dt
    
    # Transform RADEC observations into positions on the unit sphere (ICRF)
    xyz0 = tr.radec2icrfu(df['RA0'], df['DEC0'], deg=True)
    xyz1 = tr.radec2icrfu(df['RA1'], df['DEC1'], deg=True)
    
    # Those are the line of sight (LOS) vectors
    los_icrf0 = np.array([xyz0[0], xyz0[1], xyz0[2]]).T
    los_icrf1 = np.array([xyz1[0], xyz1[1], xyz1[2]]).T

    # Transform LOS vectors to frame of choice
    los0 = tr.frameChange(los_icrf0, 'icrf', 'ecliptic') 
    los1 = tr.frameChange(los_icrf1, 'icrf', 'ecliptic') 
    
    # Use the position of the observer and the LOS to project the position of
    # the asteroid onto a heliocentric great circle with radius r
    # fix that observer1 = df[['x_obs', 'y_obs', 'z_obs']].values
    observer0 = observerInterpolant(df['time'].values-dt).T
    observer1 = observerInterpolant(df['time'].values+dt).T

    # Calculate how much the heliocentric distance changes
    # during the obsevations based on assumed dr/dt    
    dt0 = df['time'].values-tref-dt
    dt1 = df['time'].values-tref+dt
    
     # calculate angular momentum h = r x v, 
    # split velocity in parallel and orthogonal components to r: v = v_o + v_r 
    # then h= r x v_o + r x v_r 
    # since r x v_r = 0 (parallel) h = r x v_o 
    # v_o = v - v_r = v - rdot, and |h| = r*(v-rdot)
    # for circular orbits rdot = 0 and v = sqrt(GM/r)  
    # with F(r) = m d^2r/dt^2 - m r w^2 = m d^2r/dt^2 - m h^2/r^3
    # dt = tref-df['time'].values
    # h = r*(sqrt(GM/r)-drdt)    
    # dr = drdt*dt + (-GM/(r*r)+h**2/r**3)*dt**2/2 
    # r_plus_dr = r+dr
    
    h = r*(sqrt(GM/r)-drdt)  
    dr0 = drdt*dt + (-GM/(r*r)+h**2/r**3)*dt0**2/2
    dr1 = drdt*dt + (-GM/(r*r)+h**2/r**3)*dt1**2/2
    
#     dr0 = (drdt-GM/(r*r)*dt0/2)*dt0
#     dr1 = (drdt-GM/(r*r)*dt1/2)*dt1
    
    r_plus_dr0 = r+dr0
    r_plus_dr1 = r+dr1

    # Heliocentric postions of the observed asteroids
    posu0 = vc.sphereLineIntercept(los0, observer0, r_plus_dr0)
    posu1 = vc.sphereLineIntercept(los1, observer1, r_plus_dr1)

    if(verbose):
        print('Heliocentric positions generated.')
        print('Building trail arrows...')
        
    
    # tracklet position for filtered pairs
    x = 0.5*(posu0+posu1)
    # tracklet time
    t = df['time'].values
    # tracklet velocity through central differencing
    va = []
    vapp = va.append
    dta = 2*dt
    dx = posu1-posu0
    
    for d in range(0,3):
        vapp(np.divide(dx[:,d],dta))
    v = np.array(va).T
    
    if (filtering):
        if(verbose):
            print('Filtering arrows by max velocity...')
        vnorm=vc.norm(v)
        v_idx=np.where(vnorm<=v_max)[0]
    
        trailidx=np.take(trailidx,v_idx,axis=0)
        x=np.take(x,v_idx,axis=0)
        v=np.take(v,v_idx,axis=0)
        t=np.take(t,v_idx,axis=0)
    
    if(verbose):
        print('Tracklets created:',len(trailidx))
    
    # correct arrows for down-leg light travel time. At the time of observation, the source has
    # moved from its previous position where it emitted light to a new position advancing on its orbit
    if(lttc):
        if(verbose):
            print('(Linear correction for light travel time aberration...')
        xo = observerInterpolant(t[trailidx]).T
        dist = vc.norm(x-xo)
        xl = x.T+dist/cn.CAUPD*v.T
        x = xl.T 
 
    return x, v, t, trailidx

    
def leapfrog(tend,r,rdot,dt,gm,h):
    """Leap frog integrator for radial Kepler problem
    d2r/dt2 = -GM/r^2 + h/r^3, where r is the scalar distance 
    to the focus (e.g. Sun), t is the time, GM the gravitational parameter,
    and h the (specific) angular momentum of the orbit.
    
    Parameters:
    -----------
    tend ... final time of propagation
    r    ... radial distance to focus (e.g. Sun)
    rdot ... dr/dt, radial velocity
    dt   ... time step
    GM   ... gravitational parameter 
    h    ... angular momentum

    """
    t=0
    
    if (tend<0):
        dt = -dt
    
    dth = 0.5*dt
    rdoth = rdot+dth*(-gm/r**2+h**2/r**3)
    
    #rtest=[]
    while (abs(t)<abs(tend-dt)):
        r = r+dt*rdoth
        rdoth = rdoth+dt*(-gm/r**2+h**2/r**3)
        t=t+dt
        # rtest.append([t,r,rdoth-dth*(-gm/r**2+h**2/r**3)])
        
    # last step
    r = r+(tend-t)*rdoth
    rdot = rdoth-dth*(-gm/r**2+h**2/r**3)
    return r, rdot #, np.array(rtest)  

vleapfrog = np.vectorize(leapfrog)
    
def deduplicateClusters(cdf):    
    """Deduplicate clusters produced by helioLinC2 
    based on similar observations (r,rdot are discarded)

    Parameters:
    -----------
    cdf ... Pandas DataFrame containing object ID (objId), observation ID (obsId)
              
    Returns:
    --------
    cdf2     ... deduplicated Pandas DataFrame 
    """
    
#   cdf2=(cdf.iloc[cdf.astype(str).drop_duplicates(
#         subset='obsId',keep="first").index]).reset_index(drop=True)    

    cdf['clusterObs'] = cdf['obsId'].astype(str) 
    cdf.sort_values(by=['clusterObs','var_pos','var_vel'],inplace=True)
    cdf.drop_duplicates(subset='clusterObs',keep="first",ignore_index=True, inplace=True)
    return cdf.drop(columns=['clusterObs'])
        
def collapseClusterSubsets(cdf):
    """Merge clusters that are subsets of each other 
    as produced by HelioLinC2.

    Parameters:
    -----------
    cdf ... Pandas DataFrame containing object ID (objId), observation ID (obsId)
              
    Returns:
    --------
    cdf2                ... collapsed Pandas DataFrame 
    subset_clusters     ... indices of input dataframe (cdf) that are subsets
    subset_cluster_ids  ... linked list of cluster ids [i,j] where j is a subset of i
    """
    
   
    #for index, row in clusters_df.iterrows():
    vals=cdf.obsId.values
    subset_clusters=[]
    subset_clusters_app=subset_clusters.append
    subset_cluster_ids=[]
    subset_cluster_ids_app=subset_cluster_ids.append

    cdf_idx=range(0,len(cdf))

    vals_set=[]
    vals_set_app=vals_set.append
    vals_min=[]
    vals_min_app=vals_min.append
    vals_max=[]
    vals_max_app=vals_max.append

    for i in cdf_idx:
        vals_set_app(set(vals[i]))          
        vals_min_app(np.min(vals[i]))
        vals_max_app(np.max(vals[i]))         

    vmin=np.array(vals_min)
    vmax=np.array(vals_max)

    for i in cdf_idx:
        for j in cdf_idx:
            if(i != j):
                    #use the fact that we have ordered sets here
                    if(vmax[i]<vmin[j]):
                        break
                    elif(vmin[i]>vmax[j]):
                        break
                    else:
                        is_subset=vals_set[i].issubset(vals_set[j])

                        if(is_subset):
                            subset_clusters_app(i)
                            subset_cluster_ids_app([i,j])
                            break
        if(np.mod(i, 1000) == 0):                
            print('Progress [%], ', i/cdf_idx[-1]*100)

    idx_to_drop = subset_clusters
    
    cdf2=cdf.drop(index=idx_to_drop)
    return cdf2, subset_clusters, subset_cluster_ids 


def meanArrowStatesInClusters(xpvp, cluster, garbage=False, trim=25):
    """Calculate mean propagated state in each cluster.
    
    Parameters:
    -----------
    xpvp      ... array containing states for each cluster
    cluster   ... output of clustering algorithm (sklearn.cluster)
    
    Keyword Arguments:
    ------------------
    garbage   ... return results for garbage cluster in dbscan (usually index 0)
    trim      ... cutoff percentage for trimmed mean, 
                  25 would slices off ‘leftmost’ and ‘rightmost’ 25% of scores.
    
    Returns:
    --------
    mean_state     ... trimmed mean state in each cluster
    unique_labels  ... cluster labels
    """
    #cluster names (beware: -1 is the cluster of all the leftovers)
    if(garbage):
        unique_labels = np.unique(cluster.labels_)
    else:
        unique_labels = np.unique(cluster.labels_)[0:]
    #number of clusters
    n_clusters = len(unique_labels)
         
    mean_states=[]
    mean_states_add=mean_states.append
    tmean=stats.trim_mean
    
    var_states=[]
    var_states_add=var_states.append
    tvar=stats.tvar
    
    #cluster contains 
    for u in unique_labels:
        idx = np.where(cluster.labels_ == u)[0]
        
        # trimmed mean state in this cluster
        mean_states_add(tmean(xpvp[idx,:], trim/100, axis=0))
        # trimmed variance in this cluster
        var_states_add(tvar(xpvp[idx,:], axis=0, ddof=1))
    return np.array(mean_states), np.array(var_states), unique_labels  


 



def observationsInArrows(df, goodpairs, *args, **kwargs):
    """Find which observations go into an arrow.
    
    Parameters:
    -----------
    df          ... pandas dataframe with observations
    goodpairs   ... filtered list of pairs of observations [obsid1,obsid2] linked into arrows
    
    Returns:
    --------
    df_obs_in_arrows ... pandas dataframe where the index is the arrow id 
                         and the first and second column are the first and second observation 
                         that go into the respective arrow.
    """
    df_obs_in_arrows = pd.DataFrame(goodpairs,**kwargs)
    return df_obs_in_arrows


