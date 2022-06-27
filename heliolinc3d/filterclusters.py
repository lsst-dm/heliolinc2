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
filterclusters

LSST Solar System Processing

Filtering of HelioLinC clusters via orbit determination
Implementation: Python 3.6, S. Eggl 20200310
"""

# Accelerators
import numpy as np

# Pandas DataFrames
import pandas as pd

# HelioLinC2 modules
from . import constants as const
from . import ephemeris as ephem

# Thor Initial Orbit Determination
#from thor.orbits import iod


__all__ = ['filterClusters', 'collapseClusterSubsets']


def filterClusters(df_obs, df_cluster, rms_max=200, filter_type='mean_state'):
    """Filter clusters of observations via Initial Orbit Determination
    
    Parameters:
    -----------
    df_obs     ... Pandas DataFrame containing obsId, epoch an RA DEC values for all observations
    df_cluster ... Pandas DataFrame containing lists of linked observations
    index      ... Pandas index / clusterId
    row        ... Pandas row containing obsIds in the respective cluster with index index
    rms_max    ... RMS O-C cutoff for filter
    
    Keyword Arguments:
    ------------------
    filter_type   ... type of filter: 'mean_state'
    
    Returns:
    --------
    idx      ... list of Pandas index / clusterIds
    drop_idx ... list of True/False: does cluster comply with RMS limit?
    rms      ... lsit of actual RMS values for the clusters
    """
    
    drop_idx=[]
    drop_idx_add=drop_idx.append
    drop_clusterId=[]
    drop_clusterId_add=drop_clusterId.append
    idx=[]
    idx_add=idx.append
    cluster_rms=[]
    cluster_rms_add=cluster_rms.append
    
    for index, row in df_cluster.iterrows():
        idx_add(index)
        df=df_obs[[x in row['obsId'] for x in df_obs['obsId']]]
           
        rms = -999
        try:
            [rms, epoch, orbit] = meanStateFilter(df, index, row, 
                                                  frame_state_asteroid='icrf',
                                                  frame_state_observer='icrf')              
            if(rms<rms_max):
                drop_idx_add(False) 
                
            else:
                drop_idx_add(True)
                drop_clusterId_add(row['clusterId'])
        except:
            drop_idx_add(True)
            drop_clusterId_add(row['clusterId'])

        cluster_rms_add(rms)
    
    df_cluster['rms']=cluster_rms
    
    df_cluster_drop=pd.DataFrame()
    df_cluster_drop['clusterId']=drop_clusterId
    
    complement  = df_cluster[~df_cluster['clusterId'].isin(df_cluster_drop['clusterId'])]
       
    return complement


def meanStateFilter(df, index, row, frame_state_asteroid='icrf', frame_state_observer='icrf', **kwargs):
    """Initial orbit determination from a set of Right Ascension and Declination observations.
    
    Parameters:
    -----------
    df                 ... Pandas DataFrame containing nightly RA and DEC [deg], time [JD, MJD] UTC,
                           heliocentric observer positions and velocities [au]
    index, row         ... Pandas DataFrame index and row containing clustered observation Ids, 
                           epoch and mean state
                           of cluster
                      

    Keyword arguments:
    ------------------
    frame_state_asteroid  ... coordinate frame of asteroid/mean state: 'ecliptic' or 'icrf' 
    frame_state_observer  ... coordinate frame of observer state: 'ecliptic' or 'icrf'
    output_units          ... units for O-C results: 'deg', 'arcsec', 'rad'
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
    epoch             ... epoch of best fit orbit
    state             ... state of best fit orbit (x,y,z,vx,vy,vz) [au, au/day]
    """
    

    epoch=row['cluster_epoch']
    state=row[['x_a','y_a','z_a','vx_a','vy_a','vz_a']].to_numpy().astype('float')
    [rms, dra, ddec] = ephem.radecResiduals(df, epoch, state, output_units='arcsec',
                                            frame_state_ast=frame_state_asteroid, 
                                            frame_state_obs=frame_state_observer, **kwargs)
    
    return rms, epoch, state


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
     
    cdf2=(cdf.iloc[cdf.astype(str).drop_duplicates(
                   subset='obsId',keep="first").index]).reset_index(drop=True)        
    return cdf2

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

    
    cdf_idx = range(0,len(cdf))
    cdf2 = cdf.reset_index(drop=True)
 
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
                        #print(i,j,is_subset)
                        if(is_subset):
                            subset_clusters_app(i)
                            subset_cluster_ids_app([i,j])
                            break

    idx_to_drop = subset_clusters
    #print(idx_to_drop)
    
    cdf3 = cdf2.drop(index=idx_to_drop)
    return cdf3, subset_clusters, subset_cluster_ids 




##############################################################################
#
# If Thor package is available filtering can be performed with Gauss method
#
##############################################################################
def filterClustersThor(df_obs, df_cluster, rms_max=200):
    """Filter clusters of observations via Initial Orbit Determination from Thor package
    
    Parameters:
    -----------
    df_obs     ... Pandas DataFrame containing obsId, epoch an RA DEC values for all observations
    df_cluster ... Pandas DataFrame containing lists of linked observations
    
    
    Keyword Arguments:
    -------------------
    rms_max    ... RMS O-C cutoff for filter
    ncores     ... number of cores
    
    Returns:
    --------
    df_cluster_filtered      ... filtered Pandas DataFrame containing lists of linked observations
    """
    
    df_c = deduplicateClusters(df_cluster)
    
    keep_cluster=[]
    keep_cluster_add=keep_cluster.append
    
    cluster_rms=[]
    cluster_rms_add=cluster_rms.append
    
    for index, row in df_c.iterrows():
        df=df_obs[[x in row['obsId'] for x in df_obs['obsId']]]
        
        #df=df_obs[df_obs['obsId'].isin(row['obsId'])]        
        rms = -999
        try:
            [rms, epoch, orbit] = iodFilter(df)

            if(rms<rms_max):
                keep_cluster_add(row['clusterId']) 
                
        except:
            pass

        cluster_rms_add(rms)
    
    df_cluster['rms']=cluster_rms
    
    df_cluster_filtered  = df_cluster[df_cluster['clusterId'].isin(keep_cluster)]
    
#     df_cluster_filstered=df_cluster[[x in l_dict for x in df_cluster['clusterId']]]    
    return df_cluster_filtered

def iodFilter(df, **kwargs):
    """Initial orbit determination from a set of Right Ascension and Declination observations.
    
    Parameters:
    -----------
    df                ... Pandas DataFrame containing nightly RA and DEC [deg], time [JD, MJD] UTC,
                          heliocentric ecliptic observer positions and velocities [au]
    
    Returns:
    --------
    rms               ... root mean square (RMS) of RA and DEC O-Cs [arcseconds]
    epoch             ... epoch of best fit orbit
    state             ... state of best fit orbit (x,y,z,vx,vy,vz) [au, au/day]
    
    Requires:
    ---------
    THOR iod
    
    """
    
    [idx, threeobs] = select3obs(df, return_df=True, method='max_arc')
    
    coords_eq_ang = threeobs[['RA','DEC']].values
    t = threeobs['time'].values
    coords_obs = threeobs[['x_obs','y_obs','z_obs']].values
    
    gauss_sol = iod.gaussIOD(coords_eq_ang, t, coords_obs, 
                           velocity_method='gibbs', light_time=True, iterate=True, 
                           mu=const.GM, max_iter=10, tol=1e-15)
       
    rms_all = np.zeros(3)
    i=0
    for sol in gauss_sol:
        [rms, dra, ddec] = ephem.radecResiduals(df, sol[0], sol[1:7], output_units='arcsec')
        rms_all[i]=rms 
        i=i+1   
    
    i_min=np.argmin(rms_all)
    
    # print best state with lowest RMS
    return rms_all[i_min], gauss_sol[i_min][0], gauss_sol[i_min][1:7]
 
def select3obs(df, method='max_arc', return_df=True, **kwargs):
    """Select three observations from an observation dataframe.
    
    Parameters:
    -----------
    df        ... Pandas DataFrame containing nightly RA and DEC [deg], time [JD, MJD] UTC,
                  heliocentric ecliptic observer positions and velocities [au]
                  
    Keyword Arguments:
    ------------------
    method    ... method for selecting three observations: 
                  'max_arc': maximise observation arc
                  'random': uniform random sampling
    return_df ... boolean (default True): return dataframe conainting data from three observations only
                
    Returns:
    --------
    idx    ... DataFrame Indices for selection 
    
    """
    if(method == 'max_arc'):
        idx=[]
        idx_add=idx.append
        tmin=df['time'].min()
        tmax=df['time'].max()
        tcenter=(tmax+tmin)/2
        
        idx_add(df[(df['time']==df['time'].min())].index[0])
        idx_add(df.iloc[(df['time']-tcenter).abs().argsort()[0:1]].index[0])
        idx_add(df[(df['time']==df['time'].max())].index[0])
    
    elif(method == 'random'):
        idx=np.random.choice(df.index,3, replace=False)
        
    else:
        raise Exception("Error in select3obs: Unknown selection method. ")
    
    if(return_df):
        return idx, df.loc[idx]
    else:
        return idx     
