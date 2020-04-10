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
performancecheck

LSST Solar System Processing

Routines for checking the performance of HelioLinC2 
Implementation: Python 3.6, S. Eggl 20200310
"""

# Accelerators
import numpy as np

# Pandas DataFrames
import pandas as pd

# Constants such as the speed of light and GM
from . import constants as cnst


__all__ = ['discoverableObjects', 'clusterPurity', 'objectsInClusters',
           'correctPairs', 'observationsInCluster', 'observationsInArrows',
           'obs2heliocentricArrows',
          ]


def discoverableObjects(df, min_n_nights, min_n_obs_per_night):
    """Which asteroids are discoverable?
    
    Parameters:
    -----------
    df                      ... Simulated observations dataframe containing observations, night, and objId
    min_n_nights            ... minimum number of nights with min_n_obs_per_night observations required for discovery
    min_n_obs_per_night     ... minimum number of observations per night required for discovery
    
    Returns:
    --------
    n_disc      ... number of objects that should be discoverable given the observational data in df
    disc_objId  ... a list of object Ids of discoverable objects 
    """
    
    obs_grouped_by_objId=df.groupby([ 'objId','night'])
    objId_night=np.array(list(obs_grouped_by_objId['night'].groups.keys()))
    objId_night_Nobs=objId_night[(obs_grouped_by_objId['obj'].count() >= min_n_obs_per_night)]
    nnights_per_obj_with_geNobs=np.bincount(objId_night_Nobs[:,0])
    
    #n_disc=len(nnights_per_obj_with_geNobs[(nnights_per_obj_with_geNobs>= min_n_nights)])
    
    objs_to_be_discovered=np.where(nnights_per_obj_with_geNobs >= min_n_nights)[0]
    # the first two obects in the JPL database are FD and NS
    disc_objId=df['obj'][(df['objId'].isin(objs_to_be_discovered[2:]))].unique()
    
    n_disc=len(disc_objId)

    return n_disc, disc_objId

def clusterPurity(objects_in_cluster_df):
    """Determine purity of clusters.
    
    Parameters:
    -----------
    objects_in_cluster_df ...Pandas DataFrame containing objID and clusterId (as index)
    
    Returns:
    --------
    n_pure     ... number of pure clusters containing only observations of one object
    percentage ... percentage with respect to number of clusters: len(objects_in_cluster_df)
    n_noise    ... number of noise clusters
    """
    #len1=[]
    #len1_add=len1.append
    n_pure=0
    n_noise=0
    for index, row in objects_in_cluster_df.iterrows():
        #print(row[0], len(row[0]))
        #if(type(row[0]) is list):
            if(len(row[0]) == 1):
                if(row[0]=='FD' or row[0]=='NS'):
                    n_noise=n_noise+1
                else:
                    n_pure=n_pure+1
            else:
                    n_noise=n_noise+1
        
    #n_pure=len(len1)
    percentage = np.round(n_pure/len(objects_in_cluster_df.index)*100,2)
    return n_pure, percentage, n_noise

def objectsInClusters(obs_df, cluster_df, fileio=False):
    """Convert observations in a cluster to unique object names.

    Parameters:
    -----------
    obs_df      ... Pandas DataFrame containing observations and object IDs
    cluster_df  ... Pandas DataFrame containing observations in each cluster

    Returns:
    --------
    ofdf.index                ... index / clusterId of cluster_df
    objects_in_cluster   ... list of object designations in cluster
    """
    objects_in_cluster = []
    objects_in_cluster_add = objects_in_cluster.append

    for index, row in cluster_df.iterrows():
        if(fileio):
            obj=obs_df['obj'][(obs_df['obsId'].isin(convertFileIo(row['obsId'])))].unique()
        else:
            obj=obs_df['obj'][obs_df['obsId'].isin(row['obsId'])].unique()
            
        objects_in_cluster_add(obj)
        
    obj_ic=np.array(objects_in_cluster)

    objects_in_cluster_df=pd.DataFrame(obj_ic, columns=['objId'])

    return len(objects_in_cluster_df.index), objects_in_cluster_df


def convertFileIo(string):
    """Convert each line of linkages read in from a file to a list.
    
    Parameters:
    -----------
    string  ... line of linkages read in from file into Pandas DataFrame
    
    Returns:
    --------
    li      ... list of linkages
    
    """
    string2=string.replace('[',"")
    string3=string2.replace(']',"")
    li = list(string3.split(" "))
    for ele in li:
        if (ele == ''):
            li.remove('')
    return li

def correctPairs(df,pairs):
    """Which pairs are actually good, i.e. which tracklets are real?
    
    Parameters:
    -----------
    df      ... pandas dataframe with observations
    pairs   ... list of pairs of observations [obsid1,obsid2] linked into arrows
    
    Returns:
    --------
    correct ... logical array (dimension of pairs[:,0])
    """
    p=np.array(pairs)
    #find out which pairs are actually good
    pair_obj=np.array([df['objId'][p[:,0]].values,df['objId'][p[:,1]].values]).T
    #print(pair_obj)
    #correct=np.where(df['obj'][p[:,0]].values == df['obj'][p[:,1]].values)
    correct=np.where(pair_obj[:,0] == pair_obj[:,1])
    return correct[0]

def observationsInCluster(df, pairs, cluster, garbage=False):
    """List observations in each cluster.
    
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
        unique_labels = np.unique(cluster.labels_)
    else:
        unique_labels = np.unique(cluster.labels_)[1:]
    #number of clusters
    n_clusters = len(unique_labels)
    
    #which objects do observations in pairs (tracklets) belong to
    p = np.array(pairs)              
    
    obs_in_cluster=[]
    obs_in_cluster_add=obs_in_cluster.append
    
    #cluster contains 
    for u in unique_labels:
        #which indices in pair array appear in a given cluster?
        idx = np.where(cluster.labels_ == u)[0]
        
        #which observations are in this cluster
        obs_in_cluster_add(np.unique(p[idx].flatten()))
    
    return obs_in_cluster, unique_labels


# def objectsInCluster(df,pairs,cluster):
#     """List observations in each cluster.
    
#     Parameters:
#     -----------
#     df      ... pandas dataframe with observations
#     pairs   ... list of pairs of observations [obsid1,obsid2] linked into arrows
#     cluster ... output of clustering algorithm (sklearn.cluster)
    
#     Returns:
#     --------
#     obs_in_cluster ... list of objects in each cluster
#     """
#     #cluster names (beware: -1 is the cluster of all the leftovers)
#     unique_labels = np.unique(cluster.labels_)
#     #number of clusters
#     n_clusters = len(unique_labels)
    
#     #which objects do observations in pairs (tracklets) belong to
#     p = np.array(pairs)              
#     pair_obj = np.array([df['obj'][p[:,0]].values,df['obj'][p[:,1]].values]).T
    
#     obj_in_cluster=[]
#     obj_in_cluster_add=obj_in_cluster.append
#     #cluster contains 
#     for u in unique_labels:
#         #which indices in pair array appear in a given cluster?
#         idx = np.where(cluster.labels_ == u)[0]
#         #find unique object ids in cluster
#         uniq_obj=np.unique(pair_obj[idx])
#         obj_in_cluster_add(uniq_obj)
        
#     return obj_in_cluster, unique_labels


def observationsInArrows(df,goodpairs,*args,**kwargs):
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
    df_obs_in_arrows=pd.DataFrame(goodpairs,**kwargs)
    return df_obs_in_arrows


def obs2heliocentricArrows(df, r, drdt, tref, lttc=False, verbose=True):
    """Create tracklets/arrows from dataframe containing nightly RADEC observations
       and observer positions.

    Parameters:
    -----------
    df       ... Pandas DataFrame containing nightly RA and DEC [deg], time [JD, MJD],
                 (x,y,z)_observer positions [au, ICRF]
    r        ... assumed radius of heliocentric sphere used for arrow creation[au]
    drdt     ... assumed radial velocity
    tref     ... reference time for arrow generation. Used to calculate how much the 
                 heliocentric distance changes between observations based on assumed dr/dt


    Keyword arguments:
    ------------------
    lttc (optional)        ... light travel time correction
    verbose (optional)     ... print verbose progress statements  

    Returns:
    --------
    x         ... tracklet/arrow position (3D) [au]
    y         ... tracklet/arrow velocity (3D) [au]
    t         ... tracklet/arrow reference epoch [JD/MJD]
    """
    
    # speed of light in au/day
    c_aupd = 173.145

    # Transform RADEC observations into positions on the unit sphere (US)
    xyz = tr.radec2icrfu(df['RA'], df['DEC'], deg=True)

    # Those are the line of sight (LOS) vectors
    los = np.array([xyz[0], xyz[1], xyz[2]]).T

    # Use the position of the observer and the LOS to project the position of
    # the asteroid onto a heliocentric great circle with radius r
    observer = df[['x_obs', 'y_obs', 'z_obs']].values

    # Calculate how much the heliocentric distance changes
    # during the obsevations based on assumed dr/dt
    dt = tref-df['time'].values
    dr = drdt*dt
    r_plus_dr = r+dr

    # Heliocentric postions of the observed asteroids
    posu = ls.sphere_line_intercept(los, observer, r_plus_dr)

    if(verbose):
        print('Heliocentric positions generated.')
 
    
    # tracklet position for filtered pairs
    x = posu[:-1,:]
    # tracklet time
    t = df['time'].values
    # tracklet velocity through forward differencing
    va = []
    vapp = va.append
    dt = t[1:]-t[0:-1]
    dx = posu[1:,:]-posu[0:-1,:]
    for d in range(0,3):
        vapp(np.divide(dx[:,d],dt))
    v = np.array(va).T
    t= df['time'].values[:-1]
    # correct arrows for light travel time
    if(lttc):
        if(verbose):
            print('(Linear correction for light travel time aberration...')
        xo = observer[:-1,:]
        dist = vec.norm(x-xo)
        xl = x.T-dist/c_aupd*v.T
        return xl.T, v, t

    else:
        return x, v, t

