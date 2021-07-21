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
difi_analysis

LSST Solar System Processing

Routines for checking the performance of HelioLinC2 using difi
Implementation: Python 3.6, J. Moeyens, S. Eggl 20200409
"""

import os
import sys
import numpy as np
import pandas as pd

import difi

# Utility functions
from . import utility as ut

# COLUMN_MAPPING = {
#     # difi column name : data column name
#     "linkage_id" : "cluster_Id",
#     "obs_id" : "obs_Id",
#     "truth" : "obj",
#     "night" : "night",
#     "time" : "time"
# }


__all__ = ['obs2difi','linkages2difi', 'runDifi']


def obs2difi(df, uniqueObsId=False, uniqueObsIdName='obsId', 
             objIdName="ObjID",fieldIdName='observationId',timeName="FieldMJD",
             raName="AstRA(deg)", decName="AstDec(deg)"):
    
    """Converts observations in pandas DataFrame to difi format.
    
    Parameters:
    df              ... pandas DataFrame containing observations 
    uniqueObsId     ... True/False: If a unique observation identifier 
                        is present in the dataframe use that, 
                        otherwise create your own
    uniqueObsIdName ... what is the name of the unique observation identifier?
    objIdName       ... name of object identifier column
    fieldIdName     ... name of observed field identifier
    timeName        ... name of time identifier column
    raName          ... name of Right Ascension column
    decName         ... name of Declination column
    """

    # Rename columns to work with difi
    column_mapping={objIdName: "obj", fieldIdName: "FieldId",
                         timeName:"time",raName:"RA", decName:"DEC"}
    df.rename(columns=column_mapping, inplace=True)
    
    # Sort by observation time
    df.sort_values(by=['time'], inplace=True, na_position='last')
    
     # Create observation ID as string (difi requirement)
    df['obsId']=df.index()
    # If a unique observation Identifier is present in the dataframe use that
    if (uniqueObsId):
        df['obs_Id'] = df[uniqueObsIdName].values.astype(str)
    # Otherwise create your own
    else:    
        df['obs_Id'] = df['obsId'].values.astype(str)
    
    df['night']=ut.lsstNight(df['time'],df['time'].min())
    
    unique_objects=list(df.groupby('obj').groups.keys())
    
    uodf=pd.DataFrame(unique_objects,columns=['obj'])
    
    uodf['objId']=uodf.index
    
    df=df.merge(uodf,on='obj',how='left')
    
    return


def linkages2difi(df,clusterId_name='cluster_Id',observationId_name='obsId',output='pandas'):
    """Map pandas dataframe with linkages to difi format"""

    o = observationId_name
    c = clusterId_name
    rw = []
    for index, row in df.iterrows():
        colo=row[o].astype(str)
        colc=np.full(len(colo),row[c],dtype=int)
        rw.append(np.array([colc,colo]).T)
    difi=np.concatenate(rw)
    
    if (output=='pandas' or output =='pd'):
        difi=pd.DataFrame(difi,columns=[c,o])
    
    return difi


def runDifi(observationdf, linkingdf, obsIdName='obsId', linkageIdName='cluster_Id',
            objIdName='obj', nightName='night', timeName='time',
            findability='tracklet',  linkage_min_obs=2, max_obs_separation=1.5/24, 
            min_linkage_nights=3, 
            metric="nightly_linkages", classes=None):
    
    """Convert observations and linkage dataframe to difi input and run difi.
    
    Parameters:
    -----------
    observationdf         ... observations DataFrame
    linkingdf             ... linking 
    obsIdName
    linkageIdName
    findability
    linkage_min_obs
    max_obs_separation
    min_linkage_nights
    metric
    classes

    
    Returns:
    --------
    all_truths                ... DataFrame with
    findable_observations     ... DataFrame with findable observations
    summary                   ... Summary DataFrame of observations  
    all_linkages_heliolinc    ... DataFrame with info on linkages 
    all_truths_heliolinc      ... DataFrame with
    summary_heliolinc         ... Summary DataFrame on linking performance
    findable_objects          ... DataFrame containing findable objects only
    missed_objects            ... DataFrame containing missed objects only
    """
    
    # create difi column mapping dict
    column_mapping={}
    column_mapping["linkage_id"] = linkageIdName
    column_mapping["obs_id"] = obsIdName
    column_mapping["truth"] = objIdName
    column_mapping["night"] = nightName
    column_mapping["time"] = timeName
    
    
    # difi column name : data column name
<<<<<<< HEAD
=======

>>>>>>> 072b56ecaedd85274c273ecf041ce5eb4c8efdaa
#     "linkage_id" : "cluster_Id",
#     "obs_id" : "obs_Id",
#     "truth" : "obj",
#     "night" : "night",
#     "time" : "time"
    
    
   
    #dfobs['obs_Id']=dfobs['obs_Id'].convert_dtypes()
    
    # Analyze observations and determine what is findable
    all_truths, findable_observations, summary = difi.analyzeObservations(
        observationdf,
        metric="nightly_linkages",
        linkage_min_obs=2,          # a tracklet should be at least 2 observations
        max_obs_separation=1.5/24,  # these observations should be within 90 minutes
        min_linkage_nights=3,       # we need a tracklet on 3 unique nights
        column_mapping=column_mapping)
    
    # Reshape HelioLinC3D cluster output to be compatible with difi
    #linkages2difi(df,clusterId_name='cluster_Id',observationId_name='obsId',output='pandas')
    ldifi=linkages2difi(linkingdf,column_mapping["linkage_id"],column_mapping['obs_id'])
    #ldifi.rename(columns={column_mapping['obs_id']:'obs_id'},inplace=True)
  
    # Run difi on linkages to determine what has been found 
    all_linkages_heliolinc, all_truths_heliolinc, summary_heliolinc = difi.analyzeLinkages(
        observationdf, 
        ldifi, 
        all_truths=all_truths,
        classes=classes,
        min_obs=6, 
        contamination_percentage=20, 
        column_mapping=column_mapping)
    
    
    findable_objects = all_truths_heliolinc[(all_truths_heliolinc['findable']==1)] 
    
    missed_objects = findable_objects[findable_objects['found']==0]
    
    
    return (all_truths, findable_observations, summary, 
           all_linkages_heliolinc, all_truths_heliolinc, summary_heliolinc,
           findable_objects, missed_objects)


class Error(Exception):
    """Module specific exception."""
    
    pass


# def _convertLinkagesToDifi(observations, linkages, columnMapping=COLUMN_MAPPING):
#     ### Siegfried: you may want to remove this once you have cluster ids unique in your code
#     linkages["clusterId"] = np.arange(0, len(linkages))
#     linkages["clusterId"] = linkages["clusterId"].astype(str)
    
#     ### This is only specific to the JPL 3 month dataset: we ideally don't want this 
#     ### bit of code to exist in this function. We should expect the user to provide 
#     ### correctly formatted observations.
#     duplicatedIDs = ["NS", "FD"]
#     for j in duplicatedIDs:
#         mask = observations[columnMapping["truth"]].isin([j]) 
#         newIDs = np.array(["{}_{:08d}".format(t, i) for i, t in enumerate(observations[mask][columnMapping["truth"]].values)])
#         observations.loc[mask, columnMapping["truth"]] = newIDs

#     observations.loc[observations["obj"].str.match("^S1"), "class"] = "MBA"
#     observations.loc[observations["obj"].str.match("^S0"), "class"] = "NEO"
#     observations.loc[observations["obj"].str.match("^FD"), "class"] = "NOISE 1"
#     observations.loc[observations["obj"].str.match("^NS"), "class"] = "NOISE 2"
#     ### END
    
#     # Extract the columns we would want in the allLinkages dataframe
#     #allLinkages = linkages[["clusterId", "r", "drdt", 'cluster_epoch', 'x_ecl', 'y_ecl', 'z_ecl','vx_ecl', 'vy_ecl', 'vz_ecl']]
#     allLinkages = linkages
    
#     # Extract the columns we would want in the linkageMembers dataframe
#     linkageMembers_temp = linkages[["clusterId", "obsId"]]
    
#     # Split each linkage into its different observation IDs
#     linkage_list = linkageMembers_temp[columnMapping["obs_id"]].str.strip("[").str.strip("]").str.split().tolist()

#     # Build initial DataFrame
#     linkageMembers = pd.DataFrame(pd.DataFrame(linkage_list, 
#                                                index=linkageMembers_temp["clusterId"].values).stack(), 
#                                   columns=[columnMapping["obs_id"]])

#     # Reset index 
#     linkageMembers.reset_index(1, drop=True, inplace=True)

#     # Make linkage_id its own column
#     linkageMembers[columnMapping["linkage_id"]] = linkageMembers.index

#     # Re-arrange column order 
#     linkageMembers = linkageMembers[[columnMapping["linkage_id"], columnMapping["obs_id"]]]

#     # Not all linkages have the same number of detections, empty detections needs to be dropped
#     linkageMembers[columnMapping["obs_id"]].replace("", np.nan, inplace=True)
#     linkageMembers.sort_values(by=["clusterId", "obsId"], inplace=True)
#     linkageMembers.dropna(inplace=True)
#     linkageMembers.reset_index(drop=True, inplace=True)
    
#     return observations, linkageMembers, allLinkages

# def calcNight(mjd, midnight=0.166):
#     night = mjd + 0.5 - midnight
#     return night.astype(int)

# def calcFindableMOPS(observations, 
#                      trackletMinObs=2, 
#                      trackMinNights=3, 
#                      columnMapping=COLUMN_MAPPING):
#     # Groupby night, then count number of occurences per night
#     night_designation_count = observations.groupby(["nid"])[columnMapping["truth"]].value_counts()
#     night_designation_count = pd.DataFrame(night_designation_count)
    
#     #db
#     print(night_designation_count)
#     #edb
    
#     night_designation_count.rename(columns={columnMapping["truth"]: "num_obs"}, inplace=True)
#     night_designation_count.reset_index(inplace=True)
    
#     #db
#     print(night_designation_count)
#     #edb

    
#     # Remove nightly detections that would not be linked into a tracklet
#     night_designation_count = night_designation_count[night_designation_count["num_obs"] >= trackletMinObs]

#     # Groupby object then count number of nights
#     try: 
#         designation_night_count = pd.DataFrame(night_designation_count.groupby([columnMapping["truth"]])["nid"].value_counts())
#     except:
#         # No objects satisfy the requirements, return empty array
#         return np.array([])
#     designation_night_count.rename(columns={"nid": "num_nights"}, inplace=True)
#     designation_night_count.reset_index(inplace=True)

#     # Grab objects that meet the night requirement
#     tracklet_nights_possible = designation_night_count[columnMapping["truth"]].value_counts()
#     return tracklet_nights_possible.index[tracklet_nights_possible >= trackMinNights].values

# def analyzePerformanceMinObs(observations, linkages, minObs=3, columnMapping=COLUMN_MAPPING):

#     # Convert data frames into the difi format
#     observations, linkageMembers, allLinkages = _convertLinkagesToDifi(observations, linkages, columnMapping=columnMapping)
    
#     # Create the allTruths and preliminary summary dataframes, flag all objects with more 
#     # than minObs observations as findable
#     allTruths, summary = difi.analyzeObservations(observations,
#                                                   minObs=minObs,
#                                                   classes="class",
#                                                   columnMapping=columnMapping)


#     # Find the objects that have enough observations to make 3 tracklets of at least 2 detections
#     # on 3 unique nights
#     observations["nid"] = calcNight(observations[columnMapping["epoch_mjd"]])
#     findable_by_mops = calcFindableMOPS(
#         observations, 
#         trackletMinObs=2, 
#         trackMinNights=3, 
#         columnMapping=columnMapping
#     )
#     allTruths_mops = allTruths.copy()
#     allTruths_mops["findable"] = np.zeros(len(allTruths), dtype=int)
#     allTruths_mops.loc[allTruths_mops[columnMapping["truth"]].isin(findable_by_mops), "findable"] = 1
    
#     # Analyze linkages with non-MOPS findability criterion
#     allLinkages, allTruths, summary = difi.analyzeLinkages(
#         observations, 
#         linkageMembers,  
#         allTruths=allTruths,
#         classes="class",
#         minObs=minObs, 
#         contaminationThreshold=0, 
#         verbose=False,
#         columnMapping=columnMapping
#     )
    
#     # Analyze linkages with MOPS findability criterion
#     allLinkages_mops, allTruths_mops, summary_mops = difi.analyzeLinkages(
#         observations, 
#         linkageMembers,  
#         allTruths=allTruths_mops,
#         classes="class",
#         minObs=minObs, 
#         contaminationThreshold=0, 
#         verbose=False,
#         columnMapping=columnMapping
#     )
    
#     # Merge results into a single summary dataframe
#     summary = summary.merge(
#     summary_mops[["class", "completeness", "findable", "findable_found",
#                   "findable_missed", "not_findable_found", "not_findable_missed"]], 
#         on="class", 
#         suffixes=("", "_mops")
#     )
#     summary = summary[[
#            'class', 'completeness', 'completeness_mops', 'findable', 'found', 'findable_found',
#            'findable_missed', 'not_findable_found', 'not_findable_missed',
#            'findable_mops', 'findable_found_mops', 'findable_missed_mops',
#            'not_findable_found_mops', 'not_findable_missed_mops',
#            'linkages', 'pure_linkages', 'pure_complete_linkages',
#            'partial_linkages', 'mixed_linkages', 'unique_in_pure',
#            'unique_in_pure_complete', 'unique_in_partial',
#            'unique_in_pure_and_partial', 'unique_in_pure_only',
#            'unique_in_partial_only', 'unique_in_mixed'
#     ]]
    
#     # Merge allTruths into a single dataframe
#     allTruths = allTruths.merge(allTruths_mops[["obj", "findable"]], on=columnMapping["truth"], suffixes=("", "_mops"))
#     allTruths = allTruths[[columnMapping["truth"], 'num_obs', 'findable', 'findable_mops', 'found_pure', 'found_partial', 'found']]
    
#     return observations, linkageMembers, allLinkages, allTruths, summary

# def analyzePerformance(observations, linkages, 
#                        sspFindability={"trackletMinObs":2, "trackMinNights":2},
#                        mopsFindability={"trackletMinObs":2, "trackMinNights":3}, minObs=4,
#                        columnMapping=COLUMN_MAPPING):

#     # Convert data frames into the difi format
#     observations, linkageMembers, allLinkages = _convertLinkagesToDifi(observations, linkages, columnMapping=columnMapping)
    
#     # Create the allTruths and preliminary summary dataframes, flag all objects with more 
#     # than minObs observations as findable
#     allTruths, summary = difi.analyzeObservations(observations,
#                                                   minObs=minObs,
#                                                   classes="class",
#                                                   columnMapping=columnMapping)
    
#     # Find the objects that have enough observations to make 2 tracklets of at least 2 detections
#     # on 3 unique nights
#     observations["nid"] = calcNight(observations[columnMapping["epoch_mjd"]])
#     findable_by_ssp = calcFindableMOPS(
#         observations, 
#         columnMapping=columnMapping,
#         **sspFindability
#     )
    
#     allTruths["findable"] = np.zeros(len(allTruths), dtype=int)
#     allTruths.loc[allTruths[columnMapping["truth"]].isin(findable_by_ssp), "findable"] = 1

#     # Find the objects that have enough observations to make 3 tracklets of at least 2 detections
#     # on 3 unique nights
#     findable_by_mops = calcFindableMOPS(
#         observations, 
#         columnMapping=columnMapping,
#         **mopsFindability
#     )
    
#     allTruths_mops = allTruths.copy()
#     allTruths_mops["findable"] = np.zeros(len(allTruths), dtype=int)
#     allTruths_mops.loc[allTruths_mops[columnMapping["truth"]].isin(findable_by_mops), "findable"] = 1
    
#     # Analyze linkages with non-MOPS findability criterion
#     allLinkages, allTruths, summary = difi.analyzeLinkages(
#         observations, 
#         linkageMembers,  
#         allTruths=allTruths,
#         classes="class",
#         minObs=minObs, 
#         contaminationThreshold=0, 
#         verbose=False,
#         columnMapping=columnMapping
#     )
    
#     # Analyze linkages with MOPS findability criterion
#     allLinkages_mops, allTruths_mops, summary_mops = difi.analyzeLinkages(
#         observations, 
#         linkageMembers,  
#         allTruths=allTruths_mops,
#         classes="class",
#         minObs=minObs, 
#         contaminationThreshold=0, 
#         verbose=False,
#         columnMapping=columnMapping
#     )
    
#     # Merge results into a single summary dataframe
#     summary = summary.merge(
#     summary_mops[["class", "completeness", "findable", "findable_found",
#                   "findable_missed", "not_findable_found", "not_findable_missed"]], 
#         on="class", 
#         suffixes=("", "_mops")
#     )
#     summary = summary[[
#            'class', 'completeness', 'completeness_mops', 'findable', 'found', 'findable_found',
#            'findable_missed', 'not_findable_found', 'not_findable_missed',
#            'findable_mops', 'findable_found_mops', 'findable_missed_mops',
#            'not_findable_found_mops', 'not_findable_missed_mops',
#            'linkages', 'pure_linkages', 'pure_complete_linkages',
#            'partial_linkages', 'mixed_linkages', 'unique_in_pure',
#            'unique_in_pure_complete', 'unique_in_partial',
#            'unique_in_pure_and_partial', 'unique_in_pure_only',
#            'unique_in_partial_only', 'unique_in_mixed'
#     ]]
    
#     # Merge allTruths into a single dataframe
#     allTruths = allTruths.merge(allTruths_mops[["obj", "findable"]], on=columnMapping["truth"], suffixes=("", "_mops"))
#     allTruths = allTruths[[columnMapping["truth"], 'num_obs', 'findable', 'findable_mops', 'found_pure', 'found_partial', 'found']]
    
#     return observations, linkageMembers, allLinkages, allTruths, summary





# def df2difi(df,index_name,value_name):
#     """Map pandas dataframe with lists of values to THOR difi format"""

#     difi=df[value_name].apply(pd.Series) \
#     .merge(df, right_index = True, left_index = True) \
#     .drop([value_name], axis = 1) \
#     .melt(id_vars = [index_name], value_name = value_name) \
#     .drop("variable", axis = 1) \
#     .dropna() \
#     .sort_values(by=[index_name]) \
#     .astype('int') \
#     .reset_index(drop=True)
    
#     return difi
<<<<<<< HEAD
=======
  
>>>>>>> 072b56ecaedd85274c273ecf041ce5eb4c8efdaa

### This is all that needs to be run

# DATA_DIR = "/epyc/projects/lsst_ssp/difi/test_data"
# # Read observations into a dataframe
# observations = pd.read_csv(
#     os.path.join(DATA_DIR, "Veres_5x5deb_14nights.csv"), 
#     sep=",", 
#     dtype={
#         "obsId" : str,
#         "obj" : str,
#         "objId" : str,
#     },
#     index_col=False
# )
# # Read SSP linkages
# linkages = pd.read_csv(
#     os.path.join(DATA_DIR, "heliolinc2_clusters_14nights_5x5_mean_state_filter_rms_3arcsec.csv"), 
#     sep=",", 
#     index_col=False
# )

# # Run analysis
# observations, linkageMembers, allLinkages, allTruths, summary = analyzePerformance(observations, linkages)