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

COLUMN_MAPPING = {
    # difi column name : data column name
    "linkage_id" : "clusterId",
    "obs_id" : "obsId",
    "truth" : "obj",
    "epoch_mjd" : "time",
}

def _convertLinkagesToDifi(observations, linkages, columnMapping=COLUMN_MAPPING):
    ### Siegfried: you may want to remove this once you have cluster ids unique in your code
    linkages["clusterId"] = np.arange(0, len(linkages))
    linkages["clusterId"] = linkages["clusterId"].astype(str)
    
    ### This is only specific to the JPL 3 month dataset: we ideally don't want this 
    ### bit of code to exist in this function. We should expect the user to provide 
    ### correctly formatted observations.
    duplicatedIDs = ["NS", "FD"]
    for j in duplicatedIDs:
        mask = observations[columnMapping["truth"]].isin([j]) 
        newIDs = np.array(["{}_{:08d}".format(t, i) for i, t in enumerate(observations[mask][columnMapping["truth"]].values)])
        observations.loc[mask, columnMapping["truth"]] = newIDs

    observations.loc[observations["obj"].str.match("^S1"), "class"] = "MBA"
    observations.loc[observations["obj"].str.match("^S0"), "class"] = "NEO"
    observations.loc[observations["obj"].str.match("^FD"), "class"] = "NOISE 1"
    observations.loc[observations["obj"].str.match("^NS"), "class"] = "NOISE 2"
    ### END
    
    # Extract the columns we would want in the allLinkages dataframe
    allLinkages = linkages[["clusterId", "r", "drdt", 'cluster_epoch', 'x_a', 'y_a', 'z_a','vx_a', 'vy_a', 'vz_a']]
    
    # Extract the columns we would want in the linkageMembers dataframe
    linkageMembers_temp = linkages[["clusterId", "obsId"]]
    
    # Split each linkage into its different observation IDs
    linkage_list = linkageMembers_temp[columnMapping["obs_id"]].str.strip("[").str.strip("]").str.split().tolist()

    # Build initial DataFrame
    linkageMembers = pd.DataFrame(pd.DataFrame(linkage_list, index=linkageMembers_temp["clusterId"].values).stack(), columns=[columnMapping["obs_id"]])

    # Reset index 
    linkageMembers.reset_index(1, drop=True, inplace=True)

    # Make linkage_id its own column
    linkageMembers[columnMapping["linkage_id"]] = linkageMembers.index

    # Re-arrange column order 
    linkageMembers = linkageMembers[[columnMapping["linkage_id"], columnMapping["obs_id"]]]

    # Not all linkages have the same number of detections, empty detections needs to be dropped
    linkageMembers[columnMapping["obs_id"]].replace("", np.nan, inplace=True)
    linkageMembers.sort_values(by=["clusterId", "obsId"], inplace=True)
    linkageMembers.dropna(inplace=True)
    linkageMembers.reset_index(drop=True, inplace=True)
    
    return observations, linkageMembers, allLinkages

def calcNight(mjd, midnight=0.166):
    night = mjd + 0.5 - midnight
    return night.astype(int)

def calcFindableMOPS(observations, 
                     trackletMinObs=2, 
                     trackMinNights=3, 
                     columnMapping=COLUMN_MAPPING):
    # Groupby night, then count number of occurences per night
    night_designation_count = observations.groupby(["nid"])[columnMapping["truth"]].value_counts()
    night_designation_count = pd.DataFrame(night_designation_count)
    night_designation_count.rename(columns={columnMapping["truth"]: "num_obs"}, inplace=True)
    night_designation_count.reset_index(inplace=True)

    # Remove nightly detections that would not be linked into a tracklet
    night_designation_count = night_designation_count[night_designation_count["num_obs"] >= trackletMinObs]

    # Groupby object then count number of nights
    try: 
        designation_night_count = pd.DataFrame(night_designation_count.groupby([columnMapping["truth"]])["nid"].value_counts())
    except:
        # No objects satisfy the requirements, return empty array
        return np.array([])
    designation_night_count.rename(columns={"nid": "num_nights"}, inplace=True)
    designation_night_count.reset_index(inplace=True)

    # Grab objects that meet the night requirement
    tracklet_nights_possible = designation_night_count[columnMapping["truth"]].value_counts()
    return tracklet_nights_possible.index[tracklet_nights_possible >= trackMinNights].values

def analyzePerformance(observations, linkages, minObs=3, columnMapping=COLUMN_MAPPING):

    # Convert data frames into the difi format
    observations, linkageMembers, allLinkages = _convertLinkagesToDifi(observations, linkages, columnMapping=columnMapping)
    
    # Create the allTruths and preliminary summary dataframes, flag all objects with more 
    # than minObs observations as findable
    allTruths, summary = difi.analyzeObservations(observations,
                                                  minObs=3,
                                                  classes="class",
                                                  columnMapping=columnMapping)


    # Find the objects that have enough observations to make 3 tracklets of at least 2 detections
    # on 3 unique nights
    observations["nid"] = calcNight(observations[columnMapping["epoch_mjd"]])
    findable_by_mops = calcFindableMOPS(
        observations, 
        trackletMinObs=2, 
        trackMinNights=3, 
        columnMapping=columnMapping
    )
    allTruths_mops = allTruths.copy()
    allTruths_mops["findable"] = np.zeros(len(allTruths), dtype=int)
    allTruths_mops.loc[allTruths_mops[columnMapping["truth"]].isin(findable_by_mops), "findable"] = 1
    
    # Analyze linkages with non-MOPS findability criterion
    allLinkages, allTruths, summary = difi.analyzeLinkages(
        observations, 
        linkageMembers,  
        allTruths=allTruths,
        classes="class",
        minObs=minObs, 
        contaminationThreshold=0, 
        verbose=False,
        columnMapping=columnMapping
    )
    
    # Analyze linkages with MOPS findability criterion
    allLinkages_mops, allTruths_mops, summary_mops = difi.analyzeLinkages(
        observations, 
        linkageMembers,  
        allTruths=allTruths_mops,
        classes="class",
        minObs=minObs, 
        contaminationThreshold=0, 
        verbose=False,
        columnMapping=columnMapping
    )
    
    # Merge results into a single summary dataframe
    summary = summary.merge(
    summary_mops[["class", "completeness", "findable", "findable_found",
                  "findable_missed", "not_findable_found", "not_findable_missed"]], 
        on="class", 
        suffixes=("", "_mops")
    )
    summary = summary[[
           'class', 'completeness', 'completeness_mops', 'findable', 'found', 'findable_found',
           'findable_missed', 'not_findable_found', 'not_findable_missed',
           'findable_mops', 'findable_found_mops', 'findable_missed_mops',
           'not_findable_found_mops', 'not_findable_missed_mops',
           'linkages', 'pure_linkages', 'pure_complete_linkages',
           'partial_linkages', 'mixed_linkages', 'unique_in_pure',
           'unique_in_pure_complete', 'unique_in_partial',
           'unique_in_pure_and_partial', 'unique_in_pure_only',
           'unique_in_partial_only', 'unique_in_mixed'
    ]]
    
    # Merge allTruths into a single dataframe
    allTruths = allTruths.merge(allTruths_mops[["obj", "findable"]], on=columnMapping["truth"], suffixes=("", "_mops"))
    allTruths = allTruths[[columnMapping["truth"], 'num_obs', 'findable', 'findable_mops', 'found_pure', 'found_partial', 'found']]
    
    return observations, linkageMembers, allLinkages, allTruths, summary


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