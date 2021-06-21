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
plot

LSST Solar System Processing

Routines for ploting input and output of HelioLinC2 
Implementation: Python 3.6, S. Eggl 20200410
"""

#Plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection

import pandas as pd

import numpy as np

############################################
# Global Variables
###########################################

# font
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 12}

matplotlib.rc('font', **font)

fontP = FontProperties()
fontP.set_size('xx-small')

# markers
markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 
           'h', 'H', 'D', 'd', 'P', 'X','.',';']
# colors
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k','orange','pink',
          'brown','springgreen','lavender','salmon']

# Label for Right Ascension 
RAlabel = 'RA [deg]'
# Label for Declination
DEClabel = 'DEC [deg]'
# Astronomical Unit in km
au = 149597870.700


############################################
# Functions
###########################################

__all__ = ['plotField', 'plotField2','plotMissedObjectObs',
           'plotCartesianClusteringRadius','plotUniqueObjectsInClusters']

def plotField(dfobs,raName='RA',decName='DEC',objName='obj',noiseName="^NS",
               FalseDetectionName="^FD", plotFalseDetections=False, plotNoise=False,
               save2file=False, filename='missedObjects.png'):
    """Plot RADEC observations in field.
    
    Parameters:
    -----------
    dfobs                 ... Pandas DataFrame containing RADEC coordinates and object designations
    raName                ... name of column in DataFrame containing Right Ascension data (deg or rad)
    decName               ... name of column in DataFrame containing Declination data (deg or rad)
    objName               ... name of column containing object Name
    noiseName             ... name signature of noise data (regexp)
    FalseDetectionName    ... name signature of False Detection data (regexp)
    plotFalseDetections   ... True/False plot False Detections if present in dfobs observation DataFrame
    plotNoise             ... True/False plot noise if present in dfobs observation DataFrame
    safe2file             ... True/False save plot to file
    filename              ... file name for output 
    
    Returns:
    --------
    displays scatter plots
    """
    # collect data from different categories
    dfobs_FD=dfobs[dfobs[objName].str.match( FalseDetectionName)]
    dfobs_NS=dfobs[dfobs[objName].str.match(noiseName)]
    dfobs_obj=dfobs.drop(dfobs_FD.index).drop(dfobs_NS.index)
    
    
    fig = plt.figure(dpi=150, figsize=(4,4))
    ax = fig.add_subplot(1,1,1)
    if (plotNoise):
        ax.scatter(dfobs_NS[raName], dfobs_NS[decName], label='Noise', s=.1, marker='.', color='g')
    if (plotFalseDetections):    
        ax.scatter(dfobs_FD[raName], dfobs_FD[decName], label='False Detections', s=.1, marker='.')
    
    i = 0
    for index, row in dfobs_obj.iterrows():
        data2 = dfobs_obj[dfobs_obj[objName].isin([row[objName]])]
        plt.scatter(data2[raName],data2[decName], s=2, c=colors[np.mod(i,13)], 
                    marker=markers[np.mod(i,16)], label=row[objName])
 
        i = i+1
        
    # ax.scatter(dfobs_obj[raName], dfobs_obj[decName], label='Objects', s=.1, marker=',')
    # plt.xlabel(RAlabel)
    # plt.ylabel(DEClabel)
    # ax.legend(bbox_to_anchor=(1.1, 1.05))
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP,ncol=4)
    plt.xlabel(RAlabel)
    plt.ylabel(DEClabel)
    
    if(save2file):
        plt.savefig(filename)
    else:
        plt.show()
    
    return


def plotMissedObjectObs(dfobs, dfmissed, objName='obj', 
                        raName='RA', decName='DEC', plotAllObs=True,
                        save2file=False, filename='missedObjects.png'):
    """ Plot missed object observations in RA vs DEC.

    Parameters:
    -----------
    dfobs                 ... Pandas DataFrame of observations 
    dfmissed              ... Pandas DataFrame of missed objects
    raName                ... name of column in DataFrame containing Right Ascension data (deg or rad)
    decName               ... name of column in DataFrame containing Declination data (deg or rad)
    objName               ... name of column containing object Name
    plotAllObs            ... True/False plot all observations or only those corresponding to missed objects
    safe2file             ... True/False save plot to file
    filename              ... file name for output 
    
    Returns:
    --------
    matplotlib plot (either saved to file or plt.show()) 
    
    """    
    fontP = FontProperties()
    fontP.set_size('xx-small')
    
    data=dfobs[dfobs[objName].isin(dfmissed[objName])]
    i=0
    
    plt.figure(dpi=300)
    if(plotAllObs):
        plt.scatter(dfobs[raName],dfobs[decName], s=3, c='#CCCCCC')
            
    for index, row in dfmissed.iterrows():
        data2=data[data[objName].isin([row[objName]])]
        plt.scatter(data2[raName],data2[decName], s=2, c=colors[np.mod(i,13)], 
                    marker=markers[np.mod(i,16)], label=row[objName])
        plt.legend(title='Missed Objects', bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP)
        plt.xlabel(RAlabel)
        plt.ylabel(DEClabel)
        i=i+1
        
    if(save2file):
        plt.savefig(filename)
    else:
        plt.show()
        
    return

def plotCartesianClusteringRadius(data, r=0.1, objName='obj', xName='Ast-Sun(J2000x)(km)', 
                                  yName='Ast-Sun(J2000y)(km)',xScale=au, yScale=au, 
                                  title='Objects', xLabel='x [au]',yLabel='y [au]',
                                  save2file=False, filename='clustering.png'):
    """Plot Cartesian position and clustering radius for each observation of every object.
    
    Parameters:
    -----------
    data                  ... Pandas DataFrame of observations 

    objName               ... str, name of column containing object Name
    xName                 ... str, name of column in DataFrame containing x coordinate
    yName                 ... str, name of column in DataFrame containing y coordinate
    xScale                ... float, scaling for x axis 
    yScale                ... float, scaling for y axis  
    xLabel               ... float, scaling for x axis 
    yLabel                ... float, scaling for y axis
    safe2file             ... True/False save plot to file
    filename              ... file name for output 
    
    Returns:
    --------
    matplotlib plot (either saved to file or .show() ) 
    
    """                                
                                  
    patches = []
    for index, row in data.iterrows():
        x1=row[xName]/xScale
        y1=row[yName]/yScale
        circle = Circle((x1, y1), r)
        patches.append(circle)

    fig, ax = plt.subplots(dpi=300)

    i=0
    for index, row in mystery.iterrows():
        data2=data[data[objName].isin([row[objName]])]
        plt.scatter(data2[xName]/xScale,data2[yName]/yScale, s=2, 
                    c=colors[np.mod(i,13)], marker=markers[np.mod(i,16)], label=row[objName])
        plt.legend(title=title, bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP)
        plt.xlabel(xLabel)
        plt.ylabel(yLabel)
        i=i+1

    p = PatchCollection(patches, alpha=0.05)
    ax.add_collection(p)
    
    if(save2file):
        plt.savefig(filename)
    else:
        plt.show()
        
    return
    
    
    
def plotField2(dfobs1, dfobs2, raName='RA',decName='DEC', label1='data1', label2='data2'):
    """Plot RADEC observations in field for two different pandas DataFrames
    
    Parameters:
    -----------
    dfobs1, dfobs2 ... Pandas DataFrames containing RADEC coordinates and object designations
    
    Keyword Arguments:
    ------------------
    raName          ... name of column in DataFrame containing Right Ascension data (deg or rad)
    decName         ... name of column in DataFrame containing Declination data (deg or rad)
    label1, label2  ... labels for DataFrames in plot

    Returns:
    --------
    displays scatter plots
    """
    
    fig = plt.figure(dpi=300, figsize=(6,6))
    ax = fig.add_subplot(1,1,1)
    ax.scatter(dfobs1[raName], dfobs1[decName], label=label1, s=.1, marker='.', color='r')
    ax.scatter(dfobs2[raName], dfobs2[decName], label=label2, s=.1, marker='.', color='b')
    plt.xlabel(RAlabel)
    plt.ylabel(DEClabel)
    ax.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show()
    
    return

def plotUniqueObjectsInClusters(df):
    """Plot histogram of the number of unique objects in each cluster
    
    Parameters:
    -----------
    df ... pandas DataFrame containing objects in cluster per row
    
    Retruns:
    --------
    displays histogram
    """ 
    
    lenof=[len(row) for row in df]
    plt.figure(dpi=150)
    plt.hist(lenof,bins=range(1,20))
    #plt.yscale('log')
    plt.ylim(1,)
    plt.xlim(1,10)
    plt.xlabel('unique objects in cluster')
    plt.ylabel('number of clusters')
    plt.show()
    
    return    