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
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 12}

matplotlib.rc('font', **font)

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

import pandas as pd


def plot_field(dfobs):
    """Plot RADEC observations in field.
    
    Parameters:
    -----------
    dfobs ... Pandas DataFrame containing RADEC coordinates and object designations

    Returns:
    --------
    displays scatter plots
    """
    # collect data from different categories
    dfobs_FD=dfobs.iloc[dfobs['obj'].values=='FD']
    dfobs_NS=dfobs.iloc[dfobs['obj'].values=='NS']
    dfobs_obj=dfobs.drop(dfobs_FD.index).drop(dfobs_NS.index)
    
    
    fig = plt.figure(dpi=300, figsize=(6,6))
    ax = fig.add_subplot(1,1,1)
    ax.scatter(dfobs_NS.RA, dfobs_NS.DEC, label='NS', s=.1, marker='.', color='g')
    ax.scatter(dfobs_FD.RA, dfobs_FD.DEC, label='FD', s=.1, marker='.')
    ax.scatter(dfobs_obj.RA, dfobs_obj.DEC, label='Objects', s=.1, marker=',')
    plt.xlabel('RA [deg]')
    plt.ylabel('DEC [deg]')
    ax.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show()
    
    return

def plot_field2(dfobs1, dfobs2, label1='data1', label2='data2'):
    """Plot RADEC observations in field for two different pandas DataFrames
    
    Parameters:
    -----------
    dfobs1, dfobs2 ... Pandas DataFrames containing RADEC coordinates and object designations
    
    Keyword Arguments:
    ------------------
    label1, label2 ... labels for DataFrames in plot

    Returns:
    --------
    displays scatter plots
    """
    
    fig = plt.figure(dpi=300, figsize=(6,6))
    ax = fig.add_subplot(1,1,1)
    ax.scatter(dfobs1.RA, dfobs1.DEC, label=label1, s=.1, marker='.', color='r')
    ax.scatter(dfobs2.RA, dfobs2.DEC, label=label2, s=.1, marker='.', color='b')
    plt.xlabel('RA [deg]')
    plt.ylabel('DEC [deg]')
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