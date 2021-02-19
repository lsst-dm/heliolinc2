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
pandas_utils

LSST Solar System Processing

Utility functions for pandas DataFrames
Implementation: Python 3.6, S. Eggl 20210218
"""

import pandas as pd
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as ius

############################################
# Global Variables
###########################################

au=149597870.700

data2['re_au']=data2['AstRange(km)']/au

data2['rs_au']=np.sqrt(data2['Ast-Sun(J2000x)(km)']**2+data2['Ast-Sun(J2000y)(km)']**2+data2['Ast-Sun(J2000z)(km)']**2)/149597870.700

############################################
# MODULE SPECIFIC EXCEPTION
###########################################

class Error(Exception):
    """Module specific exception."""
    pass


############################################
# Functions
###########################################

__all__ = ['xyz2r', 'derivativesFromSpline']


def xyz2r(df, xName='Ast-Sun(J2000x)(km)', 
          yName='Ast-Sun(J2000y)(km)', 
          zName='Ast-Sun(J2000z)(km)'):
    """Calculate Euclidean Distance from 3 component vector.
    
    Parameters:
    -----------
    df       ... input dataframe
    objName  ... name of individual object/trajectory in dataframe
    xName    ... column name for x coordinate
    yName    ... column name for y coordinate
    zName    ... column name for z coordinate
    
    Returns:
    -------- 
    r        ... numpy array, Euclidean distance
    """
    r=np.sqrt(df[xName]**2+df[yName]**2+df[zName)

    return r                                       
                                           
def derivativesFromSpline(df,objName='obj',xName='rs_au',tName='time'):
    """Calculate first and second derivatives of a colum in a DataFrame
    via splines.
    
    Parameters:
    -----------
    df       ... input dataframe
    objName  ... name of individual object/trajectory in dataframe
    xName    ... column name for dependent variable in spline interpolation
    tname    ... column name for independent variable (such as time).
    
    Returns:
    -------- 
    x,xdot,xddot ... numpy array, returns x,dx/dt,d^2x/dt^2 for all values of t in df
    """
    
    i=0
    objlist=df[objName].unique()
    xdot=[]
    xddot=[]
    x=[]
    for o in objlist:
        dat=df[df[objName]==o]
        s=ius(dat[tName].values, dat[xName].values)
        for index, row in dat.iterrows():
            x.append(s.derivatives(row[tName])[0])
            xdot.append(s.derivatives(row[tName])[1])
            xddot.append(s.derivatives(row[tName])[2])
    return np.array(x),np.array(xdot),np.array(xddot)


