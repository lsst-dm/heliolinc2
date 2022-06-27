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
utility

LSST Solar System Processing utility routines.


Implementation: S. Eggl 20191215
"""

# Accelerators
import numpy as np

# Time transforms
from astropy.time import Time

__all__ = ['whichGroup','whichNight','obs2heliolinc']


############################################
# UTILITY FUNCTIONS
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

    try:
        print(required_input_columns.values())
        for col in required_input_columns.values():
            if col not in df.columns:
                raise Exception("Error in obs2heliolinc: required input column not present:",col)
                
    
        icols = [required_input_columns[r] for r in required_input_columns]
        ocols = [required_output_columns[r] for r in required_output_columns]
    
        # Create a copy with only required inputs
        df2 = df[icols]
    
    except:
        raise Exception("Error in obs2heliolinc: not all required input columns are present:",required_input_columns)
    
    # Rename columns of observations DataFrame to work with HelioLinC3D
    column_mapping = dict(zip(icols,ocols))  
    df2.rename(columns = column_mapping, inplace=True) 
    
    # Transform time from UTC to TDB
    df2['time'] = Time(df2['time'].values,scale='utc',format='mjd').tdb.mjd
  
    # Sort by observation time
    # df2.sort_values(by=['time','RA'], inplace=True)
    
    df2.reset_index(inplace=True, drop=True)
    
    # Create observation ID         
    df2['obsId'] = df2.index
    
    # Create a unique observation identifier if none is present
    if 'obsName' not in df2.columns:
        df2['obsName'] = df2['obsId'].values.astype(str)
      
    # df2['night']=(ut.lsstNight(df2['time'],df2['time'].min())).astype(int)
    
    # Check if all required columns are present
    for col in required_output_columns:
           if col not in df2.columns:
                raise Exception("Error in obs2heliolinc: required columns not present.")
               
    
    return df2


def whichGroup(obsTimeMJD):

    """Group time sorted observations in consecutive pairs of two to construct arrows, i.e. on sky positions and velocities.
    
    Parameters:
    -----------
    obsTimeMJD ... array of observation times in Modified Julian Day format
    
    Returns:
    --------
    group   ... array of integers designating the group an observation belongs to
    
    """            
    group = []
    groupId = 0
    ttold = min(obsTimeMJD)
    j = 0
    for tt in obsTimeMJD:
        if(tt>ttold and j==0):
            j = 1
            ttold = tt
        elif(tt>ttold and j==1):
            groupId = groupId+1
            ttold = tt
            j = 0
        
        group.append(groupId)

    return np.array(group).astype(int)


def whichNight(expMJD, minexpMJD=0, local_midnight=0.166):
    """Calculate the night for a given observation epoch and a survey start date.
    
    Parameters:
    -----------
    expMJD ... epoch of observation / exposure [modified Julian date, MJD]
    minexpMJD ... start date of survey [modified Julian date, MJD]
    
    Returns:
    --------
    night ... the night of a given observation epoch with respect to the survey start date.
    
    """
    const = minexpMJD + local_midnight - 0.5
    night = np.floor(expMJD - const).astype(int)
    return night