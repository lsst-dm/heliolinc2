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

__all__ = [ 'lsstNight']


############################################
# UTILITY FUNCTIONS
###########################################

def lsstNight(expMJD, minexpMJD=0, local_midnight=0.166):
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
    night = np.floor(expMJD - const)
    return night