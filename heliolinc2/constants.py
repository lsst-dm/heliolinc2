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
constants

LSST Solar System Processing

Constants.
Implementation: Python 3.6, S. Eggl 20200220
"""

############################################
# CONSTANTS
###########################################
# speed of light in au/day
CAUPD = 173.145

# Gaussian Gravitational Constant [au^1.5/Msun^0.5/D]
GAUSSK = 0.01720209894846

# default gravitational parameter [Gaussian units]
GM = GAUSSK*GAUSSK

# Julian date of J2000 epoch
J2000_JD = 2451545.0

# Earth obliquity [deg]
EARTH_OBLIQUITY = 23.4392911111112




