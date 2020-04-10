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
vector

LSST Solar System Processing

Standard vector operations utilizing numba.
Implementation: Python 3.6, S. Eggl 20191010
"""

# Accelerators
import numpy as np
import numba


__all__ = ['norm', 'unitVector', 'dot2D','rotateVector',
           'sphereLineIntercept']

############################################
# MODULE SPECIFIC EXCEPTION
###########################################
class Error(Exception):
    """Vector module specific exception."""
    
    pass

############################################
# VECTOR OPERATIONS
###########################################
@numba.njit
def norm(v):
    """Calculate 2-norm for vectors 1D and 2D arrays.

    Parameters:
    -----------
    v ... vector or 2d array of vectors

    Returns:
    --------
    u ... norm (length) of vector(s)

    """

    if(v.ndim == 1):
        n = np.vdot(v, v)
    elif(v.ndim == 2):
        lv = len(v[:, 0])
        n = np.zeros(lv)
        for i in range(lv):
            n[i] = np.vdot(v[i, :], v[i, :])
    else:
        raise TypeError

    return np.sqrt(n)


@numba.njit
def unitVector(v):
    """Normalize vectors (1D and 2D arrays).

    Parameters:
    -----------
    v ... vector or 2d array of vectors

    Returns:
    --------
    u ... unit lenght vector or 2d array of vectors of unit lenght

    """

    if(v.ndim == 1):
        u = v/norm(v)
    elif(v.ndim == 2):
        lv = len(v[:, 0])
        dim = len(v[0, :])
        u = np.zeros((lv, dim))
        for i in range(lv):
            n = norm(v[i, :])
            for j in range(dim):
                u[i, j] = v[i, j]/n
    else:
        raise Exception('Error in function unit_vector: \
                         only 1D and 2D arrays supported')

    return u


@numba.njit
def dot2D(r,v,cols_r,cols_v):
    """Calculate r.v where both are 2D numpy arrays 
    over user specified columns.

    Parameters:
    -----------
    r      ... 2d array of vectors containing r
    v      ... 2d array of vectors containing v
    cols_r ... columns in r array to be used
    cols_v ... columns in v array to be used

    Returns:
    --------
    u ... 1d array of r.v values

    """

    lv = len(v[:, 0])
    dimr = len(cols_r)
    dimv = len(cols_v)
    
    if(dimr != dimv):
        raise Exception('Error in function dot_2d: input arrays \
                        have to have the same dimensions')
    u = np.zeros(lv)

    for i in range(lv):
        for j in range(dimr):
            u[i] = u[i] + r[i,cols_r[j]]*v[i,cols_v[j]]
    return u

def rotateVector(angle, axis, vector, deg=True):
    """Rotate vector about arbitrary axis by an angle.

    Parameters:
    -----------
    angle  ... rotation angle 
    axis   ... rotation axis: numpy array (n,3)
    vector ... vector to be rotated: numpy array(n,3)
    deg    ... True: angles given in [deg], 
               False: angles given in [rad]
               
    Returns:
    --------
    vrot   ... rotated vector
    """
    if(deg):
        angl=np.deg2rad(angle)
    else:
        angl=angle
        
    sina = np.sin(angl)
    cosa = np.cos(angl)
    
    u = unitVector(axis)
    uxv = np.cross(u, vector)
    uuv = u*(np.vdot(u, vector))
    vrot = uuv.T+cosa*np.cross(uxv, u).T+sina*uxv.T
    return vrot.T

def sphereLineIntercept(l, o, r):
    """Calculate intercept point between line y = l.x + o
       and a sphere around the center of origin of the
       coordinate system: c=0

       Parameters:
       ------------
       l ... vector of line of sight
       o ... observer position
       r ... vector of distances (radii on the heliocecntric sphere)

       Returns:
       --------
       x_intercept  ... position of intersection
                        between the line and the sphere
                        along the line of sight
    """

    x_intercept = np.full(np.shape(l), np.nan)

    ln = unitVector(l)

    for i in range(len(l[:, 0])):

        lo = np.vdot(ln[i, :], o[i, :])

        if (r.size == 1):
            r2 = r*r
        else:
            r2 = r[i]*r[i]

        discrim = lo**2 - (np.vdot(o[i, :], o[i, :]) - r2)
        if(discrim >= 0):
            # line and sphere do actually intersect
            d = -lo + np.sqrt(discrim)
            x_intercept[i, :] = o[i, :]+d*ln[i, :]

    return x_intercept