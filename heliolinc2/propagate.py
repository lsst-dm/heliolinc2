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
propagate

LSST Solar System Processing

Transformations between common coordinate and time systems
Implementation: Python 3.6, S. Eggl 20191115
"""
# Accelerators
import numpy as np
import numba
from multiprocessing import Pool

from joblib import Parallel, delayed 

# NASA NAIF Spice wrapper 
import spiceypy as sp

# Constants such as the speed of light and GM
from . import constants as cnst
from . import transforms as tr


__all__ = ['propagateState', 'propagate2body', 'propagateLinear', 'stateFromOrbit']


def propagateState(x, v, t, tp, n_jobs=1, propagator='linear'):
    """Propagation of states with choce of propagator.

    Parameters:
    -----------
    x            ... array of 3D positions
    v            ... array of 3D velocities
    t            ... array of epochs for states (x,v)
    tp           ... epoch to propagate state to
    n_jobs       ... number of worker processes for 2body propagation
    propagator   ... select propagator from 'linear, 2body, nbody'
   
    
    Returns:
    --------
    xp           ... array of propagated 3D positions
    vp           ... array of propagated 3D velocities
    dt           ... array of time deltas wrt the propatation epoch: tp-t
    """

    if(propagator == 'linear'):
        [xp, vp, dt] = propagateLinear(x, v, t, tp)
        
    elif(propagator == '2body'):
        [xp, vp, dt] = propagate2body(x, v, t, tp, n_jobs)
        
    elif(propagator == 'nbody'):
        raise Exception('N-body propagation not yet implemented.')
        
    else:
        raise Exception('Error in Propagate Arrows: select valid propagator.')
    
    return xp, vp, dt

def prop2b(x, v, dt):
    state_in = np.hstack((x,v))
    #print(state_in)
    state_out = sp.prop2b(cnst.GM, state_in, dt)
    #print(state_out)
#     print(state_out)
    return state_out

def propagate2body(x, v, t, tp, n_jobs):
    """ Propagate states to the same time using spicepy's 2body propagation.

    Parameters:
    -----------
    x            ... array of 3D heliocentric/barycentric positions
    v            ... array of 3D heliocentric/barycentric velocities
    t            ... array of epochs for states (x,v)
    tp           ... epoch to propagate state to
    n_jobs       ... number of worker processes for 2body propagation

    Returns:
    --------
    xp           ... array of propagated 3D positions
    vp           ... array of propagated 3D velocities
    dt           ... array of time deltas wrt the propatation epoch: tp-t
    """
    
    dt = np.array(tp-t)
   
#     print('x.ndim',x.ndim)
#     print('x.shape',x.shape)
#     print('x',x)
#     print('v',v)
    
    if (len(x)<1):
        # empty
        xp =[]
        vp =[]
        
    elif(x.ndim==1):
        state = sp.prop2b(cnst.GM, np.hstack((x, v)), dt)
        xp = state[0:3]
        vp = state[3:6]

    elif(x.ndim==2):
        lenx = len(x[:, 0])
        #state = np.zeros((lenx,6))
        if(n_jobs>1):
            #print('n_jobs:',n_jobs)
            #batchsize=np.rint(lenx/(n_jobs)/20).astype(int)  
            #batchsize=10000*n_jobs
            #print('batch_size:', batchsize, ' of ', lenx)
            
#             results = Parallel(n_jobs=n_jobs, 
#                                batch_size=batchsize,
#                                pre_dispatch='2*n_jobs'
#                                )(delayed(prop2b)(x[i,:],v[i,:],dt[i]) 
#                                for i in range(lenx))  

            results = Parallel(n_jobs=n_jobs, 
                               batch_size='auto',
                               pre_dispatch='2*n_jobs'
                               )(delayed(prop2b)(x[i,:],v[i,:],dt[i]) 
                               for i in range(lenx))  
            state = np.array(results)
            #print('state',state)
        else:
            state = np.zeros((lenx,6))
            for i in range(lenx):
                state[i, :] = sp.prop2b(cnst.GM, np.hstack((x[i, :], v[i, :])), dt[i])[0:6]
        
        #print('state', state)     
        xp = state[:, 0:3]
        vp = state[:, 3:6]
    else:
        raise TypeError
        
    return xp, vp, dt    
    

def propagateLinear(x, v, t, tp):
    """Linear propagattion of arrows to the same time.

    Parameters:
    -----------
    x  ... array of 3D positions
    v  ... array of 3D velocities
    t  ... array of epochs for states (x,v)
    tp ... epoch to propagate state to

    Returns:
    --------
    xp ... array of propagated 3D positions
    v  ... array of velocities
    dt ... array of time deltas wrt the propatation epoch: tp-t
    """
    dt = tp-t
    
    if (len(x)<1):
        # empty
        xp =[]
        
    elif(x.ndim==1):
        xp = x + v*dt
        
    elif(x.ndim==2):    
        xp = x + (v*np.array([dt, dt, dt]).T)
        
    return xp, v, dt


def stateFromOrbit(orbital_elements, epochs, target_epoch, element_type='Cometary', propagator='2body', frame='icrf'):
    """Calculate Cartesian state at a given target epoch from orbital elements
    
    Parameters:
    -----------
    orbital_elements  ... set (or sets) of orbital elements [n, 6]
    epochs            ... reference epochs of orbital elements 
    target_epoch      ... target epoch for propagation
    
    Keyword Arguments:
    ------------------
    element_type      ... type of orbital elements ('Cometary', 'Keplerian')
    propagator        ... type of propagation used ('2body', 'linear')
    frame             ... frame for Cartesian coordinates ('icrf', 'ecliptic')
    
    Returns:
    --------
    pstates           ... array of propagated Cartesian states (x,y,z,vx,vy,vz)
    
    External Libraries:
    -------------------
    heliolinc2 as hl
    """
    
    pstate=[]
    if(element_type == 'Cometary'):
        ele2state = tr.cometary2cartesian
    elif(element_type == 'Keplerian'):
        ele2state = tr.keplerian2cartesian
    else:
        raise Exception('Error in stateFromOrbit: unknown element type')
    
    state=[]
    for i in range(len(epochs)):
        state.append(ele2state(epochs[i],orbital_elements[i,:], frame='ecliptic'))
    
    states=np.array(state)
    pstate=propagateState(states[:,0:3], states[:,3:6], epochs, target_epoch, propagator=propagator)              
   # print(pstate)        
    psecl=np.hstack([pstate[0],pstate[1]])
    
    if(frame == 'icrf'):
        pstates = tr.ecliptic2icrf(psecl)
    elif(frame == 'ecliptic'):    
        pstates = psecl
    else:
        raise Exception('Error in stateFromOrbit: specified frame unknown')

    return pstates
# def propagate2body(x, v, t, tp):
#     """ Propagate states to the same time using spicepy's 2body propagation.

#     Parameters:
#     -----------
#     x  ... array of 3D heliocentric/barycentric positions
#     v  ... array of 3D heliocentric/barycentric velocities
#     t  ... array of epochs for states (x,v)
#     tp ... epoch to propagate state to

#     Returns:
#     --------
#     xp ... array of propagated 3D positions
#     vp ... array of propagated 3D velocities
#     dt ... array of time deltas wrt the propatation epoch: tp-t
#     """
#     dt = np.array(tp-t)
    
# #     print('x.ndim',x.ndim)
# #     print('x.shape',x.shape)
# #     print('x',x)
# #     print('v',v)
    
#     if (len(x)<1):
#         # empty
#         xp =[]
#         vp =[]
        
#     elif(x.ndim==1):
#         state = sp.prop2b(cnst.GM, np.hstack((x, v)), dt)
#         xp = state[0:3]
#         vp = state[3:6]

#     elif(x.ndim==2):
#         lenx = len(x[:, 0])
#         dimx = len(x[0, :])
#         dimv = len(v[0, :])
#         try:
#             assert(dimx == dimv)
#         except(TypeError):
#             raise TypeError

#         xp = []
#         xp_add = xp.append
#         vp = []
#         vp_add = vp.append
#         for i in range(lenx):
#             state = sp.prop2b(cnst.GM, np.hstack((x[i, :], v[i, :])), dt[i])
#             xp_add(state[0:3])
#             vp_add(state[3:6])

#     else:
#         raise TypeError
        
#     return np.array(xp), np.array(vp), dt
