#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 21:19:52 2021

@author: philipp
"""

#velocity verlet
import numpy as np
from numba import jit
import scipy.constants as sc

#Implementation of the half step velocity verlet with thermostat
def vel_half(v,F,m,dt,e):
    v = (v+F/m*dt/2)/(1+e*dt/2)
    return v

def eta(e,vhalf,m,dt,Td):
    e = e + G(vhalf,m,dt,Td)*dt
    return e

def G(v,m,dt,Td,Mfac = 50):
    dtc = Mfac*dt
    N = len(v[0,:])
    g = 3*N
    Q = g*sc.k*Td*dtc**2
    K = np.sum(v**2*m)/2
    G = 1/Q*(2*K-3*N*sc.k*Td)
    return G
#@jit
def r_full(r,vhalf,dt):
    """Velocity verlet algorythmus. The positions calculation.

    Args:
        r (float matrix 3,N): Positions of each particle
        vhalf (float matrix 3,N): Velcoity each particle
    """
    r = r+vhalf*dt #update position
    return r

#@jit
def vel_full(vhalf,Fbar,m,dt,ebar):
    """Update the particle velocity
    
    Args:
        v (float matrix 2,N): velocities of the particles
        F (float matrix 2,N): forces of the previous step
        F_bar (float matrix 2,N): force of the actual constellation
        dt (float): duration of one time step
        e (float): eta factor
    """
    v =  vhalf+ dt/2*(Fbar/m-ebar*vhalf)
    return(v)