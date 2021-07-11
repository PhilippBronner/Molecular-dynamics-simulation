"""
Created on Fri May  7 22:40:46 2021
@author: philipp Bronner

This routine performes the whole simulation of the Water molecules. 
With given initial conditions.
"""

#Simulation algorythm

#Equilibration
import numpy as np
import Forces_Distances as Force
import Velocity_Verlet_thermostat as VelVer
import scipy.constants as sc
from tqdm import tqdm #step counter
from initialize import *


def simulation_prod(r0,v0,N,L,m,rcut,dt,T0,steps,Td):
    """ Simulation of the LJ fluid process. Starting with a equilibriated system.
    
    Args:
        r0 (np.array (3,N)): initial positions [m]
        v0 (np.array (3,N)): initial velocities [m/s]
        N (int): number of particles
        L (float): length of the box [m]
        m (float): particle mass
        rcut (float): cut off distance for LJ [m]
        dt (float): duration of one time step [s]
        T0 (float): initial Temperature
        steps (int): number of equilibration steps.
    
    Returs: 
        rmat (3,N,steps): positions during the equilibration
        vmat (3,N,steps): velocities ''
        K_obs (steps): Kinetic energy
        U_obs (steps): Potential energy
    """
    dr, rij = Force.dist_rij(r0, N, L, rcut) #initial distance
    F = Force.F_TOT(dr, rij, rcut)           #initial forces
    eta = np.zeros((steps))                  #initial scaling factor in the NVT theromstat.
    eta[0] = 0
    r = r0                                   #initial positions
    v = v0                                   #initial forces
    
    crip = 10 #every 10th step is stored
    r_traj = np.zeros((3,N,int(steps/crip))) #trajectory matrix
    
    for ii in tqdm(range(1,steps)): #loop for the steps, simulation is performed
        #Half step velocity verlet implementation...
        vhalf = VelVer.vel_half(v,F,m,dt,eta[ii-1])
        r = VelVer.r_full(r,vhalf,dt)
        #r = periodicbox(r,L) #check periodic boundary, not used for the water.
        eta[ii] = VelVer.eta(eta[ii-1],vhalf,m,dt,Td)
        dr, rij = Force.dist_rij(r, N, L, rcut)
        
        Fbar = Force.F_TOT(dr, rij, rcut) #Force update

        v = VelVer.vel_full(vhalf,Fbar,m,dt,eta[ii]) #Velocity update
        vcm = np.sum(v*m, axis = 1)/np.sum(m) #top center of mass movement 
        v = v-vcm[:,np.newaxis]
        F = Fbar #Force restet
        if ii%crip == 0:
            r_traj[:,:,int(ii/crip)] = r
    return [r_traj ,r, v]

def periodicbox(pos,L):
    """periodic boundary condition. 
    A particle which leaves the box on one side, enters on the other side.
    
    Args:
        pos (float matrix 2,N) positions of the particles
    """
    #find and correct particles at the left and down border
    pos = np.where(pos <=L, pos, pos%L) 

    #find an correct particles at the upper and right border
    pos = np.where(pos >=0, pos, pos%L)

    return pos


#Optinal functions for temperature observation
def Temp_calc(v,m,N):
    #Caculate the Temperature, dependent on the velocity.
    K = np.sum(v**2*m)/2
    T_act = 2/3*K/sc.k/N
    return T_act


