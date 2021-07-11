#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 16:52:30 2021

@author: philipp
"""

import numpy as np
import scipy.constants as sc

def read_traj(name,steps, N):
    """import the trajectory
        Args:
            name (str): file name of the trajectroy
        return:
            r_traj (3,N,steps): trajectroy matrix
    """
    steps_dat  = steps
    r_traj = np.zeros((3,N,steps_dat))
    f = open(name, 'r')
    lines = f.readlines()
    counter = 0
    step = 0
    for line in lines:
        if counter == 9+N:
            counter = 0
            step += 1
        if counter >= 9:
            pos = line.split()
            particle = int(pos[0])-1
            r_traj[0,particle,step] = float(pos[2])
            r_traj[1,particle,step] = float(pos[3])
            r_traj[2,particle,step] = float(pos[4])
            #print(counter)
        counter += 1
    return r_traj

from tqdm import tqdm #step counter

def rfd(data_traj,Ndr,L):
    """calculated the radial distributoin function
        
        Args:
            data_traj (3,N,steps): matrix with the trajectroy of the particles
            Ndr (float):   Number of steps in which the rdf is calculated
            L (float):     Length of the Box
            steps (float): number of simulation steps
            N (float):     number of particles
        Return:
            r (Ndr): array with the positions where the rdf is evaluated (x axis)
            g (Ndr): array with the values of the rdf
    """
    dr_step = L/2/Ndr #step width of the radius
    steps = len(data_traj[0,0,:])
    N = len(data_traj[0,:,0])
    dist_particle = np.zeros((N,N,steps)) #dummy matrix for the particle-particle distance
    dim = len(data_traj[:,0,0]) #dimensionality of the Simulation
    for ii in tqdm(range(steps)): #calculate the thistances between each particle for each step
        rmat = np.tile(data_traj[:,:,ii],(N,1,1))
        dr = rmat - np.transpose(rmat)
        dr = np.where(abs(dr)<= L/2, dr, dr-L*np.sign(dr))
        dist_particle[:,:,ii] = np.sqrt(np.sum(dr[:,:,:]**2,axis = 1))

    r = np.linspace(dr_step,L/2,Ndr) #x axis
    g = np.zeros((len(r))) # dummy array for the rdf
    rho = N/(L**dim) # partcle density
    
    n, edge = np.histogram(dist_particle, bins = Ndr, range = (0,L/2))
    dn = n/steps/N  #mean number insider r
    DV = 2*np.pi*r**(dim-1)*dr_step #finite width of the shere shell
    g = dn/DV/rho/2 #rdf
    return r[1:],g[1:]

def rfd_Molec(data_traj,Ndr,L):
    """calculated the radial distributoin function
        
        Args:
            data_traj (3,N,steps): matrix with the trajectroy of the particles
            Ndr (float):   Number of steps in which the rdf is calculated
            L (float):     Length of the Box
            steps (float): number of simulation steps
            N (float):     number of particles
        Return:
            r (Ndr): array with the positions where the rdf is evaluated (x axis)
            g (Ndr): array with the values of the rdf
    """
    dr_step = L/2/Ndr #step width of the radius
    steps = len(data_traj[0,0,:])
    N = len(data_traj[0,:,0])
    M = int(N/2)
    dist_particle = np.zeros((M,M,steps)) #dummy matrix for the particle-particle distance
    dim = len(data_traj[:,0,0]) #dimensionality of the Simulation
    r_sp = (data_traj[:,0::2,:]+data_traj[:,1::2,:])/2
    r_sp = np.where(abs(r_sp) < L/2, r_sp, r_sp-L*np.sign(r_sp))
    for ii in tqdm(range(steps)): #calculate the thistances between each particle for each step
        rmat = np.tile(r_sp[:,:,ii],(M,1,1))
        dr = rmat - np.transpose(rmat)
        dr = np.where(abs(dr)<= L/2, dr, dr-L*np.sign(dr))
        dist_particle[:,:,ii] = np.sqrt(np.sum(dr[:,:,:]**2,axis = 1))

    r = np.linspace(dr_step,L/2,Ndr) #x axis
    g = np.zeros((len(r))) # dummy array for the rdf
    rho = M/(L**dim) # partcle density
    
    n, edge = np.histogram(dist_particle, bins = Ndr, range = (0,L/2))
    dn = n/steps/M  #mean number insider r
    DV = 2*np.pi*r**(dim-1)*dr_step #finite width of the shere shell
    g = dn/DV/rho/2 #rdf
    return r[1:],g[1:]


