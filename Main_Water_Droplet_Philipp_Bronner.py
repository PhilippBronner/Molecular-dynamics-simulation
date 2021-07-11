"""
Philipp Bronner, 11.07.2021

This file executes the Simulation of Water molecules. 
This is performed in an NVE ensemble using various interatomic forces (see Forces_Distances.py)
A file for OVITO containing the trajectory (position of every 10th step) is stored.
"""


# Importing simulation packages
import simulation10 as sim #The external simulation functions
from initialize import *  #Setup initial conditions
import scipy.constants as sc


########Initial_setup##########
r0 = initial_setup(N,L,dim) #setup the Positions
masses = initial_masses(N, M, mO, mH) #setup a mass array
v0 = initial_vel(T0,N,masses,dim) #setup velocities

write_init_pos(r0,N)  #Write the intitial Positions to a file

steps_Eq = 1000 #steps for equilibration
steps_Prod = 5000 #steps for production
Td = 300
##########Equilibration and Production run################ 
req, rend, vend =  sim.simulation_prod(r0,v0,N,L,masses,rcut,dt,T0,steps_Eq,Td) #equilibriation
rprod, rend, vend =  sim.simulation_prod(rend,vend,N,L,masses,rcut,dt,T0,steps_Prod,Td) #production
write_traj('trajectory_6P.dump', rprod)
