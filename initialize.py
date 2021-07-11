############File containing the parameters for LJ-Sheet04
"""
Created on Fri May 21 11:33:12 2021
@author: Phili
"""


#Import fundamential packages
import numpy as np
import scipy.constants as sc
from tqdm import tqdm #step counter
import matplotlib
import matplotlib.pyplot as plt
import time #To count the Time

#Setup initial parameters
#######Number of Particles############
n = 2                           # Number of Molecules along one axis
######################################

sigma = 3.21e-10                # [m]
dim = 3                         # dimensionality of the system
N = (n**3)*3                    # Total number of Particles,
M = n**3                        # Number of molecules

global fac_dens
mO = 2.6561e-26                 # Oxygen particle mass [kg]
mH = 1.6735e-2                  # Hydrogen particle mass [kg]
T0 = 300                        # Initial Temperature
tfac = 0.5                      # Scaling for the time
dt = 1e-15 * tfac               # duration of one step [s]
#####Calculated parameters######

fac_dens = 1
number_dens = fac_dens*sigma**(-dim) # number densitx [1/m^3]
V = M/number_dens             # Volume of the Box [m^3]
L = V**(1/3)                  # Side length of the cubic box [m]
facsig = 2.5                 # cut off scaling for sigma             
rcut = facsig*sigma           # cut off distance
Td = T0                     # target temperature



def initial_setup(N,L,dim):
    """ Creates a intitial Box with N particles ordered
    in a 3D grid. The N particles form water molecules.
    
    N (int): Number of particles
    L (float): Box size
    """
    r = np.zeros((3,N))
    #n = int(np.rint((N/3)**(1/dim)))
    #nz = int(n+n)
    x1 = np.linspace(0+L/n/2,L-L/n/2,n)
    d0 = 1.0e-10 #intial distance between hydrogen and oxygen
    x2 = x1 + d0
    x3 = x1 - d0/np.sqrt(2)
    z = np.zeros(3*n)
    xy = np.zeros(3*n)
    for ii in range(n):
        z[ii*3+1] = x1[ii]
        z[ii*3+2] = x2[ii]
        z[ii*3]   = x3[ii]
        xy[ii*3] = x1[ii]
        xy[ii*3+1] = x1[ii]
        xy[ii*3+2] = x1[ii]
    rx, ry,rz = np.meshgrid(x1,x1,z)
    r[0,:] = np.reshape(rx,N)
    r[1,:] = np.reshape(ry,N)
    r[2,:] = np.reshape(rz,N)
    r[0,0::3] += d0/np.sqrt(2) #add additional offset to one Oxygen
    return r

def initial_masses(N,M,mO,mH):
    """Creates initial masses
    """
    masses = np.zeros((N))
    for ii in range(M):
        masses[ii*3] = mH
        masses[ii*3+1] = mO
        masses[ii*3+2] = mH
    return masses

#setup velocity
def initial_vel(T,N,masses,dim):
    """ Creates initial velocities, corresponding to a Temperature. 
    The velocities are normal distributed and not Maxwell-Boltzmann
    
    Args:
        T (float): Temperature
        N (int): Number of particles
        m (float): Mass of the partiles
    """
    rand_v = np.random.uniform(-1,1,(N,dim)) #random velocities
    sum_v = np.sum(rand_v, axis = 0)/N #Total velocity
    rand_v = rand_v - sum_v #Subtract cm-Motion
    K = np.sum(np.sum(rand_v**2,axis = 1)*masses)/2 #Kinetic energy
    T_act = 2/3*K/sc.k/N #Temperature
    R_scale = T/T_act #Scaling factor
    v = np.transpose(rand_v*np.sqrt(R_scale)) #rescalte to reach T
    return v


def write_init_pos(r0,N):
    """Writes a snapshot of the initial configuration to 
    a file, that can be opened by ovito
    """
    types = np.linspace(0,N-1,N)
    f = open('initial_pos.dump', 'w')
    for kk in range(1):
        f.write('ITEM: TIMESTEP \n')
        f.write('{} \n'.format(dt*kk/3600/24))
        f.write('ITEM: NUMBER OF ATOMS \n')
        f.write('{} \n'.format(N))
        f.write('ITEM: BOX BOUNDS pp pp pp \n')
        f.write('{} {} \n'.format(-0,L))
        f.write('{} {} \n'.format(-0,L))
        f.write('{} {} \n'.format(-0,L))
        f.write('ITEM: ATOMS id type x y z Radius \n')
        for ii in range(N):
            f.write(' {} {} {} {} {} {} \n'.format(ii+1,types[ii],r0[0,ii],r0[1,ii],r0[2,ii], .2e-10 ))
        f.close() 
        
        
def write_traj(name,r_eq):
    """ Function which can write a trajectory to a file.
    This file can be opened by ovito
    
    Args:
        name (str): name of the file
        r_eq: positions of the particles
    """
    f = open(name, 'w') #eqilibration.dump'
    N =len(r_eq[0,:,0])
    steps = len(r_eq[0,0,:])
    types = np.linspace(0,N-1,N)
    types = np.ones(N)
    types[1::3] = 2
    for kk in tqdm(range(steps)):
        f.write('ITEM: TIMESTEP \n')
        f.write('{} \n'.format(dt*kk))
        f.write('ITEM: NUMBER OF ATOMS \n')
        f.write('{} \n'.format(N))
        f.write('ITEM: BOX BOUNDS pp pp pp\n')
        f.write('{} {} \n'.format(-0,L))
        f.write('{} {} \n'.format(-0,L))
        f.write('{} {} \n'.format(-0,L))
        f.write('ITEM: ATOMS id type x y z Radius \n')
        for ii in range(N):
            f.write(' {} {} {} {} {} {}\n'.format(ii+1,types[ii],r_eq[0,ii,kk],r_eq[1,ii,kk],r_eq[2,ii,kk], .2e-10, ))
    f.close()    
    return         
        