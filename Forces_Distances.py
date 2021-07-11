"""
Created on Tue May  4 21:12:47 2021
@author: philipp Bronner

This routine contains functions that calculate the Forces acting between Water molecules.
These are:
    - LJ 12-6 interaction between the oxygen atoms
    - Intramolecular interaction between H-O atoms
    - Intramolecular interaction dependingo on the angle between the two
    O-H bonds.
    - Electrostatic interaction between each Atom

Furthermore the distance between each particle is claculated. No PCB is applied
"""

#Importing packages
import numpy as np
import scipy.constants as sc
from numba import jit #to speed up a little bit
from scipy.special import erf
import numpy.linalg as npl #to calculate the norm
from initialize import * #import the mising simulation parameters

def dist_rij(pos, N, box, rcut):
    """Calculates the distances between the particles as vector
    and as absolute value.
    
        Returns:
            rij (N,N): Distances between each particle, with high number on the diagonal axis
            dr (N,3,N): Vector for distances between each particle
    """
    rmat = np.tile(pos,(N,1,1))
    dr = rmat - np.transpose(rmat)
    #No periodic boundary condition
    #dr = np.where(abs(dr)<= box/2, dr, dr-box*np.sign(dr))
    #dr[:,:2,:] = np.where(abs(dr[:,:2,:])<= box/2, dr[:,:2,:], dr[:,:2,:]-box*np.sign(dr[:,:2,:]))
    rij = np.where(abs(dr) < rcut,  dr, 1000000) #implement cut off1
    rij = np.sqrt(np.sum(dr[:,:,:]**2,axis = 1)) #the absolute value of each vector
    rij = rij + np.identity(N)*100000
    return dr, rij
@jit
def Molecule(rij,Fij):
    """Calculates the inramocecular bonding interaction. 
    These are forces between O-H inside each water molecule.
    
        Args:
            rij (N,N): absolute distances between each particle.
            Particle 3*n and 3*n+2 (n is natural number of M) are Hydrogen
            Particle 3*n+1 is the Oxygen. In total there are nMax = M molecules.
            Fij (N,N): absolute-value force matrix. The corresponding values are overwritten with 
            the intramolecular bonding foce,
        Returns:
            Fij (N,N): updated force matrix.
    """
    M = int(len(rij[0,:])/3) #number of Molecules
    for ii in range(M): #Look over the Molecules
        n = ii*3 
        r0 = 1.0e-10 #equilibrium bond length
        kb = 1058e3*4.184/sc.N_A*(1e10)**2 #focre constant
        Fij[n,n+1] += -kb*(rij[n,n+1] - r0) #first bond
        Fij[n+1,n] += -kb*(rij[n+1,n] - r0)
        Fij[n+1,n+2] += -kb*(rij[n+1,n+2]- r0) #second bond
        Fij[n+2,n+1] += -kb*(rij[n+2,n+1]- r0)
    return Fij
    
def Molecule_theta(dr,Fmat):
    """Calculates the intramolecular angle interaction.
    These is aquadratic potential depending on the angle between the twoo O-H bonds 
    of hydrogen.
    
        Args:
            dr (N,3,N): vector matric with the distances between each particle.
            Fmat (N,3,N): vector foce matrix 
        Return;
    """
    kangle = 75e3*4.184/sc.N_A
    theta0 = 104.5/360*2*np.pi
    for ii in range(M):
        n = ii*3
        r1 = dr[n+1,:,n] #O-H bond vector
        r2 = dr[n+1,:,n+2] #O-H bond vector
        a1 = sum(r1*r2)/sum(r2*r2)*r2 #projection of r1 on r2.
        a2 = r1-a1 #perpendicular to r2! (parallel to e_theta)
        b1 = sum(r1*r2)/sum(r1*r1)*r1 #projection of r2 on r1
        b2 = r2-b1 #perpendicular on r1 (paralllel to e_theta)
        theta = np.arccos(r1.dot(r2)/(npl.norm(r1)*npl.norm(r2)))
        Fijvec = kangle*(theta-theta0)*a2/npl.norm(a2)/npl.norm(r2) #Force on particle 2
        Fjivec = kangle*(theta-theta0)*b2/npl.norm(b2)/npl.norm(r1) #Force on particle 1
        Fmat[n,:,n+2] += Fijvec #add the angle foces 
        Fmat[n+2,:,n] += Fjivec
    return Fmat
@jit
def F_TOT(dr, rij, rcut):
    """Total Focre acting on each atom.
    These are:
    - LJ 12-6 interaction between the oxygen atoms
    - Intramolecular interaction between H-O atoms
    - Intramolecular interaction dependingo on the angle between the two
    O-H bonds.
    - Electrostatic interaction between each Atom
    
        Args:
            dr (N,3,N): vector matix that contains the distances between each particles.
            rij (N,N): absolute distance matrix between each particle. 
            Diagonal elements are set to +inf
            rcut: cut off distance for the LJ potential
        Returns: 
            Fi (3,N): Force acting on each particle, depending on its position.
            
    """
    N = len(rij[0,:])
    # First of all the foces that act along rij (er) dicrection.
    Fij = np.zeros((N,N))
    Fij = F_LJ(rij, Fij)
    Fij = Molecule(rij, Fij)/rij
    Fij += F_electrostatic(dr, rij)/rij #devide by rij that the multiplication with dr is a unit vector multplication.
    Fijmat = np.multiply(dr,Fij[:,np.newaxis,:]) #multiply with the distance vector, to get a vector
    #sum up the total Force acting on one particle.
    #the diagonal elements of the Force are interaction of the particle with it self
    #nansum leaves these "nan" elements out.
    #as ling as all forces were previously devided by rij the multiplication with dr is a unit vector.
    
    Fijmat = Molecule_theta(dr, Fijmat) #Force that acts along e_theta
    Fi = np.nansum(Fijmat, axis = 0) 
    return Fi
    return 

def F_electrostatic(dr, rij):
    """Caclulates the electrostatic culomb focre with cut off
    """
    qi = np.ones((N))*0.42*sc.e
    qi[1::3] = -0.84*sc.e
    qij = np.outer(qi,qi)
    Rc = 2.5*sigma #cut off
    alph = 1/Rc
    sum1 = 2*alph*np.exp(-1*(alph*rij)**2)/np.sqrt(np.pi)/rij
    sum2 = erf(alph*rij)/rij**2
    Fqij =-1*qij/4/np.pi/sc.epsilon_0*(sum1+sum2)
    return Fqij


def F_LJ(rij, Fij):
    """Calculates the 12-6 LJ force
    """
    sigma = 3.21e-10 #meter
    T = Td #Temperature equals target temperature [K]
    epsilon = 0.2*sc.k*T #LJ constants
    rijOO = rij[1::3,1::3] #selects the Oxygen, each third particle.
    Fij[1::3,1::3] = 4*epsilon*(12*sigma**12/rijOO**13-6*sigma**6/rijOO**7) #absolute LJ Force
    return Fij


