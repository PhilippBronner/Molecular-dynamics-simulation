#Data analysis
"""
Created on Fri Jun 18 23:19:02 2021

@author: philipp
"""
from numba import njit
from numba import jit
import numpy as np
import analysis as ana
from initialize import *
from tqdm import tqdm #step counter
import scipy.constants as sck
import scipy.stats as sc
import numpy.linalg as npl



##Import data
name0 = 'trajectory_1P.dump'

name1 = 'trajectory2.dump' #particle density 0.05 sigma**3
name2 = 'trajectory1.dump' #particle density 0.25 sigma**3
r_1P = ana.read_traj(name0, 1000, 3)
# bond length and angle
rOH1 = r_1P[:,0,:]-r_1P[:,1,:]
rOH2 = r_1P[:,2,:]-r_1P[:,1,:]
theta = np.zeros((1000))
for ii in range(1000):
    theta[ii] = np.arccos(rOH1[:,ii].dot(rOH2[:,ii])/(npl.norm(rOH1[:,ii])*npl.norm(rOH2[:,ii])))
theta = theta*360/2/np.pi

plt.hlines(np.mean(npl.norm(rOH1[:,1:300],axis = 0)), 0, 300,color = 'r', label = 'r OH1 mean', linestyle = '--')
plt.hlines(np.mean(npl.norm(rOH2[:,1:300],axis = 0)), 0, 300,color = 'b', label = 'r OH2 mean', linestyle = '--')
plt.plot(npl.norm(rOH1[:,1:300],axis = 0), label = 'r OH1', color = 'r')
plt.plot(npl.norm(rOH2[:,1:300],axis = 0), label = 'r OH2', color = 'b')
plt.xlabel('t [0.5fs]')
plt.ylabel('rOH [m]')
plt.title('distance r OH')
plt.legend()
plt.savefig('rOH.jpg', dpi = 300)
plt.show()

plt.title('Development of theta')
plt.xlabel('t [0.5fs]')
plt.plot(theta[1:1000])
plt.ylabel('Theta [°]')
plt.savefig('Theta.jpg',dpi = 300)
plt.show()


#partc)

name1 = 'trajectory_6P.dump' #particle density 0.25 sigma**3
N = int(6**3*3)
r_4P = ana.read_traj(name1, 5000, N)
M = int(N/3)

#center of mass
#ri_mean = np.mean(r_4P, axis = 2)
#masses = initial_masses(N, M, mO, mH)
Mges = np.sum(masses)
r_4PO = r_4P[:,1::3,:]
r_4P1 =np.where(abs(r_4PO) >= 5*L, L, r_4PO)
rsp = np.sum(r_4P1, axis = 1)/len(r_4P1[0,:,0])

#rsp = np.tensordot(r_4P1,masses, axes = [1,0])/Mges

rsp = rsp[:,np.newaxis,:]
rfm = r_4PO - rsp

rfm = npl.norm(rfm, axis = 0)
Ndr = 50
n, edge = np.histogram(rfm, bins = Ndr, range = (0,L))
dens = n/5000 /(4*np.pi*edge[1:]**2*(edge[1]-edge[0]))
plt.plot(edge[2:],dens[1:])
plt.title('particle density inside the water droplet')
plt.ylabel('density [1/cm³]')
plt.xlabel(r'Distance $R_{cm}$  $|r-r_cm|$ [m]')
plt.savefig('density.jpg', dpi = 300)
plt.show()

#r_traj1 = r_traj[:,:,int(2000/5):]
#r_traj = ana.read_traj(name2)
#r_traj2 = r_traj[:,:,int(2000/5):]

