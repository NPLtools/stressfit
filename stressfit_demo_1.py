# -*- coding: utf-8 -*-
# Created on Thu Apr 15 16:35:29 2021
# @author: Jan Saroun, saroun@ujf.cas.cz
"""Demonstration: strain data from macroscopic stress model.

Requires external data for stress-lattice strain conversion, e.g. from
EPSC simulations or in-situ loading test.
"""
import numpy as np
from scipy.interpolate import interp1d
import stressfit.shapes as S
import stressfit.sample as sam
import stressfit.mccfit as mc
import stressfit.graphs as gr
import stressfit.commands as comm
from stressfit.compliance import Compliance
import stressfit.dataio as dataio
deg = np.pi/180.

# Define environment
# input data folder
inpath = mc.path2win(r'.\input')
# output data folder
outpath = mc.path2win(r'.\output')
# path to sampling points data
gpath = mc.path2win(r'./input/')

#%% Sampling distribution

# Load Monte Carlo events representing the sampling distribution from a text file.
# The event table contains neutron coordinates, weights and dhkl values.  You need to specify column numbers for position, ki and kf vectors, weights and dhkl.
# 
# Imported MC events are defined in the laboratory frame, with the origin at the centre of the instrumental gauge volume. Sample and orientation thus defines zero scan position and scan direction in the sample.

# load sampling points
gpath = r'./input/'
sampling = comm.load_sampling(gpath + "events_S_1mm.dat", verbose=True)

#%% Define compliance matrix

# create Compliance from resources
Smat = Compliance.from_resource('compliance/Fe_iso', phase='ferrite')
# define phase and hkl 
phase = 'ferrite'
hkl = [2,1,1]


#%% Define beam attenuation

# There are two options, uncomment one of them: 

# Option 1: Load attenuation as a table (wavelength, mu), [1/cm]. 
exttab = dataio.load_data('Fe_mu.dat', kind='tables')
sam.setExtinction(table=exttab)

# Option 2: Set attenuation as a single coefficient [1/cm]:
# sam.setExtinction(mu=1.96) 


#%% Input data

# Define sample shape
#---------------------

# Set sample shape
# Dimensions [mm]
radius1 = 4.
radius2 = 8
height = 50.0
shape = S.ShapeShellCyl(radius1, radius2, height)


# Define experimental data
#--------------------------

# For each measured reflection and sample orientation, provide
# - strain data file: sample position (encoder value), strain, error
# - intensity data file: sample position (encoder value), intensity, error
# - scan direction and origin (encoder = 0) in sample coordinates
# - sample rotation centre (sample coordinates)
# - sample orientation (Euler angles)

    
# Load measured strains and integral peak intensities and encapsulate 
# them with corresponding metadata.

scans = {}
scans['hoop'] = comm.load_input('eps_SS_hoop.dat', 
                         intensity='int_SS_hoop.dat', 
                         path = inpath,
                         scandir=[0., 0., -1.],
                         scanorig=[0, 0, 0],
                         rotctr=[0,0,0], 
                         angles=[225, 0, 0],
                         sampling=sampling)

scans['axi'] = comm.load_input('eps_SS_axi.dat',
                         intensity='int_SS_axi.dat', 
                         path = inpath,
                         scandir=[0., 0., -1.],
                         scanorig=[0, 0, 0],
                         rotctr=[0, 0, 0], 
                         angles=[135, 90, 90],
                         sampling=sampling)

scans['rad'] = comm.load_input('eps_SS_rad.dat',
                         intensity='int_SS_rad.dat', 
                         path = inpath,
                         scandir=[0., 0., -1.],
                         scanorig=[0, 0, 0],
                         rotctr=[0, 0, 0], 
                         angles=[135, 0, 0],
                         sampling=sampling)


# Select which data to process
data = scans['hoop']


#%% Initiate dependent parameters

# Sample rotation. 
# At zero angles, the z-axis is parallel to the beam, y-axis is vertical.
# omega, chi, phi are Euler angles (rotation around Y, X, Y).
# Example for a plate symmetric reflection with incidence angle = theta:
#   set chi=0, phi=0, omega=90 + theta

# Scattering angle
take_off = sampling['tth']*deg 

# Set sample position
angles = data['angles']*deg
shape.rotate(*angles)
shape.moveTo(data['scanorig'])

# Assign shape
sam.shape = shape



# define ki and kf vectors in laboratory frame
ki = np.array([0., 0., 1.])  # default for SIMRES simulations
kf = sam.rotate(ki, 1, take_off) # rotation of ki by the take-off angle

# Set the events to the sample component
sam.setSampling(data['sampling'])

# set orientation of Q to Smat
ki0 = shape.getLocalDir(ki)
kf0 = shape.getLocalDir(kf)
Q = ki0-kf0
Smat.set_hkl(phase, hkl, Q=Q)


#%% Plot sampling scene

# number of sampling events to show
nev = 2000
# plot range in [mm]
rang = [13, 13]  
# projection plane (zy=0, xz=1, xy=2)
proj = 1  
outpng = outpath + 'scene.png'
gr.plotScene(rang, proj, shape, ki, kf, data['scandir'], sam.getSampling(nev), 
             save = True, file = outpng)

   

#%% Define stress model
# Give initial values, followed by flags (0|1), fx=1 means a free variable, 0 for fixed.  

# depth values [mm]
x =  [0., 1., 1.5, 2.0 , 2.5, 3.0 ,  3.5  , 3.7  ,   4.]
# stress values [MPa] for each principal direction
sig_xx =  [0., 0., 10., 20., 53. , 80., -200. , -395., -250.]
sig_yy =  sig_xx
sig_zz = len(x)*[0.]
sig = np.array([sig_xx, sig_yy, sig_zz]).T
fx = len(x)*[1]
fy = len(x)*[1]

# Set the method for interpolation between the nodes.
# Use one of 'natural','clamped', 'PCHIP', 'Akima'
# PCHIP = Piecewise cubic Hermitian interpolation polynomial
interpolation = 'natural'


#%% Convert stress to strain

eps = Smat.strain(sig)


#%%

# Initialize model
sfit = mc.Sfit(nev=3000, xdir=data['scandir'], ftol=1.e-3)

# data to be fitted
sfit.data = data['eps']

# choose randomly a subset of sampling events
sam.shuffleEvents()

# define strain distribution, dim=number of points for interpolation 
# par = nodes [x,y] values
# vary = corresponding flags for fixed (0) or free(1) variables.
sfit.defDistribution(par=[x, eps[:,0]], vary=[fx, fy], ndim=100, scaled=True)

# define function scaling (amplitude, strain offset, depth-shift) 
sfit.defScaling(par=[1., 0, 0], vary=[0, 0, 0])

# define interpolation method
sfit.setInterpModel(interpolation)

#mc.runFit(sfit, maxiter=1, areg=1e-7, bootstrap=False, loops=False, guess=False)
sfit.reportFit()
#gr.plotStrainModel(sfit)





