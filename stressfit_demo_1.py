# -*- coding: utf-8 -*-
# Created on Thu Apr 15 16:35:29 2021
# @author: Jan Saroun, saroun@ujf.cas.cz
"""Demonstration: strain data fit to macroscopic stress model.

Compliance tensor elements can be loaded from exteralexternal file.
For isotropic materials, compliance can be calculated from given 
Young modulus and Poisson constant.
"""
import numpy as np
import stressfit.commands as comm
import stressfit.shapes as S
import stressfit.sample as sam
import stressfit.mccfit as mc
deg = np.pi/180.


#%% USER INPUT - ENVIRONMENT

# Define environment: input/output folders (can be absolute or relative)
# Set to None for searching in package resources
input_path = None # input data
output_path = './output' # output data
tables_path = None # lookup tables: e.g. material data
sampling_path = None # directory with simulated sampling points 

# Set environment
comm.set_environment(data=input_path, output=output_path, tables=tables_path)

#%% USER INPUT - SAMPLING

# file with the sampling points  
sampling_file = 'events_S_1mm.dat' 
# number of sampling points to load from the file
nev_load = 3000
# number of sampling points to plot
nev_plot = 3000
# number of sampling points to use in concolution
nev_use = 3000

# load sampling
sampling = comm.load_sampling(sampling_file, path=sampling_path, maxn=nev_load)
sampling.print_properties()

#%% USER INPUT - MATERIAL

# Material attenuation, provide either of
# File name: A table with 2 columns: wavelength [A], attenuation [1/cm]
# Float number: attenuation [1/cm]:
att = 'Fe_mu.dat'


## Compliance tensor
# For isotropic material: provide Yound modulus and Poisson ratio
# Othervise: load from database (JSON format)
# See resources/compliance in the package data for examples and file format.)
#
# Compliance tensor is provided as 6 elements, corresponding
# to the 3rd row of the Compliance matrix in Voigt notation in the measurement 
# reference frame, where Q is parallel to z. 
#
# Parameters:
#   resource: resource file (JSON format)
#   phase: phase name as defined in the resource file
#   hkl: indices of corresponding diffraction plane
#   Q_angles: 2 angles which define the Q-vector orientation with respect  
#       to the sample reference frame (ZY Euler angles)

compliance = {'resource':'compliance/Fe_iso', 
              'phase': 'ferrite',
              'hkl': [2,1,1]
              }

# alternatively, define isotropic material as e.g.
# compliance = {'E':220, 'nu': 0.28,'hkl': [2,1,1]}



#%% USER INPUT - SAMPLE SHAPE

# ## Sample definition
# Following block defines the <b>sample shape, position and orientation.</b>
# ### Coordinates
# The laboratory frame is defined by y-axis vertical and z-axis pointing along the incident beam.
# The positioning includes rotation by YXY Euler angles ($\omega$, $\chi$, $\phi$) and a linear shift (multiple commands <code>rotate</code> and <code>moveTo</code> are possible to achieve required position).
# 
# Imported MC events are defined in the laboratory frame, with the origin at the centre of the instrumental gauge volume. Sample position and orientation thus define zero scan position and scan direction in the sample.
# 
# ### Sample shape
# Sample dimensions depend on the selected shape. Choose one of the classes defined in the Shapes folder of the package:
# 
# <code> S.ShapePlate(thickness) </code> <br />
# An infinitely large flat plate of given thickness.
# 
# <code> S.ShapePlateCurved(thickness, length, height, [r1x, r1y],  [r2x, r2y]) </code><br />
# A curved plate of given thickness (z), length (x) and height (y). <code>[r1x, r1y]</code> are curvature radii along x and y of the front surcae (z>0). <code>[r2x, r2y]</code> are the radii for the rear surface.
# 
# <code> S.ShapeCyl(radius, height) </code><br />
# A cylindrical shape with axis along y-axis.
# 
# <code> S.ShapeShellCyl(Rin, Rout, height) </code><br /> 
# A hollow cylinder with axis along y-axis. Rin, Rout are the inner and outer radii.
# 
# <code> S.ShapeSph(radius) </code> <br />
# A spherical sample.
# 
# The next block in this template defines a tube: a hollow cylinder with inner radius 4 mm, outer radius 8 mm and height 50 mm. Zero scan position corresponds to the instrumental gauge volume centered at the surface. Measured strain direction is defined by the angle <code>omega</code>.
# 

# Set sample shape (see comment above).
# Dimensions [mm]
radius1 = 4.
radius2 = 8
height = 50.0
shape = S.ShapeTube(radius1, radius2, height)

#%% USER INPUT - Input data


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
                         scandir=[0., 0., -1.],
                         scanorig=[0, 0, 0],
                         rotctr=[0,0,0], 
                         angles=[225, 0, 0],
                         sampling=sampling)

scans['axi'] = comm.load_input('eps_SS_axi.dat',
                         intensity='int_SS_axi.dat', 
                         scandir=[0., 0., -1.],
                         scanorig=[0, 0, 0],
                         rotctr=[0, 0, 0], 
                         angles=[135, 90, 90],
                         sampling=sampling)

scans['rad'] = comm.load_input('eps_SS_rad.dat',
                         intensity='int_SS_rad.dat', 
                         scandir=[0., 0., -1.],
                         scanorig=[0, 0, 0],
                         rotctr=[0, 0, 0], 
                         angles=[135, 0, 0],
                         sampling=sampling)

# select one data set for plotting 
scan = scans['axi']

#%% Initialization commands



# Set beam attenuation
comm.set_attenuation(att)    

# Set sample shape
comm.set_shape(shape)

# Define compliance
S = comm.set_compliance(**compliance)

# Set scan with associated experiment geometry 
comm.set_scan(scan)   


#%% Plote scene with sampling distribution

# 2D plot of experiment geometry:
# scene width,height in [mm]
scene_range = [16, 16]  
# projection plane (zy=0, xz=1, xy=2)
scene_projection = 1  

# Plot experiment geometry 
# (red arrow shows motion of sampling points in stationary ample)
comm.plot_scene(nev_plot, scan=scan, rang=scene_range, proj=scene_projection)
  

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





