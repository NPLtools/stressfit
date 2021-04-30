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
deg = np.pi/180.

# Define environment
# input data folder
inpath = mc.path2win(r'.\input')
# output data folder
outpath = mc.path2win(r'.\output')



#%% 1. Sample orientation:

# Sample rotation. 
# At zero angles, the z-axis is parallel to the beam, y-axis is vertical.
# omega, chi, phi are Euler angles (rotation around Y, X, Y).
# Example for a plate symmetric reflection with incidence angle = theta:
#   set chi=0, phi=0, omega=90 + theta

omega = 180*deg + 45*deg
chi = 0.0*deg
phi = 0.0*deg

# Scattering angle
take_off = 91.77*deg 

# Set sample shape (see comment above).
# Dimensions [mm]
radius1 = 4.
radius2 = 8
height = 50.0
shape = S.ShapeShellCyl(radius1, radius2, height)

# Set sample position
shape.rotate(omega, chi, phi)
shape.moveTo(np.array([0., 0., 0]))

# Assign shape
sam.shape = shape

# Define beam attenuation coefficient. Uncomment one of the two options: 
# Option 1: Set attenuation as a table (wavelength, mu), [1/cm]. The file must be in the input directory.
exttab = np.loadtxt('tables/Fe_mu.dat')
sam.setExtinction(table=exttab)  # lookup table

# Option 2: Set attenuation as a single coefficient [1/cm]:
# sam.setExtinction(mu=1.96) # single value

# define ki and kf vectors in laboratory frame
ki = np.array([0., 0., 1.])  # default for SIMRES simulations
kf = sam.rotate(ki, 1, take_off) # rotation of ki by the take-off angle

#%% Sampling distribution
# Load Monte Carlo events representing the sampling distribution from a text file.
# The event table contains neutron coordinates, weights and dhkl values.  You need to specify column numbers for position, ki and kf vectors, weights and dhkl.
# 
# Imported MC events are defined in the laboratory frame, with the origin at the centre of the instrumental gauge volume. Sample and orientation thus defines zero scan position and scan direction in the sample.

# load sampling points
gpath = r'./input/'

def load_sampling(filename, columns=[1, 4, 7, 10, 11], maxn=0, verbose=False):
    """Load Monte Carlo events representing the sampling distribution.
    
    The event table contains neutron coordinates, weights and dhkl values.  
    You need to specify column numbers for position, ki and kf vectors, 
    weights and dhkl.
    
    Imported MC events are defined in the laboratory frame, with the origin 
    at the centre of the instrumental gauge volume. Sample and orientation 
    thus defines zero scan position and scan direction in the sample.
    

    Parameters
    ----------
    filename : str
        Input file name.
    columns : list, optional
        Column indices of r[0], ki[0], kf[0], weight and dhkl (starts from 0)
    maxn : int, optional
        Maximum number or records. If zero, tak all records in the file.
    verbose: boolean
        If true, print calculated gauge parameters.

    Returns
    -------
    dict
        - data: event data
        - columns: given input parameter
        - nrec: number of records
        - ctr: sampling centre
        - dmean: mean dhkl
    """
    out = {}
    data = np.loadtxt(filename)
    if maxn:
        nrec = min(maxn, data.shape[0])
    else:
        nrec = data.shape[0]
    # Calculate centre of mass of the distribution 
    P = data[:nrec,columns[3]]/np.sum(data[:,columns[3]])
    ctr = np.zeros(3)
    ki = np.zeros(3)
    kf = np.zeros(3)
    for i in range(3):
        ctr[i] = data[:nrec, columns[0] + i].dot(P)
        ki[i] = data[:nrec, columns[1] + i].dot(P)
        kf[i] = data[:nrec, columns[2] + i].dot(P)
    
    wav = 2*np.pi/np.sqrt(ki.dot(ki))
    dmean = data[:nrec, columns[4]].dot(P)
    tth = np.arcsin(wav/2/dmean)*360/np.pi
    if verbose:
        print('Loaded event list with {:d} records'.format(nrec))    
        print('Gauge centre: [{:g}, {:g}, {:g}] '.format(*ctr))
        print('Mean wavelength: {:g}'.format(wav))
        print('2 theta: {:g}'.format(tth))
        print('d0: {:g}\n'.format(dmean))    
    out['data'] = data[:nrec,:]
    out['columns'] = columns
    out['nrec'] = nrec
    out['ctr'] = ctr
    out['dmean'] = dmean
    out['wav'] = wav
    out['tth'] = tth
    return out

sampling = load_sampling(gpath + "events_S_1mm.dat", verbose=True)

# Set the events to the sample component
sam.setSampling(sampling)


#%% Scan definition
# Define scan steps, direction and calculate corresponding depth scale. At the end, the diffraction geometry is plotted together with the sampling points. <b>The red arrow shows the direction of sample motion.</b>  Color scale of the sampling points shows the associated pseudo-strains.

# Scan direction (local coordinates)
sdir = [0., 0., -1.]
# number of sampling events to use for convolution
nev = 2000    

# plot the situation
nd = nev  # number of events to show
rang = [13, 13]  # plot range in [mm]
proj = 1  # projection plane (zy=0, xz=1, xy=2)
outpng = outpath + 'scene.png'
gr.plotScene(rang, proj, shape, ki, kf, sdir, sam.getSampling(nev), 
             save = True, file = outpng)



#%% Input data

# For each measured reflection and sample orientation, provide
# - strain data file: sample position (encoder value), strain, error
# - intensity data file: sample position (encoder value), intensity, error
# - scan direction and origin (encoder = 0) in sample coordinates
# - sample rotation centre (sample coordinates)
# - sample orientation (Euler angles)


def load_input(epsfile, intfile, 
              scandir=[0, 0, 1], 
              scanorig=[0, 0, 0], 
              rotctr=[0, 0, 0],
              angles=[180+45, 0, 0],
              sampling=None):
    """Load experimental data and metadata.
    
    Parameters
    ----------
    epsfile : str
        File name for strain data: position (encoder value), strain, error
    intfile : str, optional
        File name for intensity data: position (encoder value), strain, erro.
    scandir : list or array, optional
        Scan direction in sample coordinates
    scanorig : list or array, optional
        Sscan origin (encoder = 0) in sample coordinates
    rotctr : list or array, optional
        Sample rotation centre (sample coordinates)
    angles : list or array, optional
        Sample orientation (Euler angles YXY)
    sampling: dict
        Sampling events loaded by the function load_gauge.

    Returns
    -------
    dict
        Input data: keys are derived from argument names.

    """
    out = {}
    out['eps'] = np.loadtxt(inpath + epsfile)
    if intfile:
        out['int'] = np.loadtxt(inpath + intfile)
    else:
        out['int'] = None
    out['scandir'] = np.array(scandir)
    out['scanorig'] = np.array(scanorig)
    out['rotctr'] = np.array(rotctr)
    out['angles'] = np.array(angles)
    out['sampling'] = sampling
    return out
    
# Load measured strains and integral peak intensities and encapsulate 
# them with corresponding metadata.

scans = {}
scans['hoop'] = load_input('eps_SS_hoop.dat', 'int_SS_hoop.dat', 
                         scandir=sdir,
                         scanorig=[0, 0, 0],
                         rotctr=[0,0,0], 
                         angles=[225, 0, 0],
                         sampling=sampling)

scans['axi'] = load_input('eps_SS_axi.dat', 'int_SS_axi.dat', 
                         scandir=sdir,
                         scanorig=[0, 0, 0],
                         rotctr=[0, 0, 0], 
                         angles=[135, 90, 90],
                         sampling=sampling)

scans['rad'] = load_input('eps_SS_rad.dat', 'int_SS_rad.dat', 
                         scandir=sdir,
                         scanorig=[0, 0, 0],
                         rotctr=[0, 0, 0], 
                         angles=[135, 0, 0],
                         sampling=sampling)

#%% Load data for stress - lat. strain conversion

class LatticeResponse():
    def __init__(self, sigmaL, sigmaT):
        self.sigmaL = sigmaL
        self.sigmaT = sigmaT
        intpL = interp1d(sigmaL[:,1], sigmaL[:,0], kind='linear')
        intpT = interp1d(sigmaT[:,1], sigmaT[:,0], kind='linear')
        sigmin = max(min(sigmaL[:,1]), min(sigmaT[:,1]))
        sigmax = min(max(sigmaL[:,1]), max(sigmaT[:,1]))
        sig = np.linspace(sigmin, sigmax, num=50)
        epsL = intpL(sig)
        epsT = intpT(sig)
        self.intpEL = interp1d(sig, epsL, kind='linear')
        self.intpET = interp1d(sig, epsT, kind='linear')
    
        
    def strain(self, sigma):
        s = np.reshape(sigma, (-1, 3))
        EL = self.intpEL(s)
        ET = self.intpET(s)
        exx = EL[:,0] + ET[:,1] + ET[:,2]
        eyy = EL[:,1] + ET[:,0] + ET[:,2] 
        ezz = EL[:,2] + ET[:,0] + ET[:,1]
        return np.array([exx,eyy,ezz]).T
    

# create artificial response, just elastic ...
E = 220e-3
nu = 0.28
eps = np.linspace(-6000, 6000, num=100)
sig = np.linspace(-450, 450, num=100)
epsL = sig/E
epsT = -nu*sig/E
sigepsL = np.array([epsL, sig]).T
sigepsT = np.array([epsT, sig]).T

# create response function
lat211 = LatticeResponse(sigepsL,sigepsT)

print(lat211.strain([0,0,300]))


sdata = np.loadtxt(inpath+'stress_table.dat')
z = sdata[:,0]
sigma = sdata[:,1:4]
epsilon = lat211.strain(sigma)
np.savetxt(outpath+'tmp_epsilon.dat',epsilon, delimiter='\t', fmt='%g')



#%% Define stress model
# Give initial values, followed by flags (0|1), fx=1 means a free variable, 0 for fixed.  

# depth values [mm]
x =  [0., 1., 2.0 , 2.5, 3.0 ,  3.5  , 3.7  ,   4.]
# stress values [MPa] for each principal direction
sig_xx =  [0., 0., 10., 20., 53. , 80., -200. , -395., -250.]
sig_yy =  sig_xx
sig_zz = len(x)*[0.]
sig = np.array([sig_xx, sig_yy, sig_zz]).T

# Set the method for interpolation between the nodes.
# Use one of 'natural','clamped', 'PCHIP', 'Akima'
# PCHIP = Piecewise cubic Hermitian interpolation polynomial
interpolation = 'natural' 

#%% Calculate smeared strain data



fx = len(x)*[1]
fx[-1] = 0

# initial strain values [1e-6]
y =  len(x)*[0.]
fy = len(y)*[1]
fy[0]=0
