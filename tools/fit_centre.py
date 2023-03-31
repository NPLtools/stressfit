
# coding: utf-8

# In[192]:


import numpy as np
from matplotlib import pyplot as plt
import stressfit.shapes as S
import stressfit.sample as sam
import stressfit.graphs as gr
from stressfit.mccfit import Ifit, params2array, path2win
from stressfit.geometry import Geometry
from lmfit import Minimizer, Parameters

# Fit gauge centre from intensity scan data

# In[228]:
# input data folder
inpath = path2win(r'D:\Data\HK4-2018\April\stressfit\input')
# output data folder
outpath = path2win(r'D:\Data\HK4-2018\April\out\stressfit\output')


# Orientation:
deg = np.pi/180.
omega = -90*deg + 32.5*deg
chi = 0*deg
phi = 0.0*deg
# Dimensions [mm]
radius1 = 4.
radius2 = 8.
thickness = 5.
length = 100.0
height = 50.0
# Curvatures [1/mm]
rho1 = 0.e-2*np.array([1., 1.])
rho2 = 0.e-2*np.array([1., 1.])

# Sample shape - choose one of classes defined in the samples folder
#shape = S.ShapePlateCurved(thickness, length, height, rho1, rho2)
#shape = S.ShapeCyl(radius2, height)
shape = S.ShapeShellCyl(radius1, radius2, height)

# Set sample position
shape.rotate(omega, chi, phi)
shape.moveTo(np.array([0., 0., 0.]))
#shape.moveTo(np.array([-0.294, 0., 0.187]))

# Assign shape and attenuation to the sample
sam.shape = shape

# Set attenuation as a single coefficient or a table (wavelength, mu), [1/cm]
exttab = np.loadtxt('Fe_mu.dat')
sam.setExtinction(table=exttab)
#plt.plot(exttab[:,0],exttab[:,1],'k-')
#plt.show()
#sam.setExtinction(mu=1.15)

# define ki and kf vectors in laboratory frame
take_off = 62.555*deg  # scattering angle
ki = np.array([0., 0., 1.])  # default for SIMRES simulations
kf = Geometry.rotate(ki, 1, take_off) 


# ## Sampling distribution
# Load Monte Carlo events representing the sampling distribution from a text file.
# The event table contains neutron coordinates, weights and dhkl values.  You need to specifie column numbers for position, ki and kf vectors, weights and dhkl.

# In[204]:

# load sampling points
#data = np.loadtxt(inpath + "gauge_R_0158.dat")
data = np.loadtxt(inpath + "gaugeR180.dat")
#data = np.loadtxt("C:\\Saroun\\Simulace\\NPI Rez\HK4\\gauge1.dat")
columns = [1, 4, 7, 10, 11]  # index of r[0], ki[0], kf[0], weight and dhkl (column indexing starts from 0)
sumP = np.sum(data[:,columns[3]])
x0 = data[:,columns[0]].dot(data[:,columns[3]])/sumP
y0 = data[:,columns[0]+1].dot(data[:,columns[3]])/sumP
z0 = data[:,columns[0]+2].dot(data[:,columns[3]])/sumP
print('gauge centre: [{:g}, {:g}, {:g}] '.format(x0, y0, z0))
sam.setSamplingEvents(data, columns, ctr = [x0, y0, z0])


# ## Scan definition
# Define scan steps, direction and calculate corresponding depth scale. At the end, the diffraction geometry is plotted together with the sampling points. Color scale of the sampling points shows the associated pseudo-strains.

# In[229]:

nev = 2000  # number of sampling events to use for convolution

sdir1 = [0., 0., -1.] # scan direction 1 (local coordinates)
sdir2 = [1., 0., 0.]  # scan direction 2 (local coordinates)

# plot the situation
nd = nev  # number of events to show
rang = [12, 12]  # plot range in [mm]
proj = 1  # projection axis
r = data[0:nd, 1:4]
p = data[0:nd, 10]
dhkl = data[0:nd, 11]
gr.plotScene(rang, proj, shape, ki, kf, sdir1, r  - sam.sctr, p, dhkl)
gr.plotScene(rang, proj, shape, ki, kf, sdir2, r  - sam.sctr, p, dhkl)

# ## Input data
# Load integral peak intensities and strains measured as a function of scan depth in two separate arrays with three columns: depth [mm], intensity [any unit] or strain [$\mu\epsilon = 10^{-6}$] and error (std. deviation).

# In[230]:

# 1st direction, intensities, file name
#intfile1 = ['tubeX1_int.dat', 'tubeX2_int.dat']
intfile1 = ['tube_rad_adj.dat']
ctr1 = 30.76 # nominal scan centre

# 2nd direction, intensities, file name
#intfile2 = ['tubeY1_int.dat','tubeY2_int.dat']
intfile2 = ['tube_hoop_adj.dat']
ctr2 = 24.96 # nominal scan centre

# load data
n1 = len(intfile1)
intdata1 = None
for i in range(n1):
    fn = intfile1[i]
    dt = np.loadtxt(inpath + fn)
    if (i>0):
        intdata1 = np.concatenate((intdata1,dt),axis=0)
    else:
        intdata1 = dt

n2 = len(intfile1)
intdata2 = None
for i in range(n2):
    fn = intfile2[i]
    dt = np.loadtxt(inpath + fn)
    if (i>0):
        intdata2 = np.concatenate((intdata2,dt),axis=0)
    else:
        intdata2 = dt

# subtract centre estimate

        
intdata1[:,0] = intdata1[:,0] - ctr1
intdata2[:,0] = intdata2[:,0] - ctr2

# %%get scan points for intensities and strains
# They are normally the same, but they are allowed to be different
# intensities
nci1 = intdata1.shape[0]  # number of scan points
xci1 = intdata1[:,0]  # steps
xcsi1 = shape.scanDepth(xci1, sdir1)  # depth scale

nci2 = intdata2.shape[0]  # number of scan points
xci2 = intdata2[:,0]  # steps
xcsi2 = shape.scanDepth(xci2, sdir2)  # depth scale



# In[231]:

# MODEL PARAMETERS
# Give initial values, followed by flags (0|1) for free variables 

# depth values [mm]
x =  [0., 1., 2., 3., 4., 5.]
fx = [0, 0, 0, 0, 0, 0]

# intensity values [rel. units]
y1 = [1., 1., 1., 1., 1.05, 1.1]
y2 = [1., 1., 1., 0.95, 0.85, 0.82]
fy = [0, 0, 0, 0, 0, 0]

# Background and amplitude
A1 = 50
A2 = 40
B = 0
fA = 1
fB = 0

# zero encoder position 
zc = 0.0
fzc = 0

# set spline order (1 for linear, 3 for cubic splines)
iord = 2
#--------------------------------------------

# Initialize model
ifit1 = Ifit(nev=5000, xdir = sdir1, ftol=1.e-3)
ifit2 = Ifit(nev=5000, xdir = sdir2, ftol=1.e-3)
ifit1.defDistribution([x, y1], [fx, fy], iord=iord)
ifit1.defScaling([A1, B, zc], [0, 0, 0], minval=[0., 0., -np.inf])
ifit2.defDistribution([x, y2], [fx, fy], iord=iord)
ifit2.defScaling([A2, B, zc], [0, 0, 0], minval=[0., 0., -np.inf])

# make convolution and plot
[yf1, ef1, pos] = ifit1.getSmearedFnc(xci1)
[yf2, ef2, pos] = ifit2.getSmearedFnc(xci2)
plt.title('Measured vs. model intensities ')
plt.xlabel('Scan position, mm')
plt.ylabel('Intensity')
ymin=np.min(np.concatenate((yf1, yf2,intdata1[:,1], intdata2[:,1]), axis=0))
ymax=np.max(np.concatenate((yf1, yf2,intdata1[:,1], intdata2[:,1]), axis=0))
dy = ymax - ymin
plt.ylim(ymin, ymax+0.05*dy)
plt.errorbar(xci1, yf1, yerr=ef1, fmt='k-', label='model guess')
plt.errorbar(xci2, yf2, yerr=ef2, fmt='b-', label='model guess')
plt.errorbar(xci1, intdata1[:,1], yerr=intdata1[:,2], fmt='ko', label='data')
plt.errorbar(xci2, intdata2[:,1], yerr=intdata2[:,2], fmt='bo', label='data')
plt.legend(loc='best', frameon=False)
plt.show()


# ### Run fit
# Is the above estimate good? Then execute the following box to run fitting procedure and plot results.

# In[201]:

def costFnc(params, scan1, scan2):
    """ Function to be passed to the minimizer.
    """
    xc = params['xc'].value
    yc = params['yc'].value
    zc = params['zc'].value
    A1 = params['A1'].value
    A2 = params['A2'].value
    shape.moveToAbs(np.array([xc, yc, zc]))
    scan1.params.get('A').value=A1
    scan2.params.get('A').value=A2
    
    [yf1, ef1, pos1] = scan1.getSmearedFnc(scan1.data[:,0])
    [yf2, ef2, pos2] = scan2.getSmearedFnc(scan2.data[:,0])
    er1 = np.sqrt(ef1**2 + scan1.data[:,2]**2)
    er2 = np.sqrt(ef2**2 + scan2.data[:,2]**2)
    
    yval = np.concatenate((scan1.data[:,1], scan2.data[:,1]), axis=0)
    yfit = np.concatenate((yf1, yf2), axis=0)
    err = np.concatenate((er1, er2), axis=0)
    
    return (yfit - yval)/err

def callb(params, iter, resid, *fcn_args, **fcn_kws):
    """ Callback to provide fit progress info.
    """
    global chi
    ivar = 0
    s = ''
    for key in params:
        p = params[key]
        ivar += p.vary
        if (p.vary):
            s += ' {:g}'.format(p.value)
    sumres = np.sum(resid**2)/(resid.size - ivar)
    if (iter>1):
        chi0=chi
    else:
        chi0=1e+10    
    if (sumres < chi0):
        chi = sumres
        print('iter='+str(iter), 'chi={:g}'.format(chi), 'par='+s)

# assign fit data
ifit1.data = intdata1
ifit2.data = intdata2

# initial position
r0 = [0., 0., 0.]
# initialize parameters
par = Parameters()
par.add('xc', value=r0[0], vary=1, min=-2, max=1)
par.add('yc', value=r0[1], vary=0)
par.add('zc', value=r0[2], vary=1, min=-1.0, max=1.0)
par.add('A1', value=A1, vary=1, min=0.0, max=3*A1)
par.add('A2', value=A2, vary=1, min=0.0, max=3*A2)
maxiter=100
# run fit
minner = Minimizer(costFnc, par, fcn_args=(ifit1, ifit2), iter_cb=callb)
result = minner.minimize(method='leastsq', ftol=1.0e-6, epsfcn=0.1, maxfev=maxiter)
par1 = result.params
rfit = params2array(par1)[:3,:]
fmt = 'Fitted centre correction:\nx = {:g} +- {:g}\nz = {:g} +- {:g}\n'
print(fmt.format(rfit[0, 0], rfit[0, 1], rfit[2, 0], rfit[2, 1]))
# move to the new position
shape.moveToAbs(rfit[:3,0])


# %% Plot result
# make convolution and plot
[yf1, ef1, pos] = ifit1.getSmearedFnc(xci1)
[yf2, ef2, pos] = ifit2.getSmearedFnc(xci2)
plt.title('Measured vs. model intensities ')
plt.xlabel('Scan position, mm')
plt.ylabel('Intensity')
ymin=np.min(np.concatenate((yf1, yf2,intdata1[:,1], intdata2[:,1]), axis=0))
ymax=np.max(np.concatenate((yf1, yf2,intdata1[:,1], intdata2[:,1]), axis=0))
dy = ymax - ymin
plt.ylim(ymin, ymax+0.05*dy)
plt.errorbar(xci1, yf1, yerr=ef1, fmt='k-', label='model guess')
plt.errorbar(xci2, yf2, yerr=ef2, fmt='b-', label='model guess')
plt.errorbar(xci1, intdata1[:,1], yerr=intdata1[:,2], fmt='ko', label='data')
plt.errorbar(xci2, intdata2[:,1], yerr=intdata2[:,2], fmt='bo', label='data')
plt.legend(loc='best', frameon=False)
plt.show()



# %% centred  geometry
gr.plotScene(rang, proj, shape, ki, kf, sdir1, r  - sam.sctr, p, dhkl)
gr.plotScene(rang, proj, shape, ki, kf, sdir2, r  - sam.sctr, p, dhkl)

pos = rfit[:3,0]*[1, 1, -1] + [ctr2, 0, ctr1]
err = rfit[:3,1]
fmt = 'Sample centre:\nx = {:g} +- {:g}\nz = {:g} +- {:g}\n'
print(fmt.format(pos[0], err[0], pos[2], err[2]))
