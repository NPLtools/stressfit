# -*- coding: utf-8 -*-
"""
Conversion strain to stress
Created on Tue May 15 19:44:07 2018

@author: Jan Saroun, saroun@ujf.cas.cz
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import splrep, splev

def path2win(name):
    out = name.replace('\\','/')
    if (out[-1] != '/'):
        out += '/'
    return out

def plotResult(d, sig, esig, title='', fname='', ytitle='stress, MPa'):
# plot result
    file = outpath+pfx+fname+sfx
    plt.xlabel('Information depth, mm')
    plt.ylabel(ytitle)
    plt.title(title)
    plt.grid()
    n=nx
    leg = ['hoop','axial','radial']
    fm = ['k-','b-','r-']
    for i in range(3):
        plt.errorbar(d[:n], sig[:n,i], yerr=esig[:n,i], fmt=fm[i], label=leg[i])
    plt.legend(loc='best', frameon=False)
    plt.savefig(file+'.png', bbox_inches='tight')
    plt.show()
    
    # save result
    hdr = "# Reconstructed profiles, {} , "+3*"{} ".format(ytitle, *leg)+"\n"
    hdr += "# x"
    fmt = "\ty{:d}\te_y{:d}"
    for i in range(3):
        hdr += fmt.format(i,i)
    out = np.array([d[:],sig[:,0],esig[:,0],sig[:,1],esig[:,1],sig[:,2],esig[:,2]]).T
    print("Saving {}".format(file+".dat"))
    np.savetxt(file+'.dat',out,fmt='%g', delimiter='\t', header=hdr)
    

def ForcePlate(d, sig, C):
    """Calculate strains to calncel force integrals.
    For planparallel plate, assumes d spanning full thickness and equidistant d[i] points.
    """
    dd = d[1]-d[0]
    F = np.sum(sig*dd, axis=0)    
    return F/(C.dot(np.ones(3))*(max(d)-min(d)))

# %% INPUT
    
inpath= path2win(r'D:\Saroun\Publikace\2019\ECNS\py\output')  # input folder
outpath=path2win(r'D:\Saroun\Publikace\2019\ECNS\py')  # output folder
pfx = 'B_'

# Define stiffness matrix
EY = 220*1e-3;
nu = 0.28;
M = [[1-nu, nu, nu],[nu, 1-nu, nu],[nu, nu, 1-nu]]
C = EY/(1+nu)/(1-2*nu)*np.array(M)

# Load data for 3 components e_xx, e_yy, e_zz
exx = np.loadtxt(inpath+'eps_'+pfx+'hoop'+'_model.dat', usecols = (0,1,2))
eyy = np.loadtxt(inpath+'eps_'+pfx+'axi'+'_model.dat', usecols = (0,1,2))
ezz = np.loadtxt(inpath+'eps_'+pfx+'rad'+'_model.dat', usecols = (0,1,2))

# Subtract eps0 to reach sig_zz=0 ?
subeps0 = False

# %% Calculate data on common interpolated information depth scale 


# Data interpolation
txx = splrep(exx[:,0], exx[:,1], s=0, k=1)
tyy = splrep(eyy[:,0], eyy[:,1], s=0, k=1)
tzz = splrep(ezz[:,0], ezz[:,1], s=0, k=1)
etxx = splrep(exx[:,0], exx[:,2], s=0, k=1)
etyy = splrep(eyy[:,0], eyy[:,2], s=0, k=1)
etzz = splrep(ezz[:,0], ezz[:,2], s=0, k=1)
xmin = max(exx[0,0],eyy[0,0],ezz[0,0])
xmax = min(exx[-1,0],eyy[-1,0],ezz[-1,0])
nx = min(exx.shape[0], eyy.shape[0], ezz.shape[0])
d = np.linspace(xmin, xmax,num=nx)

# interpolate strains
eps = np.zeros((nx,3))
eeps = np.zeros((nx,3))
eps[:,0] = splev(d, txx, ext=1)
eps[:,1] = splev(d, tyy, ext=1)
eps[:,2] = splev(d, tzz, ext=1)
eeps[:,0] = splev(d, etxx, ext=1)
eeps[:,1] = splev(d, etyy, ext=1)
eeps[:,2] = splev(d, etzz, ext=1)

# calculate eps0 from the condition sig_zz(0) = 0
eps0 = (nu*eps[0,0] + nu*eps[0,1] + (1-nu)*eps[0,2])/(1+nu)
eeps0 = nu**2*eeps[0,0]**2 + nu**2*eeps[0,1]**2 + (1-nu)**2*eps[0,2]**2
eeps0 = np.sqrt(np.absolute(eeps0))/(1+nu)
print('\nCondition sig_zz(0) = 0 requires eps_0 = {:g} +- {:g}'.format(eps0,eeps0))

# Subtract eps0?
if subeps0:
    sfx = '-eps0'
    ie0 = 1
else:
    sfx = ''
    ie0 = 0

# calculate stress
sig = np.zeros((nx,3))
esig = np.zeros((nx,3))
for i in range(nx):
    sig[i,:] = C.dot(eps[i,:]-eps0*ie0)
    esig[i,:] = np.sqrt(np.absolute((C*C).dot(eeps[i,:]**2)))


plotResult(d, eps, eeps, title='Strain distribution', fname='eps_final', ytitle='Strain, 1e-6')
plotResult(d, sig, esig, title='Stress distribution', fname='sigma_final')

# %% Apply force equilibrium
"""
EPS0 =  ForcePlate(d, sig, C)
# Correct for EPS0:
sig1 = np.zeros((nx,3))
eps1 = np.zeros((nx,3))
for i in range(nx):
    eps1[i,:] = eps[i,:]-EPS0
    sig1[i,:] = C.dot(eps1[i,:])
print('\nSubtracted strains to force equilibrium: [{:g}, {:g}, {:g}]'.format(*EPS0))
    
# report results
plotResult(d, eps1, eeps, title='Strain distribution, force=0', fname='eps_final_force0', ytitle='Strain, 1e-6')    
plotResult(d, sig1, esig, title='Stress distribution, force=0', fname='sigma_final_force0')    
"""



