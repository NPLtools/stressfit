# -*- coding: utf-8 -*-
# Created on Mon Dec  3 08:54:52 2018
# @author: Jan Saroun, saroun@ujf.cas.cz
"""
Gaussian model for calculation of spurious strains.

Implements functions for calculation of pseudo-strains 
(also called spurious strains) which result from perturbation of sampling 
gauge volume (e.g. surface effect).

Calculates pseudo-strains and gauge senter as a function of scan depth 
for given gauge volume parameters:
    - `DSE` = pseudo-strains sensitivity [mm-1]
    - `beta` = gauge volume size parameter [mm-1]
    - `zeta` = attenuation parameter (depth dependence)
 Use the matrix model (like MX_model_tof) to get these parameters
 analytically.

Usage
-----
    1. create an instance, X=PSModel()
    2. obtain the gauge parameters DSE, beta, zeta for given setup
        a) analytically: DSE, beta, zeta are calculated by MXmodel
        b) experimentally: beta and zeta can be fitted on intensity 
        vs. depth scan
        c) simulation: DSE, beta can be obtained by fitting simulated 
        depth scans
    3. use functions get_shift and get_ctr to calculate spurious stran 
    and sampling centre of mass for given DSE, beta, zeta and an array 
    of depth values.
 
Notes
-----
Use DSE, beta, zeta as arrays to calculate multiple shift functions 
simultaneously.

"""

import numpy as np

# gaussian distribution
g_norm=1/np.sqrt(4*np.log(2))


def meanp(A,p):
    """
    Calculate mean of values in A.
    
    Parameters
    ----------
    A: ndarray
        Values to be averaged.
    p: ndarray
        Weights for averaging.
    
    """
    return A.T.dot(p) 


# Fn(u,d) functions are n-th moments of the gauge volume distribution:
# Integral{z^n exp(-(z-u)^2)}dz from 0 to d
# in relative length units:
# h=depth*beta, u = h - zeta
# zeta=a/2/beta
# d=thickness*beta
#

def F0(u, d):
    """Calculate 0-th moment of gauge volume distribution."""
    return np.sqrt(np.pi)/2*(np.erf(u)-np.erf(u-d))

def F1(u, d):
    """Calculate 1-st moment of gauge volume distribution."""
    return 0.5*(np.exp(-u**2)-np.exp(-(u-d)**2)) + u*F0(u,d)

def F2(u, d):
    """Calculate 2-nd moment of gauge volume distribution."""
    return 0.5*F0(u,d) + u*F1(u,d) - d*0.5*np.exp(-(u-d)**2) 
 
def F3(u, d):
    """Calculate 3-rd moment of gauge volume distribution."""
    return 0.5*u*F0(u,d) + (u**2+1)*F1(u,d) - 0.5*d*(d+u)*np.exp(-(u-d)**2) 


def gwidth(u,d):
    """Calculate width (FWHM) of the sampled gauge volume.
    
    Parameters
    ----------
    u: float
         u = beta*h - zeta, where h is distance [mm] from the front surface in (< 0 outside the sample), 
         zeta is the attenuation coefficient (incl. geometry factor)
    d: float
         sample thickness, [1/beta] units
    """
    f2 = F2(u,d)
    f1 = F1(u,d)
    f0 = F0(u,d)
    sig2 = f2/f0-(f1/f0)**2
    res = np.sqrt(sig2*2)/g_norm
    return res
  

def gvol(u, d, zeta, b, c):
    """Volume of the sampled gauge volume ~ intensity.
    
    Parameters
    ----------
    u: float
         u = beta*h - zeta, where h is distance [mm] from the front surface in (< 0 outside the sample), 
         zeta is the attenuation coefficient (incl. geometry factor)
    d: float
         sample thickness, [1/beta] units
    zeta: float
         attenuation parameter [dimensionless]
    b, c: float
        scattering probability coefficients: p = 1 + b*u + c*u^2
    """
    f0 = F0(u,d)
    val = f0
    if (b != 0):
        f1 = F1(u,d)
        val = val + b*f1
    if (c!=0):
        f2 = F2(u,d) 
        val = val + c*f2
    res = np.exp(-zeta*(2*u+zeta))*val        
    return res
     

def gcentre( u, d, b, c ):
    """Centre of mass of the sampling volume.
    
    Includes inhomogeneity correction. Assumes relative length units, 1/beta:   
    
    Parameters
    ----------
    u: float
         u = beta*h - zeta, where h is distance [mm] from the front surface in (< 0 outside the sample), 
         zeta is the attenuation coefficient (incl. geometry factor)
    d: float
         sample thickness, [1/beta] units
    b, c: float
        scattering probability coefficients: p = 1 + b*u + c*u^2
    """
    # approximation for large u
    if (u-d > 5):
        val = d
    elif (u < -5):
        val = 0
    else:        
        f0 = F0(u,d)
        f1 = F1(u,d)            
        if (b == 0 and c == 0):
            val = f1/f0
        else:
            f2 = F2(u,d)
            f3 = F3(u,d)
            val = (f1+b*f2+c*f3)/(f0+b*f1+c*f2)
    return val


class PSModel():
    """Class encapsulating top-level functions for analytical pseudo-strain model."""
    
    def setParam(self,dsam,b,c):
        """
        Set sample parameters.
        
        Parameters
        ----------
        d: float
            thickness [mm]
        b,c: float
            inhomogeneity parameters, assumig quadratic weight function: f(u) = 1 + b*u + c*u^2 
        """
        self.dsam = dsam
        self.b = b
        self.c = c

    def get_shift(self, h, DSE, beta, zeta):
        """Spurious strain as a function of scan depth.
        
        (h>0 is under front surface)
        
        Parameters
        ----------
        h: list or array
            depth values [mm]
        DSE: list or array
            values of surface effect amplitude [1e-6/mm]
        beta: list or array
            values of gauge volume size parameter [1/mm]
        zeta:  list or array
            values of attenuation parameter (dimensionless)
        
        DSE,beta,zeta must have the same length.
        
        
        Returns
        -------
        A 2D array with columns z, shift[0], shift[1], ...
        
        The shift values are in 1e-6 units and correspond to the values of DSE, beta, zeta.
        """
        n = np.size(h)
        nb = np.size(DSE)
        if (nb==1 and isinstance(DSE, float)):
            xDSE = [DSE]
            xbeta = [beta]
            xZeta = [zeta]
        else:
            xDSE = DSE
            xbeta = beta
            xZeta = zeta
        eshift = np.zeros((n,nb))
        for i in range(n):
            for j in range(nb):
                hh = h[i]*xbeta[j]-xZeta[j]
                d = self.dsam*xbeta[j]
                eshift[i,j] = xDSE[j]*(gcentre(hh,d,self.b,self.c)/xbeta[j] - h[i])
        return eshift

    def get_shift_ave(self,h,DSE,beta,zeta,p):
        """
        Average of get_shift result using the weight factors p.
        
        Calculates spurious strain as a function of scan depth 
        (h>0 is under front surface).
        
        Parameters
        ----------
        h: list or array
            depth values [mm]
        DSE: list or array
            values of surface effect amplitude [1e-6/mm]
        beta: list or array
            values of gauge volume size parameter [1/mm]
        zeta:  list or array
            values of attenuation parameter (dimensionless)
        p: list or array
        
        DSE, beta, zeta, p must have the same length.
        
        
        Returns
        -------
        A 2D array with 2 columns z, shift
        
        The shift values are in 1e-6 units and correspond to the pseudo-strain 
        averaged over the given values of DSE, beta, zeta with weight factors p.
        """
        n = len(h)
        es = self.get_shift(h,DSE,beta,zeta)
        eshift = np.zeros(n)
        for i in range(n):
            eshift[i] = es[i,:].dot(p)
        return eshift/np.sum(p)
    
    def get_ctr(self,h,beta,zeta):
        """Gauge centre as a function of scan depth.
        
        (h>0 is under front surface)
        
        Parameters
        ----------
        h: list or array
            depth values [mm]
        beta: list or array
            values of gauge volume size parameter [1/mm]
        zeta:  list or array
            values of attenuation parameter (dimensionless)
        
        DSE,beta,zeta must have the same length.
        
        
        Returns
        -------
        A 2D array with columns z, centre[0], centre[1], ...
        
        The gauge centre values are in [mm] and correspond to the values of beta, zeta.
        """
        n = len(h)
        nb = len(beta)
        ctr = np.zeros((n,nb))
        for i in range(n):
            for j in range(nb):
                hh = h[i]*beta[j]-zeta[j]
                d = self.dsam*beta[j]
                ctr[i,j] = gcentre(hh,d,self.b,self.c)/beta[j]
        return ctr
    
    
    def get_ctr_ave(self,h,beta,zeta,p):
        """
        Average of get_ctr function of scan depth.
        
        (h>0 is under front surface)
        
        Parameters
        ----------
        h: list or array
            depth values [mm]
        beta: list or array
            values of gauge volume size parameter [1/mm]
        zeta:  list or array
            values of attenuation parameter (dimensionless)
        p: list or array
        
        beta, zeta, p must have the same length.
        
        
        Returns
        -------
        A 2D array with 2 columns z, gauge centre
        
        The gauge centre values are in [mm] and correspond to the gauge centre
        averaged over the given values of beta, zeta with weight factors p.
        """
        n = len(h)
        xctr = self.get_ctr(h,beta,zeta)
        ctr = np.zeros(n)
        for i in range(n):
            ctr[i] = xctr[i,:].dot(p)
        return ctr/np.sum(p)
        
    def get_gwidth(self,h,beta,zeta):
        """
        Gauge width as a function of scan depth.
        
        (h>0 is under front surface)
        
        Parameters
        ----------
        h: list or array
            depth values [mm]
        beta: list or array
            values of gauge volume size parameter [1/mm]
        zeta:  list or array
            values of attenuation parameter (dimensionless)
        
        DSE,beta,zeta must have the same length.
        
        
        Returns
        -------
        A 2D array with columns z, centre[0], centre[1], ...
        
        The gauge width values are in [mm]  and correspond to the values of beta, zeta.
        """   
        n = len(h)
        nb = len(beta)
        gw = np.zeros((n,nb))
        for i in range(n):
            for j in range(nb):
                hh = h[i]*beta[j]-zeta[j]
                d = self.dsam*beta[j]
                gw[i,j] = gwidth(hh,d)/beta[j]
        return gw
    
    def get_gwidth_ave(self,h,beta,zeta,p):
        """
        Average of get_gwidth function of scan depth.
        
        (h>0 is under front surface)
        
        Parameters
        ----------
        h: list or array
            depth values [mm]
        beta: list or array
            values of gauge volume size parameter [1/mm]
        zeta:  list or array
            values of attenuation parameter (dimensionless)
        p: list or array
        
        beta, zeta, p must have the same length.
        
        
        Returns
        -------
        A 2D array with 2 columns z, gauge FWHM
        
        The gauge FWHM values are in [mm] and correspond to the FWHM
        averaged over the given values of beta, zeta with weight factors p.
        """ 
        n = len(h)
        gw = self.get_gwidth(h,beta,zeta)
        gww = np.zeros(n)
        for i in range(n):
            gww[i] = gw[i,:].dot(p)
        return gww/np.sum(p)

