# -*- coding: utf-8 -*-
# `Written by`: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
"""
Analytical instrument model for pseudostrain calculations - ToF version.

Implements the class DiffTOF for a time-of-flight strain diffractometer.

Uses
----
stressfit.matrix.mxtof:
    Matrix model of the instrument.
stressfit.matrix.psmodel:
    Analytical calculations of pseudo-strains, using Gaussian sampling model.
"""
import numpy as np
#from math import sqrt, log, tan

import stressfit.matrix.mxtof as mxtof
import stressfit.matrix.psmodel as psm
import stressfit.matrix.components as com
import stressfit.dataio as io

deg = np.pi/180
g_norm=1/np.sqrt(4*np.log(2))  # gaussian distribution
g_uni=np.sqrt(1/6)  # uniform distribution



def getDetBinsCyl(radius=2000.0, angle=[90.0], phi=None):
    """
    Calculate detector bin data - cylindrical surface.
    
    Parameters
    ----------
    radius : float
        detection radius [mm]
    angle : list of float
        horizontal take-off angles [deg]
    phi : list of float (optional)
        vertical take-off angles [deg] 

        
    Returns
    -------
    list
        A list of bin data as [alpha, dL], where alpha [rad] is the take-off angle 
        and L is the distance difference with respect to the detector centre [mm]. 
    
    """
    bins = []
    if (phi is None):
        nv = 1
        #cp = [1.]
    else:
        nv = len(phi)
        #cp = np.cos(np.multiply(phi,deg))
    alpha = np.multiply(angle,deg)
    #ca = np.cos(alpha)
    nh = len(alpha)
    for iv in range(nv):
        for ih in range(nh):               
            #x = ca[ih]*cp[iv]
            #tanth = (1.0 - x)/sqrt(1.0 - x**2)        
            #theta = atan(tanth)
            ch = [alpha[ih], 0.0]
            bins.append(ch)
    return bins


def getDetBinsFlat(distance=2000.0, angle=90, bx=[0.0], by=None):
    """
    Calculate detector bin data - flat surface.
    
    Parameters
    ----------
    distance: float
        distance of the detector centre from teh sample  [mm]
    angle: float
        take-off angle [deg] of the detector centre
    bx: list of float
        bins horizontal positions relative to the centre [mm] 
    by: list of float (optional)
        bins vertical positions relative to the centre [mm] 
               
    Returns
    -------
    list
    
        A list of bin data as [alpha, dL], where alpha [rad] is the take-off angle 
        and L is the distance difference with respect to the detector centre [mm].
    
    """
    a0 = angle*deg
    bins = []
    if (by is None):
        nv = 1.
        by =[0.0]
    else:
        nv = len(by)
    nh = len(bx)
    for iv in range(nv):
        for ih in range(nh):   
            LH2 = distance**2 + bx[ih]**2
            #L2 = LH2 + by[iv]          
            #cp = sqrt(LH2/L2)
            L = np.sqrt(LH2)
            sa = bx[ih]/np.sqrt(LH2)
            da = np.asin(sa)
            alpha = a0 + da
            #ca = cos(alpha)
            #x = ca*cp            
            # tanth = (1.0 - x)/sqrt(1.0 - x**2)        
            #theta = atan(tanth)
            ch = [alpha, L - distance]
            bins.append(ch)
    return bins


class DiffTOF():
    """
    Encapsulates a model for time-of-flight strain diffractometer.
    
    The instrument can be initialized from a parameter file (see resources/ENGINX.par
    for an example). The instrumentg is composed from following components:
    
    - pulsed source
    - divergence slit or a neutron guide
    - primary slit or radial collimator
    - polycrystalline sample, planparallel plate
    - secondary slit or radial collimator
    - position sensitive ToF detector
    
    Parameters
    ----------
    flux: str
        File with flux table (wavelength [AA], flux [rel/ units])
    pulse: str
        File with pulse shape table (time [us], flux [rel/ units])
    extinction: float
        Instance of the Extinction() class (see extinction.py)    
    param: str
        File with instrument description data. Format is [name=value]
        
    Notes
    -----
    Lookup tables can be provided for pulse shape and spectrum 
    (see the examples pulse.dat and spectrum.dat in stressfit.resources).
    
    The instrument is initialized from a file with parameters 
    (see resources/ENGINX.par as an example).
        
    Call initialize(...) to construct the instrument matrix model and 
    dependencies before starting calculations. 
        
    Call setParam(...) to update parameters and the matrix model.
    
    """ 
    
    def __init__(self, flux=None, pulse=None, extinction=None, param=None):
        pulse_table = io.load_data(pulse)
        flux_table = io.load_data(flux)
        self.par = io.read_dict(param)
        self.parfile = param        
        self.MX = mxtof.MX(flux=flux_table, pulse=pulse_table, ext=extinction)

    def initialize(self, sample, bins=None, **kwargs):  
        """
        Instrument setup initialization.
        
        Parameters
        ----------
        sample: hash map
            Sample parameters in a hash map, containing at least:
            
            omega: float
                sample rotation [deg] (=0 for symmetric reflection)
                
            dhkl: float
                sample dhkl [A]
                
            thickness: float
                sample thickness [mm]
                
        bins: list
            A list of detector channels data. 
            Use getDetBinsCyl or getDetBinsFlat from MXmodel_tof to generate the list.                       
        kwargs:
            Other named parameters: overrides settings loaded form the parameter file
            on construction. 
            It handles also an argument slits = [s0w, s0d, s1w, s1d, s2w, 0.] 
            with collimation parameters passed as a list.
        """
        omega = sample['omega']
        dhkl = sample['dhkl']
        thickness = sample['thickness']
        b = 0
        c = 0
        if 'b' in sample.keys():
            b = sample['b']
        if 'c' in sample.keys():
            c = sample['c']
        
        # create peak shift model for given sample data
        self.PS = psm.PSModel()
        self.PS.setParam(thickness, b, c)

        if (bins == None):
            # by default, assume flat detector with 11 horizontal channels, 1 vertical
            dist = self.par['d_dist']
            da = 0.5*abs(self.par['d_amax']-self.par['d_amin'])
            angle = 0.5*(self.par['d_amin']+self.par['d_amax'])
            dx = dist*np.tan(da*deg)
            nx = 11
            bx = np.linspace(-dx, dx, num=nx)           
            self.bins = getDetBinsFlat(distance=self.par['d_dist'], angle=angle, bx=bx)
        else:
            self.bins = bins
        
        if kwargs:
            for key in kwargs.keys():
                # print("{} == {}".format(key,value))
                value = kwargs[key]
                if (key == 'slits' and np.size(value)>5):
                    self.par['s0_width'] = value[0]
                    self.par['s0_dist'] = value[1]
                    self.par['s1_width'] = value[2]
                    self.par['s1_dist'] = value[3]
                    self.par['s2_width'] = value[4]
                    self.par['s2_dist'] = value[5]
                else:
                    self.par[key] = value

        # Configuration variables:
        self.omega = omega*deg
        self.dhkl = dhkl
        self.thickness = thickness
        # initialize the model - calculate matrices
        self.setParam(self.par)

       
    def setParam(self, par, i_bin=None):
        """Update the instrument model for given parameters and bin index.""" 
        self.par = par
        # SOURCE
        # parameters: distance [mm], wavelength [A], tau [us], shape (see NOTE)
        if (par['s0_on']):
            srcL=-par['src_dist']+par['s0_dist']
        else:
            srcL=-par['src_dist']+par['gd_dist']  
    
        p1=par['pulse_width']
        p2=par['pulse_shape']*g_norm
        sr = []
        sr.append(com.Pulse(srcL,p1, p2))
        # DIVERGENCE SLIT
        # parameters: distance [mm], wavelength [A], width[mm]
        if (self.par['s0_on']): 
            sr.append(com.Slit(-par['s0_dist']+par['s1_dist'],par['s0_width']))
        else:
            sr.append(com.Guide(-par['gd_dist']+par['s1_dist'],par['gd_width'],par['gd_m']))  
    
        # INPUT SLIT or RADIAL COLLIMATOR
        # parameters:  distance [mm], wavelength [A], width[mm]
        sr.append(com.Slit(-par['s1_dist'],par['s1_width']))

        # SAMPLE
        # parameters: dhkl [A], surface angle [rad], thickness [mm]
        sr.append(com.Sample(self.dhkl, self.thickness))
        self.i_sam = len(sr)-1 

        # OUTPUT SLIT or RADIAL COLLIMATOR
        # parameters: distance [mm], wavelength [A], width[mm]
        sr.append(com.Slit(par['s2_dist'], par['s2_width']))

        # DETECTOR
        det = com.TofDetector(
                distance = par['d_dist'] - self.par['s2_dist'],
                binwidth = par['d_binwidth'],
                binshape = par['d_binshape']*g_uni,
                dtau = par['d_resol'],
                bins = self.bins)
        if (i_bin==None):
            self.i_bin = det.isel
        else:
            det.setChannel(i_bin)
            self.i_bin = det.isel       
        sr.append(det)
        self.i_det=len(sr)-1
        self.src=sr
           
        # update the matrix model
        self.initMX()

    def setDhkl(self,dhkl):
        """Set dhkl [Ang]."""
        self.dhkl = dhkl
        self.src[self.i_sam].dhkl = dhkl
        self.initMX()
    
    def setOmega(self,omega):
        """Set omega angle [deg]."""
        self.omega=omega*deg
        self.initMX()      
        
    def setBin(self, i_bin):
        """Set default detecotr bin index."""
        det = self.src[self.i_det]
        det.setChannel(i_bin)
        self.i_bin = det.isel
        self.initMX()
        
    def initMX(self):
        """Update matrix model for current instrument configuration."""
        self.MX.setParam(self.src, self.i_bin, self.i_sam, self.omega) 
        self.lam0 = self.MX.lam0
        self.alpha = self.MX.alpha
        self.theta = self.MX.theta

###    Matrix calculations.

    def execFnc(self, fnc, i_bin=None, *args, **kwargs):
        """
        Wrap a coll to the methods fnc_fwhm, fnc_gauge etc.
        
        It runs a function either for given detector bin or averaged over 
        all bins, depending on the value of i_bin.
        
        Parameters
        ----------    
        fnc: function
            fnc() must return a hash map including 'p' with weight factor and 
            some numerical values. See fnc_fwhm or fnc_gauge for an example.        
        i_bin: int
            A detector bin index. If i_bin==None, calculate average over 
            all detector bins.
        args, kwargs:
            Positional and keyword arguments passed to fnc.
        """
        if (i_bin==None):
            ibin0 = self.i_bin
            bins = self.src[self.i_det].bins
            nb = len(bins)
            p = np.zeros(nb)
            XX = []
            for ib in range(nb):
                self.setBin(ib)
                X = fnc(*args, **kwargs)
                p[ib] = X['p']
                XX.append(X)
            w = p/np.sum(p)
            res={}
            keys = XX[0].keys()           
            for key in keys:
                # this will preserve type of  res[key]: float or array
                res[key] = 0*XX[0][key]
                for i in range(len(XX)):
                    res[key] += XX[i][key]*w[i]
            self.setBin(ibin0)
        else:
            ibin0 = self.i_bin
            self.setBin(i_bin)
            res = fnc(*args, **kwargs)
            self.setBin(ibin0)
        return res
    
###    Shortcuts calling execFnc

    def getGaugeParams(self, i_bin=None):
        """Get gauge parameters DSE, beta, zeta."""    
        def fnc_gauge():
            res = self.MX.gaugeParams()
            return res       
        res = self.execFnc(fnc_gauge, i_bin=i_bin)
        return res
    
    def getFWHM(self, i_bin=None):
        """Return FWHM of the diffraction line (uperturbed)."""
        def fnc_fwhm():
            res = {}
            res['fwhm'] = self.MX.peakWidth()
            res['p'] = self.MX.flux
            return res
        res = self.execFnc(fnc_fwhm, i_bin=i_bin)
        return res
    
    def scanDepth(self, depth=None, keys=None, i_bin=None):
        """Get peak shift, gauge centre and gauge width.

        Parameters
        ----------
            depth: array
                Depth values.
            keys: list
                Choose what to calculate. It can contain: 'shift', 'zscale' 
                and 'width'. If empty or none, calculate all three functions.
                
        Returns
        -------
            DSE, beta, zeta and required arrays as a hash map.
        """
        def fnc_scan():
            par = self.MX.gaugeParams()
            DSE = par['DSE']
            beta = par['beta']
            zeta = par['zeta']
            if  ((keys==None) or (len(keys)<1)):
                kset = ['shift', 'zscale', 'width']
            else:
                kset = keys
            res= {'DSE': DSE, 'beta': beta, 'zeta': zeta}
            for key in kset:
                if (key=='shift'):             
                    res['shift'] = self.PS.get_shift(depth, [DSE], [beta], [zeta])[:,0]
                if (key=='zscale'):             
                    res['zscale'] = self.PS.get_ctr(depth, [beta], [zeta])[:,0]
                if (key=='width'):             
                    res['width'] = self.PS.get_gwidth(depth, [beta], [zeta])[:,0]
            res['p'] = self.MX.flux
            return res
        res = self.execFnc(fnc=fnc_scan, i_bin=i_bin)
        return res
    
    def getResolution(self, dhkl, i_bin=None):
        """Resolution curve, delta_d/d as a function of dhkl."""
        d0 = self.dhkl
        n = np.size(dhkl)
        res = np.zeros((n,2))
        for i in range(n):
            self.setDhkl(dhkl[i])
            x = self.getFWHM(i_bin=i_bin)
            res[i,0] = dhkl[i]
            res[i,1] = x['fwhm']
        self.setDhkl(d0)
        return res

    def scanOmega(self, omega, i_bin=None):
        """Gauge parameters as a function of sample orientation.
        
        Returns
        -------
        ndarray
            Gauge parameters for each omega value as a 2D array.
            The 1st index links to the omega values, the 2nd index 
            denotes the variable: 0: omega, 1: DSE, 2: beta, 3: zeta.
        """
        om0 = self.omega
        n = np.size(omega)
        res = np.zeros((n,4))
        for i in range(n):
            self.setOmega(omega[i])
            x = self.getGaugeParams(i_bin=i_bin)
            res[i,0] = omega[i]
            res[i,1] = x['DSE']
            res[i,2] = x['beta']
            res[i,3] = x['zeta']  
        self.setOmega(om0)
        return res
    
    def scanDhkl(self,dhkl, i_bin=None):
        """Calculate gauge parameters for each dhkl value.
        
        Returns
        -------
        ndarray
            Gauge parameters for each dhkl value as a 2D array.
            The 1st index links to the dhkl values, the 2nd index 
            denotes the variable: 0: dhkl, 1: DSE, 2: beta, 3: zeta.
        """
        d0 = self.dhkl
        n = np.size(dhkl)
        res = np.zeros((n,4))
        for i in range(n):
            self.setDhkl(dhkl[i])
            x = self.getGaugeParams(i_bin=i_bin)
            res[i,0] = dhkl[i]
            res[i,1] = x['DSE']
            res[i,2] = x['beta']
            res[i,3] = x['zeta']  
        self.setDhkl(d0)
        return res



