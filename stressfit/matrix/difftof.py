# -*- coding: utf-8 -*-
"""
STRESSFIT
Package for pseudostrain calculations and strain data analysis. 

Written by: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz

Created on Sat Dec  1 21:59:30 2018
@author: Jan Saroun, saroun@ujf.cas.cz
"""

import numpy as np
from .mxtof import MXmodel_tof as MX
from .scans import PSModel
from .components import Sample, Slit, Guide, Pulse, TofDetector
from stressfit import readHash, loadData
from stressfit.mx import getDetBinsCyl, getDetBinsFlat
from math import sqrt, log, tan

deg = np.pi/180
g_norm=1/sqrt(4*log(2))  # gaussian distribution
g_uni=sqrt(1/6)  # uniform distribution
        
class INSTR():
    """
    Instrument object: time-of-flight strain diffractometer.
    
    Uses: 
        
    MXmodel_tof.py:
        matrix model of the instrument
    
    PSModel.py: 
        pseudostrain calculations
    
    The instrument can be initialized from a parameter file (see resources/ENGINX.par
    for an example). The instrumentg is composed from following components:
    
    - pulsed source
    - divergence slit or a neutron guide
    - primary slit or radial collimator
    - polycrystalline sample, planparallel plate
    - secondary slit or radial collimator
    - position sensitive ToF detector
    
    Lookup tables can be provided for pulse shape and spectrum 
    (see the examples pulse.dat and spectrum.dat in ./resources).
    """    
    def __init__(self, flux=None, pulse=None, extinction=None, param=None):
        """
        Constructor of the instrument. The instrument is initialized from 
        a file with parameters (see resources/ENGINX.par as an example).
        
        Call initialize(...) to construct the instrument matrix model and 
        dependencies before starting calculations. 
        
        Call setParam(...) to update parameters and the matrix model.
        
        Arguments:
        -----------
        
        flux: str
            file with flux table (wavelength [AA], flux [rel/ units])
        pulse: str
            file with pulse shape table (time [us], flux [rel/ units])
        extinction: float
            Instance of the Extinction() class (see extinction.py)    
        param: str
            file with instrument description data. Format is [name value]
            
        """
        pulse_table = loadData(pulse)
        flux_table = loadData(flux)
        self.par = readHash(param)
        self.parfile = param        
        self.MX = MX(flux=flux_table, pulse=pulse_table, ext=extinction)

    def initialize(self, sample, bins=None, **kwargs):  
        """ Instrument setup initialization.
        
        Parameters:
        -----------
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
        self.PS = PSModel()
        self.PS.setParam(thickness, b, c)

        if (bins == None):
            # by default, assume flat detector with 11 horizontal channels, 1 vertical
            dist = self.par['d_dist']
            da = 0.5*abs(self.par['d_amax']-self.par['d_amin'])
            angle = 0.5*(self.par['d_amin']+self.par['d_amax'])
            dx = dist*tan(da*deg)
            nx = 11
            bx = np.linspace(-dx, dx, num=nx)           
            self.bins=getDetBinsFlat(distance=self.par['d_dist'], angle=angle, bx=bx)
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
        """ Updates the instrument model for given parameters and bin index.
        """ 
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
        sr.append(Pulse(srcL,p1, p2))
        # DIVERGENCE SLIT
        # parameters: distance [mm], wavelength [A], width[mm]
        if (self.par['s0_on']): 
            sr.append(Slit(-par['s0_dist']+par['s1_dist'],par['s0_width']))
        else:
            sr.append(Guide(-par['gd_dist']+par['s1_dist'],par['gd_width'],par['gd_m']))  
    
        # INPUT SLIT or RADIAL COLLIMATOR
        # parameters:  distance [mm], wavelength [A], width[mm]
        sr.append(Slit(-par['s1_dist'],par['s1_width']))

        # parameters: dhkl [A], surface angle [rad], thickness [mm]
        sr.append(Sample(self.dhkl, self.thickness))
        self.i_sam = len(sr)-1 

        # OUTPUT SLIT or RADIAL COLLIMATOR
        # parameters: distance [mm], wavelength [A], width[mm]
        sr.append(Slit(par['s2_dist'], par['s2_width']))

        # DETECTOR
        det = TofDetector(
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
        self.dhkl = dhkl
        self.src[self.i_sam].dhkl = dhkl
        self.initMX()
    
    def setOmega(self,omega):
        """set omega angle [deg]"""
        self.omega=omega*deg
        self.initMX()      
        
    def setBin(self, i_bin):
        det = self.src[self.i_det]
        det.setChannel(i_bin)
        self.i_bin = det.isel
        self.initMX()
        
    def initMX(self):
        self.MX.setParam(self.src, self.i_bin, self.i_sam, self.omega) 
        self.lam0 = self.MX.lam0
        self.alpha = self.MX.alpha
        self.theta = self.MX.theta

    """
    ----------------------------------
    Matrix calculations.
    ----------------------------------
    """     

    def execFnc(self, fnc, i_bin=None, *args, **kwargs):
        """ A wrapper which will run a function either for given detector bin or averaged over all bins. 
        If i_bin==None, average over all bins is received.
        fnc() must return a hash map including 'p' with weight factor and 
        some numerical values. See fnc_fwhm or fnc_gauge.
        *args, **kwargs are arguments passed to fnc. 
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
    
    """
    --------------------------------------------------
    Shortcuts calling execFnc
    ------------------------------------------------------
    """

    def getGaugeParams(self, i_bin=None):
        """get gauge parameters DSE, beta, zeta"""    
        def fnc_gauge():
            res = self.MX.gaugeParams()
            return res       
        res = self.execFnc(fnc_gauge, i_bin=i_bin)
        return res
    
    def getFWHM(self, i_bin=None):
        """Return FWHM of the diffraction line (uperturbed)"""
        def fnc_fwhm():
            res = {}
            res['fwhm'] = self.MX.peakWidth()
            res['p'] = self.MX.flux
            return res
        res = self.execFnc(fnc_fwhm, i_bin=i_bin)
        return res
    
    def scanDepth(self, depth=None, keys=None, i_bin=None):
        """Get peak shift,gauge centre and gauge width as a function of scan depth.
        
        Arguments
        --------
            depth: array
                depth values
            keys: list
                choose what to calculate. It can contain: 'shift', 'zscale' and 'width'.
                If empty or none, calculate all three functions
        Returns
        -------
            DSE, beta, zeta and required arrays as a hash map
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
        """ Resolution curve, delta_d/d as a function of dhkl"""
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
        """ Gauge parameters as a function of sample orientation
        return gauge parameters for each omega value
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
        """ Return gauge parameters for each dhkl value 
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
           
    def display(self):
        self.MX.display(self.src)

