# Created on Sat Dec  1 19:32:56 2018
# Written by: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
"""
Matrix model for a neutron ToF diffractometer.
Provides gauge volume parameters needed for calculation of spurious strains.


Usage
-----

    1) Create an instance, X=MXmodel_tof()
        
    2) Initialize for given setup parameters, X.setParam(inst,theta,i_sam,i_det)
       inst = array of instrument components
       theta = nominal Bragg angle [rad]
       i_sam,i_det = indexes of the Sample and TofDetector components in
       
    3) Use functions like scan_theta, scan_dhkl etc to calculate the
        gauge volume parameters as a function of scattering angle, dhkl or sample orientation:
            
        DSE = spurious strain sensitivity 
        beta = gauge volume size parameter
        zeta = attenuation parameter (depth dependence)
        fwhm = instrument resolution (delta_d/d)
 
Notes
-----
    The matrix elements are ordered as [x, z, alpha, d_lambda/lambda, time], where
    x, z:    coordinates of the scattering event (z // primary beam)
    alpha:    angular divergence of the incomming ray
    d_lambda/lambda:    wavelength deviation from the mean value
    t:    neutron delay at the scattering point with respect to the nominal time t0,
    
    t0 = Linst*m_n*lambda/h, Linst = length of the instrument (source-sample).
    
    For given diffraction line and mean wavelength, t0 is related to the Bragg angle,
    t0 = 2*Linst*m_n*dhkl/h*sin(theta)
    
    Hence the indices m and k in MXmodel_tof should be m=4 (time) and k=1 (z). 

 
"""

import numpy as np
from numpy import sqrt
from .components import matrix2str, prnMat
from math import atan, asin, cos, sin

# uniform distribution                 
g_uni = np.sqrt(1/6)
# circular distribution                 
g_circ = sqrt(1/8)
# triangular distribution
g_tri = sqrt(1/3)
# gaussian distribution
g_norm = 1/sqrt(4*np.log(2))
r8ln2 = sqrt(8*np.log(2))
# degree
deg = np.pi/180
# h/m [Ang*m/ms]
hovm = 3.95603402


def getDetBinsCyl(radius=2000.0, angle=[90.0], phi=None):
    """
    Calculate detector bin data.
    
    Parameters
    ----------
    radius : float
        detection radius [mm]
    angle : list of float
        horizontal take-off angles [deg]
    phi : list of float (optional)
        vertical take-off angles [deg] 

        
    Returns
    --------
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
    Calculate detector bin data. Flat surface.
    
    Parameters
    ----------
    distance: float
        distance of the detector centre from teh sample  [mm]
    angle: float
        take-off angle [deg] of the detector centre
    bx: list of float
        bins horizontal positions relative to the centre [mm] 
    by: list of float (optional)
        bin horizontal positions relative to the centre [mm] 
        
        
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
            L = sqrt(LH2)
            sa = bx[ih]/sqrt(LH2)
            da = asin(sa)
            alpha = a0 + da
            #ca = cos(alpha)
            #x = ca*cp            
            # tanth = (1.0 - x)/sqrt(1.0 - x**2)        
            #theta = atan(tanth)
            ch = [alpha, L - distance]
            bins.append(ch)
    return bins


class MXmodel_tof():
    """
    Matrix model of a ToF diffractometer.
    
    Parameters
    ----------
    
    flux: array
        Lookup tabl with neutron flux (2-columns with wavelength [A] and flux [rel. units])
    pulse: array
        Lookup table with source pulse widths (2-columns with wavelength [A] and pulse width [us]).
        The width is expected to express FWHM recalculated to an equivalent Gaussian variance. 
        For example, for a triangular pulse of FWHM, the width in the table should be 
        sqrt(4*ln(2))/sqrt(3) * FWHM =~ 0.96 * FWHM
        For an exponential pulse exp(-t/tau), the width should be
        sqrt(8*ln(2)) * tau =~ 2.355 * tau
        
    ext: Extinction()
        Instance of the Extinction model from stressfit
        
    Methods
    --------
    
    display:
        Dump matrices for all components and setup parameters.
    setParam:
        Set new instrument sequence and construct the T and C matrices.
    initialize:
        Matrix model initialization.
        
    """
    def __init__(self, flux=None, pulse=None, ext=None):
        self.m=3 # row index corresponding to the time-coordinate in the transport matrix (C)
        self.k=1 # row/column index corresponding to the z-coordinate in the transmision matrix (T)
        self.pulsetab=pulse
        self.fluxtab=flux
        self.ext=ext  
     
    def setParam(self, inst, i_bin, i_sam, omega):       
        """ Set new instrument sequence and construct the T and C matrices.
        
        Parameters
        ----------
        
        inst: list
            ordered list o instrument components (source to detector)
        i_bin: int
            selected detector bin index
        i_sam: int
            index to the Sample component in the list
        i_det: in
            index to the Detector component in the list
        omega: float
            sample surface angle [rad], zero for symmetric reflection
            
        """
        
        self.isam = i_sam
        self.idet = np.size(inst)-1   
        source = inst[0]
        detector = inst[self.idet]
        sample = inst[i_sam]
#TODO :        
        # select required channel
        detector.setChannel(i_bin)
        # retrieve channel data
        [alpha, theta] = detector.getChannel()
        
        
        # initialize variables
        self.omega = omega
        self.alpha = alpha
        self.theta = theta
        self.lam0 = 2*sample.dhkl*sin(abs(theta))
        self.v0 = hovm/self.lam0
        self.Ltof = 0.0
        self.T = np.zeros((5,5))
          
        # get pulse width for given wavelength
        source.tau = np.interp(self.lam0, self.pulsetab[:,0], self.pulsetab[:,1])*g_norm
        # get source flux
        self.flux = np.interp(self.lam0, self.fluxtab[:,0], self.fluxtab[:,1])
        
        # start initialization from sample:
        sample.initialize(omega, alpha, theta, self.v0)
        # get initial transport matrices for up and down stream from the sample
        # sample.display()
        C1 = sample.getCin()
        C2 = sample.getCout()
        
        # construct matrices 
        # sample -> detector
        for i in range(i_sam+1,self.idet+1):
            inst[i].initialize(C2, self.v0)
            TT = inst[i].getT()
            self.T = self.T + TT
            C2 =  inst[i].getCout()
            self.Ltof = self.Ltof + inst[i].L
        # sample -> source
        for i in reversed(range(0, i_sam)):
            inst[i].initialize(C1,self.v0)
            TT=inst[i].getT()
            self.T = self.T + TT
            C1 =  inst[i].getCout()
            self.Ltof = self.Ltof-inst[i].L

        self.C=C2
        self.a=inst[self.isam].a
        self.tof0=self.Ltof/self.v0
        self.dtau=inst[self.idet].dtau
        self.initialize()

    def initialize(self):
        """ Check transmission and transport matrices. Return true if they are
        valid (det|T|>0). 
        Calculate dependent arrays and width parameters
        it = index of time-coordinate
        iz = index of z-coordinate
        Calculate beta and the factor at <x_D>, eqs. (9) and (14)
        """
        
        DT =np.linalg.det(self.T)
        if (DT<1e-8):
           msg = matrix2str(self.T)
           raise Exception(msg+"T-matrix not invertible, det(T)<=0\n")
        res=True
        # eq. (7), define W and V as submatrices of T
        W = np.delete(np.delete(self.T,self.k,0),self.k,1)
        V = np.delete(self.T[:,self.k],self.k)
        Winv = np.linalg.inv(W)
        WV = Winv.dot(V)  # = inv(W)*V
        beta2 = self.T[self.k,self.k] - V.T.dot(WV)
        self.beta = sqrt(beta2)
        Cm = np.delete(self.C[self.m,:],self.k)
        self.DeltaX = -1.0e6*(self.C[self.m,self.k] - Cm.dot(WV))/self.tof0
        self.mu = self.ext.getMu(self.lam0)
        self.zeta = self.a[0]*self.mu/(2*self.beta)
        return res


    def gaugeWidth(self):
        """ width (fwhm) of the sampled gauge volume [mm]"""
        return 1/self.beta/g_norm
    
    def peakWidth(self):
        """peak width  delta_d/d in [1e-6] units"""
        
        Tinv = np.linalg.inv(self.T)
        TC = self.C.dot(Tinv.dot(self.C.T))
        # add detector resolution
        res = 1.0e6*sqrt(TC[self.m,self.m] + self.dtau**2)/self.tof0/g_norm
        return res

    def gaugeParams(self):
        """ Return gauge parameters as a hash map. no calculation is done, only 
        returns data calculated in setParams()
        
        Returns
        --------
        
        DSE: float
            pseudo-strain coefficients
        beta: float
            gauge size parameters
        zeta: float
            attenuation parameters
        ag: float
            Linear term for the attenuation parameter a_tot=ag*d + cg
        cg: float
            Constant term for the attenuation parameter a_tot=ag*d + cg    
        mu: float
            attenuation coefficient [1/mm]               
        p: float
            weight factor derived from the flux table
                
        """
        
        res = {}
        res['DSE'] = self.DeltaX
        res['beta'] = self.beta
        res['zeta'] = self.zeta
        res['ag'] = self.a[0]
        res['cg'] = self.a[1]
        res['mu'] = self.mu
        res['p'] = self.flux
        return res   



    def display(self,inst):       
        """ Dump matrices for all components and setup parameters.
        
        Parameters
        ----------
        
        inst: list
            ordered list o instrument components (source to detector)
        """
        
        n = np.size(inst)
        # sample -> detector
        for i in range(self.isam+1,n):
            inst[i].display()
        # sample -> source
        for i in reversed(range(0, self.isam)):
            inst[i].display()
        prnMat(self.T, "instrument.T")
        prnMat(self.C, "instrument.C")
        print('alpha = {:g}'.format(self.alpha/deg))
        print('theta = {:g}'.format(self.theta/deg))
        print('omega = {:g}'.format(self.omega/deg))
        print('Ltof = {:g}'.format(self.Ltof))
        print('a = [{:g}, {:g}]'.format(* self.a))
        print('mu = {:g}'.format(self.mu))
        print('dtau = {:g}'.format(self.dtau))
        print('lam0 = {:g}'.format(self.lam0))
        print('tof0 = {:g}'.format(self.tof0))
        print('v0 = {:g}'.format(self.v0))
        print('DeltaX = {:g}'.format(self.DeltaX))
        print('zeta = {:g}'.format(self.zeta))
        print('beta = {:g}'.format(self.beta))
