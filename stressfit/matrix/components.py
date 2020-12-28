# -*- coding: utf-8 -*-
"""
Classes describing instrument components in matrix representation.

Created on Sun Dec  2 00:01:54 2018
@author: Jan Saroun, saroun@ujf.cas.cz
"""

import numpy as np
from math import cos, sin, sqrt, log, tan

g_uni=sqrt(1/6)
# circular distribution                 
g_circ=sqrt(1/8)
# triangular distribution
g_tri=sqrt(1/3)
# gaussian distribution
g_norm=1/sqrt(4*log(2))
# exponential distribution
g_exp=sqrt(2)
# degree
deg = np.pi/180
# h/m [Ang*m/ms]
hovm = 3.95603402
GammaNi = 1.731e-3



def matrix2str(M):
    sh = M.shape
    s = "\n"
    if (len(sh)>1):
        fmt = "{:12.5g}\t"*(sh[1]-1)+"{:g}\n"
        for i in range(sh[0]):
            s += fmt.format(* M[i,:])
    else:
        fmt = "{:12.5g}\n"
        for i in range(sh[0]): 
            s += fmt.format(M[i])
    return s

def prnMat(M, id=""):
    s = matrix2str(M)
    print("\ndump {}:".format(id)+s)
    
def sign(x):
    # redefine numpy.sign so that sign(0)=1, not 0
    if (x<0):
        return -1
    else:
        return 1


def CT(L, v ):
    """Ray propagation matrix, free flight between two parallel planes
    normal to z-axis. The array elements are ordered as 
    [x, alpha, d_lambda/lambda, time]
    
    Arguments
    ----------
    
    L: float
        flight path from the sample [mm]
    v: float
        nominal velocity (for the mean wavelength)
    
    Returns
    --------
    CT: array(4,4)
        The propagation matrix.     
    """
    C = [
            [1, L, 0, 0],
            [0, 1, 0, 0 ],
            [0, 0, 1, 0],
            [0, 0, L/v, 1]
        ]        
    return np.array(C)

def CU(L, alpha, v ):
    """Ray propagation matrix from the sample up-stream.
    
    Arguments
    ----------
    
    L: float
        flight path from the sample [mm]
    alpha: float
        primary beam axis angle [rad] relative to the sample surface [rad] 
        (=theta for symmetric reflection geometry)
        geometry with incident beam on the front surface (depth=0)
    v: float
        neutron velocity
    
    Returns
    --------
    CU: array(4,5)
        The propagation matrix.
    
    """
    si = sin(alpha)
    co = cos(alpha)
    C = [
            [-si, co, -L, 0, 0],
            [0, 0, 1, 0, 0 ],
            [0, 0, 0, 1, 0],
            [co/v, si/v, 0, L/v, 1]
        ]        
    return np.array(C)
    
def CD(L, alpha, theta, v):
    """Ray propagation matrix from the sample down-stream.
    
    Arguments
    ----------
    
    L: float
        Flight path from the sample [mm]
    alpha: float
        secondary beam axis angle [rad] relative to the sample surface [rad] 
        (=theta for symmetric reflection geometry)
        geometry with incident beam on the front surface (depth=0)
    theta: float
        Bragg angle [rad] 
    v: float
        neutron velocity
    Returns
    --------
    CD: array(4,5)
        The propagation matrix.     
    
    
    """ 
    si = sin(alpha)
    co = cos(alpha)
    tt = tan(theta)
    C = [
            [si, co, L, 2*L*tt, 0],
            [0, 0, 1, 2*tt, 0 ],
            [0, 0, 0, 1, 0],
            [co/v, -si/v, 0, L/v, 1]
        ]        
    return np.array(C)


class Comp():
    """Matrix model of an abstract component, ToF version.
    
    Phase space variables are indexed in following order:
    x, z, incident angle, d_lambda/lambda, time.
    This component represents just a transparent plane. 
    Actual components like slits or collimators should be derived from Comp.
    
    Arguments
    ---------
    distance: float
        distance from the sample [mm]
    
    """
      
    def __init__(self, distance=0.0):
        self.Cin = None
        self.Cpre = None
        self.L = distance 
  
    def initialize(self, C, v0): 
        """Recalculates dependent data after change of instrument setting
        (but not distances).
        
        Arguments
        ---------
            C: array
                Transport matrix from the sample to the z=0 plane 
                of the preceding component.
            v0: float
                Mean neutron velocity [mm/us].
        """
        self.Cpre = C
        self.v0 = v0
        self.Cin = self.getCin()
        self.Cout = self.getCout()
        self.T = np.zeros((5,5))
        
    def getCin(self):        
        """Returns transport matrix from the sample to the z=0 plane of this component
        at the input. 
        """
        if (self.Cpre is None):
            C=CT(self.L, self.v0)
        else:
            C=CT(self.L, self.v0).dot(self.Cpre)
        return C
       
    def getT(self):
        """Transmission matrix"""
        return self.T

    def getCout(self):
        """Returns transport matrix from the sample to the z=0 plane of this component
        at the output. 
        """
        return self.Cin
        
    def display(self):
        #print('\nComponent: {}\n'.format(self.id))
        #print('Input matrix:\n')
        #print(self.Cin)
        #print('Output matrix:\n')
        #print(self.getCout())
        #print('Transmision matrix:\n')
        #print(self.T)
        prnMat(self.T, self.id+".T")
        

class Sample(Comp):
    """Matrix model of the sample, ToF version.
    
    Phase space variables are indexed in following order:
    x, z, incident angle, d_lambda/lambda, time.
    This is the central component of the instrument, representing an ideal pane-parallel
    polycrystalline sample. Instrument matrix is derived with respect ti its local
    coordinate system (z // depth). 
    
    Arguments
    ---------
    dhkl: float
        selected diffraction plane dhkl [AA]
    thickness: float
        sample thickness [mm]
    
    """
    def __init__(self, dhkl=1.17, thickness=10):
        self.dhkl = dhkl
        self.thickness = thickness            
        super().__init__(0.0)
        self.id='Sample'
        
    def initialize(self, omega, alpha, theta, v0):
        """Recalculates dependent data after change of input parameters
        Arguments
        ---------
            omega: float
                Surface angle [rad] reative to symmetric reflection position.
            alpha: float
                Take-off angle [rad] of the secondary beam axis.
            theta: float
                Bragg angle [rad]
            v0: float
                Horizontal component of mean neutron velocity [m/ms]
        """
        self.alpha1 = 0.5*alpha + omega
        self.alpha2 = 0.5*alpha - omega
        self.theta = theta
        self.lam0 = 2*self.dhkl*sin(abs(theta))
        self.v0 = v0
        sin1 = sin(self.alpha1)
        sin2 = sin(self.alpha2)
        ag = 1/sin1+1/sin2
        cg = 0.5*(sign(sin2)/sin2 + sign(sin1)/sin1 - ag)
        self.a = [ag, cg]           
        super().initialize(None, v0)
       
    def getCin(self):     
        """Returns transport matrix from the sample to the z=0 plane of this component
        at the input. 
        """
        C = CU(0.0, self.alpha1, self.v0)
        
        return C
        
    def getCout(self):
        """Returns transport matrix from the sample to the z=0 plane of this component
        at the output. 
        """
        C = CD(0.0, self.alpha2, self.theta, self.v0)
        return C


class Slit(Comp):
    """Matrix model of a slit, ToF version. 
    
    Arguments
    ---------
    distance: float
        distance from the sample to the slit [mm]. 
        Set a negative value if the component is from the sample up-stream.
    width: float
        slit width [mm]
    
    """             
    def __init__(self, distance=0, width=2.0):
        super().__init__(distance)
        if (distance == 0): 
            self.w = width*g_tri
        else:
            self.w = width*g_uni
        self.id = 'Slit'
        
    def getT(self):
        V = self.Cin[0,:]
        self.TT = super().getT()
        n = np.size(V)
        if (self.w > 0):
            TT = (self.w**-2)*np.kron(V,V).reshape(n,n)
            self.T = self.T + TT
        return self.T


class Guide(Comp):
    """Matrix model of a neutron guide, ToF version. 
    Width represents the constraint o beam size at the guide exit,
    m-value defines the constraint on beam divergence, assuming uniform
    distribution with FWHM = 2*m*GammaNi*wavelength.
    Set width = None to consider only the divergence constraint.
    Set m = None to consider only the beam width constraint.
    
    Arguments
    ---------
    distance: float
        Distance from the sample to the guide exit [mm]. 
        Set a negative value if the component is from the sample up-stream.
    width: float
        width of the guide exit [mm]
    m: float
        m-value of the coating
    
    """        
    def __init__(self, distance=0, width=30, m=2):
        super().__init__(distance)
        self.w = width     # width [mm] 
        self.m = m # m-value
        self.id = 'Guide'      

    def initialize(self, C, v0): 
        """Recalculates dependent data after change of instrument setting
        (but not distances).
        
        Arguments
        ---------
            C: array
                Transport matrix from the sample to the z=0 plane 
                of the preceding component.
            v0: float
                Mean neutron velocity [mm/us].
        """
        super().initialize(C,v0)
        self.lam0 = hovm/v0
        self.div = 2*self.m*GammaNi*self.lam0*g_uni
        self.Cpre = C
        self.v0 = v0
        self.Cin = self.getCin()
        self.Cout = self.getCout()
        self.T = np.zeros((5,5))
        
    def getT(self):
        V = self.Cin[0,:]
        self.T = super().getT()
        n = np.size(V)
        if (self.w):
            self.T = self.T + ((self.w*g_uni)**-2)*np.kron(V,V).reshape(n,n)
        if (self.m):
            V = self.Cin[1,:]
            self.T = self.T + (self.div**-2)*np.kron(V,V).reshape(n,n)
        return self.T


class Pulse(Comp):
    """Matrix model of a pulse source, ToF version.
    
    Arguments
    ---------
    distance: float
        distance from the sample to the source [mm]. 
        Set a negative value if the component is from the sample up-stream.    
    tau: float
        Pulse width (FWHM) [us]
    shape: float
        Pulse shape factor. (tau*shape)^2 should yield 2*sigma^2, where sigma^2 
        is the 2nd moment of the pulse distribution.
        If None, a Gaussian shape is assumed. For example:
        
        Use shape=sqrt(2) if width = decay time of an exponential pulse.
        
        Use shape=1/sqrt(3) if width = FWHM of a triangular pulse.
    
    """     
    def __init__(self, distance=0, tau=None, shape=None):
        super().__init__(distance)
        if (tau is not None):
            self.tau=tau*shape
        else:
            self.tau=tau*g_norm    # time width [us]
        self.id='Pulse'
   
    def getT(self):
        V = self.Cin[3,:]
        self.T = super().getT()
        n = np.size(V)
        if (self.tau > 0):
            self.T = self.T + (self.tau**-2)*np.kron(V,V).reshape(n,n)
        return self.T

class TofDetector(Slit):
    """Matrix model of a detector, ToF version. Defines spatial and time resolution.
    amin, amax parameters define the detector unit angular range. This is used by MXmodel to 
    average results over wavelengths corresponding to this range for given dhkl value.
    
    Arguments
    ---------
    distance: float
        distance from the sample to the detector [mm].    
    binwidth: float
        bin width [mm]
    binshape: float
        Shape factor for the bin. Use g_norm if binwidth represents FWHM of a Gaussian resolution. 
        Use g_uni for a uniform bin of this width. 
    dtau: float
        time resolution (Gaussian FWHM) [us]
    bins: list
        A list of detector channels data. Use getDetBinsCyl or getDetBinsFlat to generate the list.
    """                
    def __init__(self, distance=0, binwidth=None, binshape=None, dtau=None, bins=None):
        super().__init__(distance, binwidth*binshape/g_uni)
        # detector channels
        if (not isinstance(bins, list) or len(bins)<1):
            raise Exception("TofDetector requires that 'bins' is a non-empty list.")
        self.dist = self.L
        # time resolution (gaussian fwhm)
        self.dtau=dtau*g_norm
        self.id='TofDetector'
        # detector channels
        self.bins = bins
        # number of channels
        self.nbin = len(bins)
        # set the middle bin as the default
        n = int(0.5*(self.nbin-1))
        self.setChannel(n)


    def setChannel(self, ic):        
        """Set current channel: assign bin data for it.
        """
        self.alpha = self.bins[ic][0]
        self.theta = 0.5*self.alpha
        self.L = self.dist + self.bins[ic][1]
        self.isel = ic
    
    def getChannel(self):
        """Return alpha, theta and cos(phi) for actually selected channel"""
        return [self.alpha, self.theta]


        