# -*- coding: utf-8 -*-
# Created on Sun Dec  2 00:01:54 2018
# @author: Jan Saroun, saroun@ujf.cas.cz
"""Classes describing instrument components in matrix representation."""

import numpy as np
import math

#TODO: export all constants to a shared module stressfit.const

#TODO: upgrade to allow for 2D detector bin arrays

g_uni = math.sqrt(1/6)
# circular distribution                 
g_circ = math.sqrt(1/8)
# triangular distribution
g_tri = math.sqrt(1/3)
# gaussian distribution
g_norm = 1/math.sqrt(4*math.log(2))
# exponential distribution
g_exp = math.sqrt(2)
# degree
deg = np.pi/180
# h/m [Ang*m/ms]
hovm = 3.95603402
GammaNi = 1.731e-3


def _sign(x):
    """Redefine numpy.sign so that sign(0)=1, not 0."""
    if (x<0):
        return -1
    else:
        return 1


def _matrix2str(M):
    """Convert numpy matrix to a string, columne separated by ``|``.""" 
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
    """Print matrix on console."""
    s = _matrix2str(M)
    print("\ndump {}:".format(id)+s)

    
def CT(L, v ):
    """
    Create ray propagation matrix.
    
    Describes free flight between two parallel planes,
    normal to z-axis. The array elements are ordered as 
    [x, alpha, d_lambda/lambda, time]
    
    Parameters
    ----------
    L: float
        flight path from the sample [mm]
    v: float
        nominal velocity (for the mean wavelength)
    
    Returns
    -------
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
    """
    Ray propagation matrix from the sample up-stream.
    
    Rows correspond to [x, alpha, d_lambda/lambda, time]
    Columns correspond to [x, z, alpha, d_lambda/lambda, time]
    
    Parameters
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
    -------
    CU: array(4,5)
        The propagation matrix.
    """
    si = math.sin(alpha)
    co = math.cos(alpha)
    C = [
            [-si, co, -L, 0, 0],
            [0, 0, 1, 0, 0 ],
            [0, 0, 0, 1, 0],
            [co/v, si/v, 0, L/v, 1]
        ]        
    return np.array(C)
 
   
def CD(L, alpha, theta, v):
    """
    Ray propagation matrix from the sample down-stream.
    
    Rows correspond to [x, alpha, d_lambda/lambda, time]
    Columns correspond to [x, z, alpha, d_lambda/lambda, time]
    
    Parameters
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
    -------
    CD: array(4,5)
        The propagation matrix.     
    
    """ 
    si = math.sin(alpha)
    co = math.cos(alpha)
    tt = math.tan(theta)
    C = [
            [si, co, L, 2*L*tt, 0],
            [0, 0, 1, 2*tt, 0 ],
            [0, 0, 0, 1, 0],
            [co/v, -si/v, 0, L/v, 1]
        ]        
    return np.array(C)


class Sample():
    """
    Matrix model of the sample, ToF version.
    
    Phase space variables are indexed in following order:
    x, z, incident angle, d_lambda/lambda, time.
    This is the central component of the instrument, representing an ideal pane-parallel
    polycrystalline sample. Instrument matrix is derived with respect to its local
    coordinate system (z // depth). 
    
    Parameters
    ----------
    dhkl: float
        selected diffraction plane dhkl [AA]
    
    Notes
    -----
    Always call ``initialize`` after construction, before using any property. 
    
    """
    
    def __init__(self, dhkl=1.17):
        self.dhkl = dhkl            
        self.id='Sample'
        self.__isvalid = False 
        
    def initialize(self, omega, alpha, theta, v0):
        """
        Recalculate dependent data after change of input parameters.
        
        Parameters
        ----------
        omega: float
            Surface angle [rad] reative to symmetric reflection position.
        alpha: float
            Take-off angle [rad] of the secondary beam axis.
        theta: float
            Bragg angle [rad]
        v0: float
            Horizontal component of mean neutron velocity [m/ms]
        """
        alpha1 = 0.5*alpha + omega
        alpha2 = 0.5*alpha - omega
        sin1 = math.sin(alpha1)
        sin2 = math.sin(alpha2)
        ag = 1/sin1+1/sin2
        cg = 0.5*(_sign(sin2)/sin2 + _sign(sin1)/sin1 - ag)
        self.__a = [ag, cg]  
        self.__Cin = CU(0.0, alpha1, v0)    
        self.__Cout = CD(0.0, alpha2, theta, v0) 
        # indicate that initialization has been called
        self.__isvalid = True

       
    @property
    def Cin(self):     
        """
        Input transport matrix.
        
        Transport from the sample to the z=0 plane of this component
        at the input. 
        """      
        return self.__Cin 
    
    @property  
    def Cout(self):
        """
        Output transport matrix.
        
        Transport from the sample to the z=0 plane of this component 
        at the output. 
        """
        return self.__Cout 

    @property  
    def a(self):
        """Coefficients [ag, cg] for calculation of beam attenuation."""
        return self.__a 
 
    @property  
    def isvalid(self):
        """Return True if dependences have been calculated (call initialize)."""
        return self.__isvalid  


class Comp():
    """
    Matrix model of an abstract component, ToF version.
    
    Phase space variables are indexed in following order:
    x, z, incident angle, d_lambda/lambda, time.
    This component represents just a transparent plane. 
    Actual components like slits or collimators should be derived from Comp.
    
    Parameters
    ----------
    distance: float
        distance from the sample [mm]
    
    Notes
    -----
    Always call ``initialize`` after construction, before using any property. 
    
    """
     
    def __init__(self, distance=0.0):
        self.__L = distance 
        self.__isvalid = False
        self.id = 'Comp'
  
    def initialize(self, C, v0): 
        """
        Recalculate dependent data after change of instrument setting.
                  
        Parameters
        ----------
        C: array
            Transport matrix from the sample to the z=0 plane 
            of the preceding component.
        v0: float
            Mean neutron velocity [mm/us].
        """  
        #self.__Cpre = C
        self.__Cin = self.cal_Cin(C, v0)
        self.__Cout = self.cal_Cout(C, v0)
        # NOTE: T must be calculated after Cin, Cout !!
        self.__T = self.cal_T()
        self.__isvalid = True

    def display(self):
        """Print debug information."""        
        prnMat(self.T, id=self.id+".T")        
        
    def cal_T(self):
        """Calculate transmission matrix."""
        return np.zeros((5,5))
    
    def cal_Cin(self, C, v0):
        """Calculate input transport matrix."""
        return CT(self.__L, v0).dot(C)
        
    def cal_Cout(self, C, v0):
        """Calculate output transport matrix."""
        return self.cal_Cin(C, v0)
    
    @property
    def T(self):
        """Transmission matrix."""      
        return self.__T
        
    @property
    def Cin(self):        
        """
        Output transport matrix.
        
        Transport from the sample to the z=0 plane of this component, entry.
        
        """       
        return self.__Cin

    @property    
    def Cout(self):
        """
        Input transport matrix.
        
        Transport from the sample to the z=0 plane of this component, exit.
        
        """        
        return self.__Cout

    @property  
    def isvalid(self):
        """Return True if dependences have been calculated (call initialize)."""
        return self.__isvalid
        

class Slit(Comp):
    """
    Matrix model of a slit, ToF version.
    
    Parameters
    ----------
    distance: float
        distance from the sample to the slit [mm]. 
        Set a negative value if the component is from the sample up-stream.
    width: float
        slit width [mm]
    
    Notes
    -----
    Always call ``initialize`` after construction, before using any property. 
    
    """
    
    def __init__(self, distance=0, width=2.0):
        super().__init__(distance)
        if (distance == 0): 
            self.w = width*g_tri
        else:
            self.w = width*g_uni
        self.id = 'Slit'

    def cal_T(self):
        """Override Comp.cal_T for a Slit."""
        V = self.Cin[0,:]
        n = np.size(V)    
        if (self.w > 0):
            TT = (self.w**-2)*np.kron(V,V).reshape(n,n)
        else:
            TT = super().cal_T()
        return TT


class Guide(Comp):
    """
    Matrix model of a neutron guide, ToF version.
    
    Width represents the constraint on beam size at the guide exit,
    m-value defines the constraint on beam divergence, assuming uniform
    distribution with FWHM = 2*m*GammaNi*wavelength.
    Set ``width = None`` to consider only the divergence constraint.
    Set ``m = None`` to consider only the beam width constraint.
    
    Parameters
    ----------
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
        self.m = m         # m-value
        self.id = 'Guide'      

    def initialize(self, C, v0): 
        """Override ``initialize`` for a Guide."""
        super().initialize(C,v0)
        self.lam0 = hovm/v0
        self.div = 2*self.m*GammaNi*self.lam0*g_uni
        
    def cal_T(self):
        """Override Comp.cal_T for a Guide."""
        TT = super().cal_T()
        V = self.Cin[0,:]
        n = np.size(V)
        if (self.w):
            TT = TT + ((self.w*g_uni)**-2)*np.kron(V,V).reshape(n,n)
        if (self.m):
            V = self.Cin[1,:]
            TT = TT + (self.div**-2)*np.kron(V,V).reshape(n,n)
        return TT


class Pulse(Comp):
    """
    Matrix model of a pulse source, ToF version.
    
    Parameters
    ----------
    distance: float
        distance from the sample to the source [mm]. 
        Set a negative value if the component is from the sample up-stream.    
    tau: float
        Pulse width (FWHM) [us]
    shape: float
        Pulse shape factor. (tau*shape)^2 should yield 2*sigma^2, where sigma^2 
        is the 2nd moment of the pulse distribution.
        If None, a Gaussian shape is assumed. For example:
        
        Use ``shape=sqrt(2)`` if tau == decay time of an exponential pulse.
        
        Use ``shape=1/sqrt(3)`` if tau == FWHM of a triangular pulse.
    
    """
    
    def __init__(self, distance=0, tau=None, shape=None):
        super().__init__(distance)
        if (tau is not None):
            self.tau=tau*shape
        else:
            self.tau=tau*g_norm    # time width [us]
        self.id='Pulse'
   
    def cal_T(self):
        """Override Comp.cal_T for a Pulse."""
        TT = super.cal_T()
        V = self.Cin[3,:]
        n = np.size(V)
        if (self.tau > 0):
            TT = TT + (self.tau**-2)*np.kron(V,V).reshape(n,n)
        return TT


class TofDetector(Slit):
    """
    Matrix model of a detector, ToF version.
    
    Defines spatial and time resolution.
    The amin, amax parameters define the detector unit angular range. 
    [amin, amax] is used by the matrix model to average results over 
    wavelengths corresponding to this angular range for given dhkl value.
    
    Parameters
    ----------
    distance: float
        Distance from the sample to the detector [mm].    
    binwidth: float
        Bin width [mm]
    binshape: float
        Shape factor for the bin. Use g_norm if binwidth represents FWHM of 
        a Gaussian resolution. Use g_uni for a uniform bin of this width. 
    dtau: float
        Time resolution (Gaussian FWHM) [us]
    bins: list
        A list of detector channels data. Use getDetBinsCyl 
        or getDetBinsFlat to generate the list.
    
    """  
              
    def __init__(self, distance=0, binwidth=None, binshape=None, dtau=None, 
                 bins=None):
        super().__init__(distance, binwidth*binshape/g_uni)
        # detector channels
        if (not isinstance(bins, list) or len(bins)<1):
            raise Exception("TofDetector requires that 'bins' is a non-empty list.")
        self.dist = self.__L
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
        """Set current channel: assign bin data for it."""
        self.alpha = self.bins[ic][0]
        self.theta = 0.5*self.alpha
        self.L = self.dist + self.bins[ic][1]
        self.isel = ic
    
    def getChannel(self):
        """Return alpha, theta and cos(phi) for actually selected channel."""
        return [self.alpha, self.theta]


        