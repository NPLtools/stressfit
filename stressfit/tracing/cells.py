# -*- coding: utf-8 -*-
"""
Defines cells with material properties.

Created on Fri Sep 15 17:08:44 2023
@author: Jan Saroun, saroun@ujf.cas.cz
"""

import numpy as np
from .cuda import cp, asnumpy
from .primitives import Surface, Group


class Extinction():
    """Beam extinction function.
    
    Defines wavelength dependent beam extinction as a table or 
    polynomial function.

    Depending on the value of `kind` parameter, it can be either scalar,
    polynomial, or a lookup table. Corresponding `**kwargs` arguments are then:
    
    Extinction.SCALAR:
        mu : float
            Extinction cross-section in [1/cm]
    Extinction.POLY:
        p : list(3)
            Coefficients in [1/cm, 1/cm/A, 1/cm/A**2]
    Extinction.TABLE:
       table : array_like
           Array of shape (2,:) with wavelengths and mu in [A, 1/cm].

    Parameters
    ----------
    kind : int
        One of Extinction constants (SCALAR,POLY,TABLE)
    **kwargs :
        Arguments depending on `kind`
    
    Example
    -------
    
    .. highlight:: python
    .. code-block:: python
        
        # scalar
        ex = Extinction(Extinction.SCALAR, mu=1.1)
        # polynomial
        ex = Extinction(Extinction.POLY, p=[0.5, 0.5, 0.0])
        # table
        data = np.loadtext(filename)
        ex = Extinction(Extinction.TABE, table=data)
        
    """
    
    SCALAR = 0
    POLY = 1
    TABLE = 2
    
    def __init__(self, kind=0, **kwargs):
        if kind == Extinction.SCALAR:
            self.define_scalar(**kwargs)
        elif kind == Extinction.POLY:
            self.define_poly(**kwargs)
        elif kind == Extinction.TABLE:
            self.define_table(**kwargs)
        else:
            raise Exception('Wrong kind argument. Use SCALAR, POLY or TABLE')

    def define_scalar(self, mu=1.15):
        """Set beam extinction cross-section as a scalar value.
        
        Parameters
        ----------
        mu : float
            Extinction cross-section in [1/cm]
        """
        self._mu = 0.1*mu # convert to 1/mm
        self._coeff = None
        self._table = None
        self._ifn = Extinction.SCALAR   # use scalar value
    
    def define_poly(self, p=[0.5, 0.5, 0.0]):
        """Set beam extinction cross-section as quadratic function.
        
        mu = p[0] + p[1]*lambda + p[2]*lambda**2
        
        Parameters
        ----------
        p : list(3)
            Coefficients in [1/cm, 1/cm/A, 1/cm/A**2]
        """
        self._mu = None
        self._coeff = np.multiply(0.1, p) # convert to 1/mm
        self._table = None
        self._ifn = Extinction.POLY   # use polynomial

    def define_table(self, table=None):
        """Set beam extinction cross-section as quadratic function.
        
        mu = p[0] + p[1]*lambda + p[2]*lambda**2
        
        Parameters
        ----------
        table : array_like
            Array of shape (2,:) with wavelengths and mu in [A, 1/cm].
        """
        msg = 'Extinction table format error: \n'
        msg += 'Two-column array with lambda[A] and mu[1/cm] values is expected.'
        if table is None:
            raise Exception(msg)
        sh = cp.shape(table)
        if sh[0] == 2:
            self._table = cp.asarray([table[0,:], 0.1*table[1,:]])
        elif sh[1] == 2:
            self._table = cp.asarray([table[:,0], 0.1*table[:,1]])
        else:
            raise Exception(msg)
        self._mu = None
        self._coeff = None
        self._ifn = Extinction.TABLE   # use table

    def get_mu(self, wavelength):
        """Calculate attenuation coefficient for given wavelengths.
        
        Parameters
        ----------
        wavelength : float
            A scalar or array of wavelength values in [A].
        """
        if self._ifn==Extinction.SCALAR:
            mu = self._mu*wavelength
        elif self._ifn==Extinction.POLY:
            mu = self._p[0] + self._p[1]*wavelength + self._p[2]*wavelength**2
        elif self._ifn==Extinction.TABLE:
            mu = cp.interp(wavelength,
                             self._table[0,:], self._table[1,:],
                             left=self._table[1,0], right=self._table[1,-1])
        return mu    


class Material():
    """Sample material.
    
    Describe material properties, such as beam attenuation cross-section.
    
    """
    
    def __init__(self):
        self._ext = Extinction()
        
    def set_extinction(self, obj, **kwargs):
        """Define extinction.
        
        Beam extinction is defined by the class :class:`Exctinction`.
        
        Parameters
        ----------
        obj : obj
            If `obj` is an instance of :class:`Exctinction`, replace 
            the current extinction object with `obj`. Otherwise, 
            create it as Exctinction(kind=obj, **kwargs).
        """
        if isinstance(obj, Extinction):
            self._ext = obj
        else:
            self._ext = Extinction(kind=obj, **kwargs)  

    def xatt(self, wavelength):
        """Get beam attenuation cross-section [1/mm] for given wavelength(s).
        
        Parameters
        ----------
        wavelength : float
            A scalar or array of wavelength values in [A].
        """
        return self._ext.get_mu(wavelength)    


class Cell():
    """Sample cell.
    
    Parameters
    ----------
    surface : :class:`primitives.Surface`
        Surface object describing sample shape and tracing methods.
    material : :class:`Material`
        material object, describes cell material properties.
    
    """
    
    def __init__(self, surface:Surface, material:Material):
        if isinstance(surface, Group):
            self._surf = surface
        else:
            self._surf = Group()
            self._surf.add(surface)
            surface.isref=True
        self._mat = material
        self._refs = None # reference surface for depth definition
        self._scan = {'orig': [0,0,0], 'dir': [0,0,1]}
    
    
    @property
    def extinction(self):
        """Extinction object, an instance of :class:`Extinction`."""
        return self._mat._ext
    
    @property
    def material(self):
        """Cell material as an instance of :class:`Material`."""
        return self._mat
    
    @property
    def wavelengths(self):
        """Wavelengths for incident and output rays."""
        return self._wav
   
    @property
    def inside(self):
        """Get 1 for all events inside the cell, otherwise 0."""
        return self._surf.is_in
    
    @property
    def reference(self):
        """Get/Set reference surface object.
        
        Reference surface must be a primitive surface, not a Group, and 
        it must be a member of the cell's surface group. 
        
        Setting is done by an object reference, or an index pointing
        to the root surface group of this cell.
        
        """
        return self._refs

    @reference.setter
    def reference(self, value):
        # set by index from the root surface group
        if isinstance(value,int):
            if abs(value)<len(self._surf.surfaces):
                self._refs = self._surf.surfaces[abs(value)]
            else:
                raise Exception('No surface with index {}'.format(value))
        # set by Surface instance
        elif isinstance(value, Surface) and not isinstance(value, Group):
            if not self._surf.is_member(value):
                raise Exception('Reference surface does not belong to the group')
            else:
                 self._refs = value         
        elif value==None:
            self._refs = value 
        else:
            raise Exception('Type {} cannot be a reference surface'.format(type(value)))
        
    
    def set_events(self, r, ki, kf):
        """Set positions and wave vectors for all events to be processed.
        
        Parameters
        ----------
        r : array_like
            Positions passed as (3,:) arrays.
        ki, kf : array_like
            Input and output wave vectors passed as (3,:) arrays.
        p : array_like
            Weights of the events.
        """        
        self._surf.set_events(r, ki, kf)
        self._surf.evaluate()
        ki0 = cp.linalg.norm(ki, axis=0)
        kf0 = cp.linalg.norm(kf, axis=0)
        self.path_in = ki0*self._surf.time_in
        self.path_out = kf0*self._surf.time_out
        wi = 2*cp.pi/ki0
        wf = 2*cp.pi/kf0
        self._wav = [wi, wf]

    def cal_depth(self):
        """Calculate depths under reference surface."""
        depth = None
        if self._refs is not None:
            depth = self._refs.cal_depth()
        return depth

    def cal_qref(self):
        """Calculate q-vectors in reference surface  coordinates."""
        qref = None
        if self._refs is not None:
            qref = self._refs.cal_qref()
        return qref
    
    def get_extinction(self):
        """Return extinction factors for input and output beams."""
        e1 = self.path_in*self._mat.xatt(self._wav[0])
        e2 = self.path_out*self._mat.xatt(self._wav[1])
        return [cp.exp(-e1), cp.exp(-e2)]
        
    
    
        