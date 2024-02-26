# -*- coding: utf-8 -*-
"""
Handling of sampling event lists.

Includes file input with various formats, random reshufling and getters for 
required event subsets.

Created on Tue Sep 19 10:09:45 2023
@author: Jan Saroun, saroun@ujf.cas.cz
"""

import numpy as np
from .scan import ScanGeometry
from .cuda import cp

class Sampling():
    """Base class for objects describing sampling of scattering events.
    
    Descendants should implement :meth:`import_events` (called internally 
    by the constructor) for handling particular sampling data formats.
    
    Use StrainSampling for diffraction measurements of lattice strain.

    """

    def __init__(self):
        # initiate info - stores calculated sampling properties
        self._info = {}
        # reference to a file 
        self._info['file'] = '' 
        # mean ki, kf and centre of gravity
        for key in ['ki0','kf0','cog']:
            self._info[key] = np.zeros(3)
        # number of events
        self._info['nrec'] = 0 
        # initiate sdata for sampling arrays.
        self.sdata = {}
        # position, initial and final wave vectors, weight, scattering vector
        for key in ['r','ki','kf','p']:
            self.sdata[key] = None

    def import_events(self, source, subset=None, data_format=None, **kwargs):
        """Import event arrays from the source data.
        
        Called by the constructor. Nothing is done by the base Sampling class.
        Descendants should implement import method and store results in 
        the `sdata` and `info` fields.  

        Parameters
        ----------
        source: dict
            Source data given as dictionary.
        subset: int or range
            Subset of events to be taken from data.
        data_format : obj
            Format of source data. Descendants may implement import methods
            for multiple formats.           
        **kwargs : dict
            Other parameters used by the class methods for data import.
        
        """
        msg = 'Not implemented. Subclass must implement abstract methods.'
        raise NotImplementedError(msg) 

    @property
    def info(self):
        """Properties of the event list."""
        return self._info
    
    @property
    def ki(self):
        """Mean incident wave vector."""
        return self._info['ki']

    @property
    def kf(self):
        """Mean final wave vector."""
        return self._info['kf']

    @property
    def cog(self):
        """Centre of gravity."""
        return self._info['cog']
    
    @property
    def nrec(self):
        """Number of events in the list."""
        return self._info['nrec']
   
    def get_events(self, n=None, random=False, items=None):
        """Return a sample of events.
        
        Return event variables from self.sdata.
        
        Parameters
        ----------
        n : int
            Required number of events. Return all available events if None.
        
        random : bool
            Make a random selection of n events. If False, use first n elements 
            from the event list.
        items : list
            List of variable names to export. If non, exports all available 
            sampling variables.
        
        Returns
        -------
        dict
            Arrays for sampling event variables.
        
        """
        out = {}
        
        if items is None:
            keys = list(self.sdata.keys())
        else:
            # always export r,ki,kf,p
            keys = ['r','ki','kf','p']
            for k in items:
                if k in self.sdata and not k in keys:
                    keys.append(k)
        nr = min(n,self.nrec)
        # random selection of n events
        if random:
            idx = np.linspace(0, self.nrec, num=self.nrec, endpoint=False, 
                              dtype=int)
            np.random.shuffle(idx)
            isel = idx[0:nr]
            for key in keys:
                var = self.sdata[key]
                if len(np.shape(var))>1:
                    out[key] = var[:,isel]
                else:
                    out[key] = var[isel]
        # first n events
        else:
            for key in keys:
                var = self.sdata[key]
                if len(np.shape(var))>1:
                    out[key] = var[:,0:nr]
                else:
                    out[key] = var[0:nr]
        return out
   
    
    def get_events_ex(self, scan:ScanGeometry, nev=5000, random=False, 
                      items=None):
        """Get sampling data as arrays expanded by the number of scan points.
        
        Each event is transformed by scan translations and rotations.
        """
        nev = min(nev, self.nrec)
        # get event data
        events = self.get_events(n=nev, random=random, items=items)
        # convert to cupy arrays
        ev = {}
        for k in events:
            ev[k] = cp.asarray(events[k])
        
        # apply global scan transformation to the events
        if scan.transform is not None and not scan.transform.is_identity:
            ev['r'] = scan.transform.r_to_loc(ev['r'])
            ev['ki'] = scan.transform.r_to_loc(ev['ki'])
            ev['kf'] = scan.transform.r_to_loc(ev['kf'])
        # get scan transformations/positions
        scan_pos = scan.positions()
        ns = scan.nsteps
        # get events for each scan position and align all along one dimension
        evx = {}
        for k in events:
            sh = events[k].shape
            if len(sh)==1:
                evx[k] = cp.zeros(ns*sh[0], dtype=float)
            else:
                evx[k] = cp.zeros((sh[0], ns*sh[1]), dtype=float)
        for i in range(ns):
            j1 = i*nev
            j2 = j1+nev
            t = scan_pos['tr'][i]
            # transform global event coordinates to cell root basis
            evx['r'][:,j1:j2] = t.r_to_loc(ev['r'])
            evx['ki'][:,j1:j2] = t.v_to_loc(ev['ki'])
            evx['kf'][:,j1:j2] = t.v_to_loc(ev['kf'])
            for k in filter(lambda el: el not in ['r','ki','kf'], evx):
                sh = events[k].shape
                if len(sh)==1:
                    evx[k][j1:j2] = ev[k]
                else:
                    evx[k][:,j1:j2] = ev[k]
        return evx    
   
    def print_properties(self):
        """Print selected properties of the event list."""
        print('This is the base Sampling class with no data.') 
        print('Use descendants of Sampling to store sampling data.')
    

class StrainSampling(Sampling):
    """Sampling data for diffraction measurements of lattice strain.

    The sampling events can be obtained by MC ray-tracing simulation
    of the instrument at given setup. Generation of such a file
    is implemented in SIMRES ver. 6.3.5 and above 
    (https://github.com/saroun/simres)

    Parameters
    ----------
    source: dict
        Source data given as dictionary (See :meth:`import_events`).
    subset: int or range
        Subset of events to be taken from data.
    data_format : str
        Format identifier for source data. See :meth:`import_events` 
        for details.
    **kwargs : dict
        Other parameters passed to :meth:`_calc_properties`
    
    Example
    -------
    .. highlight:: python
    .. code-block:: python
        
        import numpy as np
        my_array = np.loadtxt('myevents.dat')
        my_data = {'data':my_array,'columns':[1, 4, 7, 10, 11]}
        sampling = StrainSampling(mydata, subset=(0,5000))
        sampling.print_properties()
    
    """
    
    def __init__(self, source, subset=None, data_format='array', **kwargs):
        super().__init__()
        self.import_events(source, subset=subset, data_format=data_format,
                           **kwargs)

    def import_events(self, source, subset=None, data_format='array', **kwargs):
        """Import event arrays from the source data.
        
        Called by the constructor.   

        Parameters
        ----------
        source: dict
            Source data given as dictionary (see below).
        subset: int or range
            Subset of events to be taken from data.
        data_format : str
            Format identifier for source data. Implemented formats are
            
            - array: an array with columns corresponding to scattering
              event variables (r, ki, kf, p, dhkl).
        **kwargs : dict
            Other parameters passed to :meth:`_calc_properties`

        Source formats:
        ---------------

        array:
        
            `source` includes following items:
                
            data : array_like
                An array with columns corresponding to scattering
                event variables (r, ki, kf, p, dhkl).
            columns : list of int
                Column indices corresponding to [r, ki, kf, p, dhkl]

        """
        # extract sampling data
        if data_format=='array':
            # verify input
            for key in ['data','columns']:
                if not key in source:
                    raise Exception('Missing item in source data: {}'.format(key))
            self._import_array(subset=subset, **source)
        else:
            raise Exception('Unknown source format: {}'.format(data_format))
        
        self._info['nrec'] = self.sdata['r'].shape[1]
        # calculate integral sampling properties
        self._calc_properties(**kwargs)
        # add strain values to event list
        d0 = self._info['d0']
        self.sdata['eps'] = (self.sdata['dhkl'] - d0)/d0

    def _import_array(self, data=None, columns=[1, 4, 7, 10, 11], subset=None):
        """Extract given subset of sampling events from source data.
        
        Parameters
        ----------
        data : array_like
            An array with columns corresponding to scattering
            event variables (r, ki, kf, p, dhkl).
        columns : list of int
            Column indices corresponding to [r, ki, kf, p, dhkl]
        subset: int or range
            Subset of events to be taken from data.
        """
        q1 = data.shape[1] > np.max(columns)
        q2 = len(columns)>=5
        if not q1 and q2:
            raise Exception('Incorrect column index list: {}.'.format(columns))
        
         # get named list of columns
        keys = ['r','ki','kf','p','dhkl']
        # column span for each variable
        L = [3,3,3,1,1] 
        self.ids = {keys[i]:columns[i] for i in range(len(keys))}
        # subset range of events
        if subset is not None:
            if isinstance(subset, int):
                rang = slice(0, subset)
            else:
                rang = slice(*subset)
        else:
            rang = slice(0, -1)
        # store variables in dict
        self.sdata = {}
        for i in range(len(keys)):
            key = keys[i]
            ic = self.ids[key]
            # Transpose data ! 
            self.sdata[key] = data[rang, ic:ic+L[i]].T       
        # normalize weights
        psum = np.sum(self.sdata['p'])
        self.sdata['p'] = self.sdata['p']/psum
    
    def _calc_properties(self, **kwargs):
        """Calculate integral properties of the event list.
        
        See the Parameters list for the calculated quantities.
        In addition, add calculated arrays to sdata: strain (eps) 
        and scattering vector (q). 
        
        Parameters
        ----------
        Any of the following properties, which then overwrite the calculated 
        ones:
            
        cog : array_like
            Distribution centre of gravity in mm
        ki : array_like
            Mean incident wave vector in A^-1
        kf : array_like
            Mean final wave vector in A^-1
        q : array_like
            Mean scattering vector in A^-1
        width : array_like
            Distribution widths along x,y,z [mm]
        d0 : float
            Mean dhkl in A^-1
        wav : float
            Mean wavelength in A
        tth : float
            Mean scattering angle in deg
        """
        out = {}
        # Calculate centre of mass of the distribution 
        p = self.sdata['p']
        # centre of mass
        out['cog'] = np.sum(p*self.sdata['r'], axis=1)
        # mean ki, kf
        out['ki'] = np.sum(p*self.sdata['ki'], axis=1)
        out['kf'] = np.sum(p*self.sdata['kf'], axis=1)
        # resolution widths along x,y,z
        r2 = self.sdata['r']**2
        std = np.sum(p*r2, axis=1) - out['cog']**2
        out['width'] = np.sqrt(std)*np.sqrt(8*np.log(2))
        # mean dhkl
        out['d0'] = np.sum(p*self.sdata['dhkl'])
        # mean wavelength
        k0 = np.linalg.norm(out['ki'],axis=0)
        out['wav'] = 2*np.pi/k0
        # mean scattering angle
        out['tth'] = np.arcsin(out['wav']/2/out['d0'] )*360/np.pi
        # override by values defined in kwargs, convert to arrays or floats
        for key in out:
            if key in kwargs:
                if np.size(out[key])>1:
                    out[key] = np.array(kwargs[key], dtype=float)
                else:
                    out[key] = float(kwargs[key])
        
        self._info.update(out)
   
    def print_properties(self):
        """Print selected properties of the event list."""
        print('Number of loaded sampling points: {:d}'.format(self.nrec)) 
        print('Gauge centre: [{:g}, {:g}, {:g}] '.format(*self._info['cog']))
        print('Mean wavelength: {:g}'.format(self._info['wav']))
        print('2 theta: {:g}'.format(self._info['tth']))
        print('d0: {:g}\n'.format(self._info['d0'])) 

