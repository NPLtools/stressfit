# -*- coding: utf-8 -*-
"""
Convolution integrals with event lists.

This module defines classes for Monte Carlo convolution of a cloud 
of scattering events with a cell shape, considering material properties such as
beam attenuation or intrinsic distributions of scattering intensity and 
lattice strain. 

Assumptions:
    - Material properties are functions of a single scalar position variable, 
      which is either the position along the scan trajectory, or the depth 
      under a reference surface (defined by the Cell object).
    - Beam attenuation is uniform and isotropic, only wavelength ependence is 
      considered.

Use :class:`GaugeIntegral` to calculate gauge effects such as 
scattering intensity, pseudo-strain, sampling centre of gravity or sampling 
width as a function of scan position.

To be implemented:    
    - Intrinsic scalar strain distribution: strain is independent on deviations
      of the q-vector direction from its mean value. 
    - Intrinsic strain tensor distribution. Requires 1D distribution of strain
      tensor on input. Cases to be considered: 
          a. Constant principal directions and 3 position-dependent principal 
             strain components. 
          b. 6 position-dependent strain cmponents
          c. Planar strain (2 equal lateral and 1 normal strain components)
          d. Axial strain (1 axial and 1 normal component)
    - Intrinsic stress tensor distribution. Like previous, but strain is 
      calculated from compliance tensor provided by user.

Wish list:
    - Texture anisotropy: account for the effect of texture anisotropy on
      beam attenuation and scattering intensity. Develop scripts for MTEX to
      calculate required lookup tables.
    - Deal with multi-phase materials: define composite cells with 2 or more 
      different materials.
    - Deal with ToF data: use sampling events for multiple reflections/materials.
    

Created on Tue Sep 19 10:04:12 2023
@author: Jan Saroun, saroun@ujf.cas.cz
"""
from .cuda import cp
from .cells import Cell
from .events import Sampling
from .scan import ScanGeometry


class GaugeIntegral():
    """Convolution of event lists with sample properties.
    
    Parameters
    ----------
    cell : Cell
        Object describing sample shape and material properties.
    sampling : Sampling
        Object with sampling event list.
    scan : ScanGeometry
        Definition of scan providing coordinate transformations.
    
    """
    
    def __init__(self, cell:Cell, sampling:Sampling, scan:ScanGeometry):
        self._need_init = True # call to initiate is needed
        self._cell = cell
        self._sampling = sampling
        self._scan = scan
        self._nev = sampling.nrec 
        self._ev = None
        self._gauge = {}
 
    
    def _check_scan_variable(self, varid):
        if self._ev is None:
            msg = 'GaugeIntegral not initiated.'
            raise Exception(msg)
        if not varid in ['position', 'depth']:
            msg = 'Only position and depth are allowed '
            msg += 'as scan variable name. {} was given.'
            raise Exception(msg.format(varid))
        if varid == 'depth' and not varid in self._ev:
            msg = 'Cannot use depth as scan variable.' 
            raise Exception(msg)
    
    @property
    def nev(self):
        """Number of events used for integration."""
        return self._nev
    
    def reset(self):
        """Request for new initialization."""
        self._need_init = True 
    
    def initiate(self, nev=5000, random=False):
        """Prepare dependent fields for fast integration procedure.
        
        Select nev events from sampling data, convert to CuPy if required,
        and calculate dependent arrays. 
        """
        nev = min(nev, self._sampling.nrec)
        # get event data
        self._ev = self._sampling.get_events_ex(self._scan, nev=nev, random=random)
        self._nev = nev
        self._cell.set_events(self._ev['r'],self._ev['ki'],self._ev['kf'])
        self._ev['q'] = self._ev['kf']-self._ev['ki']
        self._ev['att'] = self._cell.get_extinction()
        self._ev['position'] = self._scan.to_scan_coord(self._ev['r'])['pos']
        d = self._cell.depth
        if d is not None:
            self._ev['depth'] = d
        q = self._cell.qref
        if q is not None:
            self._ev['qref'] = q
        # initiate variable for gauge data, to be filled by conv_gauge
        self._gauge.clear()
        self._need_init = False

    def conv_gauge(self, wfunc=1.0, varid='position'):
        """Convolution of sampling distribution with the sample.
        
        Parameters
        ----------
        wfunc : obj
            Optional weight function for sampling (e.g. intrinsic scattering
            intensity).
        varid : str
            Variable name - spatial variable which `wfunc` depends on. 
            It can be either `depth` for depth under the reference surface,
            or `position` for position along the scan.
            
        Returns
        -------
        dict
            Integrated scan arrays, such as actual sampling volume 
            sampling centre of mass or pseudo-strain as functions of scan 
            position.
        
        Output arrays
        -------------
        x : array_like
            Scan positions in step uinits (deg or mm) 
        pos : array_like
            Scan positions in step units (deg or mm). 
        err_pos : array_like
            Error of scan position (stdev).      
        width : array_like
            Width (FWHM-equivalent) of the sampling volume. It is measured 
            along the scan (in scan step units).
        cnts : array_like
            Scattering intensity. It is calculated as the sampling volume 
            multiplied by beam attenuation factor and the optional 
            weight function.
        err_cnts : array_like
            Error of `cnts` (stdev).
        cog : array_like
            Centre of gravity of actual sampling distribution given as x,y,z 
            coordinates (sample reference frame). The array shape is (3,:).
        depth : array_like
            The information depth under the reference surface (if defined).
        err_depth : array_like
            Error of `depth` (stdev).
        eps : array_like
            Pseudo-strain
        err_eps : array_like
            Error of `eps` (stdev).
        """  
        self._check_scan_variable(varid)        
        if self._need_init:
            self.initiate(nev=self._nev)
        ns = self._scan.nsteps
        x = self._ev[varid]
        if isinstance(wfunc, (float,int)):
            weight = wfunc
        else:
            weight = wfunc(x)
            
        # event weight
        ps = self._ev['p'].reshape(ns,-1)
        sump = cp.sum(ps, axis=1)
        
        # inside counter
        isin = self._cell.inside        
        isin_s = isin.reshape(ns,-1)
        sumn = cp.sum(isin_s,axis=1)
       
        # total weight
        att = self._ev['att'][0] + self._ev['att'][1]
        w = self._ev['p']*isin*att*weight  
        ws = w.reshape(ns,-1)
        sumw = cp.sum(ws,axis=1)
        
        # sum weights
        yi = sumw
        yi2 =  cp.sum(ws**2,axis=1)
        
        # By this we avoid division-by-zero warnings:
        sg = cp.array((sumn > 0) & (sumw > 0) & (sump > 0), dtype=int)
        sumn += (1-sg)*1
        sump += (1-sg)*1.
        sumw += (1-sg)*1
        
        # save weighting data for later use by other convolution procedures
        self._gauge['w'] = w
        self._gauge['sumw'] = sumw
        self._gauge['sumn'] = sumn
        self._gauge['sg'] = sg
        self._gauge['varid'] = varid

        # centre of gravity
        r = self._ev['r']
        rw = w*r
        rws = rw.reshape(3,ns,-1)
        cog = cp.sum(rws, axis=2)
        
        
        
        
        # position along scan
        wpo = w*self._ev['position']
        wpo2 = w*self._ev['position']**2
        po = wpo.reshape(ns,-1)
        po2 = wpo2.reshape(ns,-1) 
        
        # intensity
        cnts = sg*yi/sump
        cnts2 = sg*yi2/sump
        err_cnt = sg*cp.sqrt(cp.absolute(cnts2 - cnts**2)/sumn)
        
        # cog
        cog = cog*sg/sumw
        
        # position along scan
        pos = sg*cp.sum(po,axis=1)/sumw
        pos2 = sg*cp.sum(po2,axis=1)/sumw
        width = sg*cp.sqrt(cp.absolute(pos2 - pos**2))
        
        # collect output
        out = {}
        out['x'] = self._scan.steps
        out['pos'] = pos
        out['err_pos'] = width/cp.sqrt(sumn)
        out['width'] = 2.3548*width
        out['cnts'] = cnts
        out['err_cnts'] = err_cnt
        out['cog'] = cog
        # depth under reference surface
        if self._ev['depth'] is not None:
            d = self._ev['depth']
            d2 = d**2
            dw = cp.reshape(w*d, (ns,-1))
            d2w = cp.reshape(w*d2, (ns,-1))             
            depth = sg*cp.sum(dw,axis=1)/sumw
            depth2 = sg*cp.sum(d2w,axis=1)/sumw
            err_depth = sg*cp.sqrt(cp.absolute(depth2 - depth**2))/cp.sqrt(sumn)
            out['depth'] = depth
            out['err_depth'] = err_depth
        # pseudo-strain
        if 'eps' in self._ev:
            e = self._ev['eps']
            e2 = e**2
            ew = cp.reshape(w*e, (ns,-1))
            e2w = cp.reshape(w*e2, (ns,-1)) 
            eps = sg*cp.sum(ew,axis=1)/sumw
            eps2 = sg*cp.sum(e2w,axis=1)/sumw
            err = sg*cp.sqrt(cp.absolute(eps2 - eps**2)/sumn) #+ (1-sg)*cp.average(eps)
            out['eps'] = eps
            out['err_eps'] = err
            out['mask'] = sg
            # save in gauge data for future use
            self._gauge['eps'] = eps
            self._gauge['err_eps'] = err
        if 'qref' in self._ev:
            q = self._ev['qref']
            rw = w*q
            rws = rw.reshape(3,ns,-1)
            qref = sg*cp.sum(rws, axis=2)/sumw
            out['qref'] = qref
        return out

    def conv_strain(self, efunc=0.0):
        """Convolution of sampling distribution with scalar strain distribution.
        
        **NOTE**: Previous call of :meth:`conv_gauge` is assumed. 
        
        The strain function (`efunc`) should describe intrinsic strain 
        distribution as a function of scan position (varid='position')
        or depth under reference surface (varid = 'depth'), where
        `varid` is the parameter passed to :meth:`conv_gauge`. The strain 
        distribution is considered uniform in directions normal to the scan. 
        A dependence on the scattering vector is also not considered. 
        This should be a reasonable approximation if sampling is done for 
        a single reflecting plane and a single detector with small spread 
        of scattering vectors.  
        
        Parameters
        ----------
        efunc : obj
            Strain distribution function.
            
        Returns
        -------
        dict
            Integrated scan arrays with results.
        
        Output arrays
        -------------
        x : array_like
            Scan positions in step units (deg or mm) 
        eps : array_like
            Strain (as observed), including pseudo-strain.
        err_eps : array_like
            Error of `eps` (stdev).
        """        
        # check that gauge data have been calculated
        if not 'w' in self._gauge:
            msg = 'Event weigting data not available. '
            msg += 'Call conv_gauge() first'
            raise Exception(msg)
   
        # get gauge data
        w = self._gauge['w'] 
        sumw = self._gauge['sumw']
        sumn = self._gauge['sumn']
        sg = self._gauge['sg']
        varid = self._gauge['varid']    
        ns = self._scan.nsteps
        x = self._ev[varid]
        if isinstance(efunc, (float,int)):
            strain = efunc
        else:
            strain = efunc(x)            

        # collect output
        out = {}
        out['x'] = self._scan.steps
        e2 = strain**2
        ew = cp.reshape(w*strain, (ns,-1))
        e2w = cp.reshape(w*e2, (ns,-1)) 
        eps = sg*cp.sum(ew,axis=1)/sumw
        eps2 = sg*cp.sum(e2w,axis=1)/sumw
        err2 = cp.absolute(eps2 - eps**2)/sumn
        # add pseudo-strain
        if 'eps' in self._gauge:
            out['eps'] = eps + self._gauge['eps'] 
            out['err_eps'] = cp.sqrt(err2 + self._gauge['err_eps'] **2)
        else:
            out['eps'] = eps
            out['err_eps'] = cp.sqrt(err2)
        out['mask'] = sg
        return out
    