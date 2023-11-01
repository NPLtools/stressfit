# -*- coding: utf-8 -*-
"""
Convolution integrals with event lists.

Each integration is over inner area of a Surface object and includes 
attenuation factors for input and output rays. 

This module provides methods for following MC integrations:
    - Scalar function (e.g. intrinsic scattering intensity) 
    - Scalar strain distribution (1D). Assumes constant sampled strain at given 
      depth (relative to reference surface). Returns measured strain (scalar).
    - Tensor strain distribution (1D). Assumes strain tensor is constant at 
      given depth (relative to reference surface). Principal directions do
      not change. Returns up to 3 strain components for principal directions.
      Less components are needed for planar strain option. 
    - Tensor stress distribution (1D). Like previous, but strain is calculated
      from compliance tensor provided by user.

Created on Tue Sep 19 10:04:12 2023
@author: Jan Saroun, saroun@ujf.cas.cz
"""
from .cuda import cp
from .cells import Cell
from .events import Sampling
from .scan import ScanGeometry


class GaugeIntegral():
    """Class encapsulating data and methods for fast convolution integrals.
    
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
        self._func = None
    
    
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
        self._ev['position'] = self.cal_scan_positions()
        d = self._cell.cal_depth()
        if d is not None:
            self._ev['depth'] = d
        self._need_init = False

    def cal_scan_positions(self):
        """Calculate event positions along scan."""
        pos = self._scan.to_scan_coord(self._ev['r'])['pos']
        return pos 

    def conv_func(self, func, varid='depth'):
        """Convolution with scalar funcion.
        
        Parameters
        ----------
        func : obj
            Weight function to be used in convolution with the sampling 
            distribution.
        varid : str
            Variable name - spatial variable which func depends on. 
            It can be either `depth` for depth under reference surface,
            or `position` for position along scan.
            
        """    
        
        if not varid in ['position', 'depth']:
            msg = 'Only position and depth are allowed '
            msg += 'as scan variable name. {} was given.'
            raise Exception(msg.format(varid))
        if varid == 'depth' and not varid in self._ev:
            msg = 'Cannot use depth as scan variable.' 
            raise Exception(msg)
        
        if self._need_init:
            self.initiate(nev=self._nev)
        ns = self._scan.nsteps
        
        
        
        x = self._ev[varid]
        if func is None:
            fval = 1.0
        elif isinstance(func, float):
            fval = func
        else:
            fval = func(x)
            
        # event weight
        ps = self._ev['p'].reshape(ns,-1)
        sump = cp.sum(ps, axis=1)
        
        # inside counter
        isin = self._cell.inside        
        isin_s = isin.reshape(ns,-1)
        sumn = cp.sum(isin_s,axis=1)
       
        # total weight
        att = self._ev['att'][0] + self._ev['att'][1]
        w = self._ev['p']*isin*att*fval  
        ws = w.reshape(ns,-1)
        sumw = cp.sum(ws,axis=1)
        
        # sum weights
        yi = sumw
        yi2 =  cp.sum(ws**2,axis=1)
        
        # This will avoid division by zero warnings:
        sg = cp.array((sumn > 0) & (sumw > 0) & (sump > 0), dtype=int)
        sumn += (1-sg)*1
        sump += (1-sg)*1.
        sumw += (1-sg)*1

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
        return out

