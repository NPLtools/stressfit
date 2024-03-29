# -*- coding: utf-8 -*-
"""
Module for convolution of an interpolated function
with sampling distribution defined by an array of scattering events.

The interpolated function can define an intrinsic distribution of scattering
probability or distribution of stress.

Assumes 1D function defined by a set of [depth, p(depth)] points
interpolated by splines of given order.
The depth is related to the actual model of the sample shape, which is
an instance of SampleAbstract class.
The depths should be always defined as normal distance from the main surface
(typically the outer or top surface of the shape).
See the getLocalDepth function of given shape model.

Created on Tue Aug 15 19:12:50 2017
@author: Jan Saroun, saroun@ujf.cas.cz
"""
import numpy as np
from numpy.linalg import norm

shape = None
ext_coeff = 0.115  # absorption coefficient in 1/mm
ext_lambda = None
ext_mu = None

_sampling = None

def _verify_keys(keys, data):
    res = all (k in data for k in keys)
    return res

class Sampling():
    """Set array of sampling events.

    The sampling events can be obtained by MC ray-tracing simulation
    of the instrument at given setup. Generation of such a file
    is implemented in SIMRES ver. 6.3.5 and above

    Parameters
    ----------
        src: dict
            Source data given as dictionary:
                - data: replaces events parameter
                - columns: replaces columns parameter
                - ctr: replaces ctr parameter
                - ki: mean incident k-vector
                - kf: mean final k-vector
        events: array[:,:]
            list of sampling points coordinates and weights
        columns: array of int
            column indices for r, ki, kf, p, dhkl
                - r : event position
                - ki, k : initial and final wave vector
                - p : probability
                - dhkl : dhkl sampled by this event
        ctr: array(3)
            Gauge centre (default is calculated from sampling positions)
       
    """
    
    def __init__(self, src=None, events=None, columns=None, ctr=None):       
        self.sdata = events
        self.idata = columns
        self.file = ''
        if src is not None:
            if not _verify_keys(('data','columns'), src):
                raise Exception('Required fields not found in src.')
            self.src = src
            self.sdata = src['data']
            self.idata = src['columns']
            if 'file' in src:
                self.file = src['file']
        else:
            self.src = {}
        self.src.update(self.calc_properties())
        if (ctr is not None):
            self.sctr = np.array(ctr)
        else:
            self.sctr = self.src['ctr'] 
        self.ki =self.src['ki']
        self.kf =self.src['kf']
        if self.sdata is not None:
            self.setRange(self.sdata.shape[0])

    @classmethod
    def fromdict(cls, source, ctr=None):
        """Create from dictionary.
        
        Requires at least following keys in source:
            - data: coordinates of the sampling events
            - columns: column positions for r, ki, kf, p, dhkl
        
        Parameters
        ----------
        source: dict
            Dictionary with sampling data.
        ctr: list(3)
            gauge centre coordinates (optional)
        
        """
        if 'data' in source and 'columns' in source:
            s = cls(src=source, ctr=ctr)
        else:
            raise Exception('Required fields not found in source.')
        return s
    
    
    def calc_properties(self):
        out = {}
        nrec = self.sdata.shape[0]
        data = self.sdata
        columns = self.idata
        # Calculate centre of mass of the distribution 
        P = data[:nrec,columns[3]]/np.sum(data[:nrec,columns[3]])
        ctr = np.zeros(3)
        w = np.zeros(3)
        ki = np.zeros(3)
        kf = np.zeros(3)
        for i in range(3):
            c =  data[:nrec, columns[0] + i]
            c2 = c*c
            ctr[i] = c.dot(P)
            w[i] = np.sqrt(c2.dot(P)-ctr[i]**2)*np.sqrt(8*np.log(2))
            ki[i] = data[:nrec, columns[1] + i].dot(P)
            kf[i] = data[:nrec, columns[2] + i].dot(P)
        
        wav = 2*np.pi/np.sqrt(ki.dot(ki))
        dmean = data[:nrec, columns[4]].dot(P)
        tth = np.arcsin(wav/2/dmean)*360/np.pi   
        out['data'] = data[:nrec,:]
        out['columns'] = columns
        out['nrec'] = nrec
        out['ctr'] = ctr
        out['width'] = w
        out['dmean'] = dmean
        out['wav'] = wav
        out['tth'] = tth
        out['ki'] = ki
        out['kf'] = kf    
        return out
    
    def print_properties(self):
        nrec = self.sdata.shape[0]
        print('Number of loaded sampling points: {:d}'.format(nrec)) 
        print('Gauge centre: [{:g}, {:g}, {:g}] '.format(*self.src['ctr']))
        print('Mean wavelength: {:g}'.format(self.src['wav']))
        print('2 theta: {:g}'.format(self.src['tth']))
        print('d0: {:g}\n'.format(self.src['dmean'])) 
    
    def setRange(self, nev):
        """Define subset of events and calclate mean dhkl."""
        [jr, jki, jkf, jp, jd] = self.idata[0:5]
        self.nev = min(nev, self.src['nrec'])
        sumd = 0
        sump = 0
        for ir in range(self.nev):
            p = self.sdata[ir, jp]
            dhkl = self.sdata[ir, jd]
            sump += p
            sumd += p*dhkl
        self.d0 = sumd/sump  
        self.sump = sump

def getSampling(nev=0):
    """Return Sampling object.
    
    Parameters
    ----------
    nev: int
        Number if sampling points to use (optional).
    """
    if nev>0:
        _sampling.setRange(nev)
    return _sampling

def setSampling(sampling):
    """Assign sampling for use by convolution methods. 
    
    Parameters
    ----------
    sampling: :obj:`stressfit.sample.Sampling`
        Object with sampling points.

    """
    global _sampling
    _sampling = sampling

def setExtinction(mu = None, table=None):
    """Set beam attenuation.
    
    Set either as scalar (mu) in 1/cm, or as a lookup table
    with mu [1/cm] as a function of wavelength [AA].
    
    """
    global ext_coeff, ext_lambda, ext_mu
    if (mu is None) & (table is None):
        raise RuntimeError('Define one of the arguments: mu or table.')
    if (mu is not None) & (table is not None):
        raise RuntimeError('Define only one of the arguments.')
    
    if (mu is not None):
        ext_coeff = 0.1*mu*np.array([1., 1.])  # convert to 1/mm
        ext_lambda = None
        ext_mu = None
    elif (table is not None):
        if (table.shape[1] < 2):
            msg = 'Extinction table format error: \n'
            msg += 'Two-column array with lambda and mu values is expected.'
            raise RuntimeError(msg)
        ext_lambda = table[:, 0]
        ext_mu = 0.1*table[:, 1]  # convert to 1/mm


def setSamplingEvents(events, columns, ctr = None):
    """ Set array of sampling events

    The sampling events can be obtained by MC ray-tracing simulation
    of the instrument at given setup. Generation of such a file
    is implemented in SIMRES ver. 6.3.5 and above

    Arguments:
        events -- array[:,:]
        columns -- defines column indices for r, ki, kf, p, dhkl
                    r -- event position
                    ki, kf -- initial and final wave vector
                    p -- probability
                    dhkl -- dhkl sampled by this event
    """
    global _sampling
    _sampling = Sampling(events, columns, ctr=ctr)


def getMu(wavelength):
    """ Calculate attenuation coefficient for given wavelength

    """
    x = np.array([1., 1.])
    if (ext_mu is not None):
        mu = x*np.interp(wavelength,
                         ext_lambda, ext_mu, left=ext_mu[0], right=ext_mu[-1])
    else:
        mu = ext_coeff*x
    return mu

def getExtinction(r, ki, kf):
    """ Calculate extinction factor for given event positions and ki,kf vectors

        Using approximation to small radii for higher speed
    """
    paths = shape.rayPaths(r, ki, kf)
    # path =  [depthi, depthf, pathi, pathf]
    p = np.abs(np.array([paths[0], paths[1]]))
    mu = getMu(2.*np.pi/norm(ki))
    arg = -p.T.dot(mu.T)
    res = np.exp(arg)
    return res


def getExtinction2(r, ki, kf):
    """ Calculate extinction factor for given event positions and ki,kf vectors

        Exact, but slower solution for flight paths.
    """
    ext = shape.pathInOut(r, ki, kf)  # get [inside, pin, pout]
    aex = np.array(ext)
    ins = aex[:, 0]
    p = aex[:, 1:3]
    mu = getMu(2.*np.pi/norm(ki))
    arg = -p.dot(mu.T)
    res = ins*np.exp(arg)
    return res

def shuffleEvents():
    """get random sample of nev events from sdata"""
    idx = np.arange(_sampling.sdata.shape[0])
    np.random.shuffle(idx)
    sel = _sampling.sdata[idx,:]
    _sampling.sdata = sel
    if _sampling.nev:
        _sampling.setRange(_sampling.nev)
    else:
        _sampling.setRange(_sampling.sdata.shape[0])


def convResol(x, xdir, iext, nev, ifunc=None):
    """Calculate resolution parameters for given scan.

    Parameters
    ----------
    x :  array
        scan positions [mm] along xdir
    xdir : array
        scan direction in local coordinates
    iext :  int
        model for extinction
        0 = small radius approximation (fast)
        1 = exact (slow)
    nev :  int
        number of events from sdata to integrate
    ifunc : function
        A function which returns scattering probability distribution
        as a function of the information depth.
    
    Returns
    -------
    dict:
        keys = [x, pos, width, ctr, cov]
    
    Where:
        - x : scan positions [mm]
        - pos : information depth (depends on sample shape definition)
        - width : information width (FWHM of depth)
        - ctr :  sampling centre (x,y,z) in local coordinates
        - cov : covariance matrix of sampling distribution 
          (6 elements in Voigt notation)
    """
    global shape, _sampling
    # Voigt indexing
    iv = [[0,0], [1,1], [2,2], [1,2], [0,2], [0,1]]
    # initialize local variables
    nx = x.shape[0]
    yd = np.zeros(nx) # depth
    yd2 = np.zeros(nx) # depth^2
    rc = np.zeros((nx,3)) # centre
    cov = np.zeros((nx,6)) # covariance
    xs = np.matrix(x)
    # where to find r, ki, kf, p, dhkl
    if (_sampling.nev != nev):
        _sampling.setRange(nev)
    rnd = range(0, _sampling.nev)
    [jr, jki, jkf, jp, jd] = _sampling.idata[0:5]
    # scan positions relative to r
    dr = np.array(xs.T*xdir)

    # loop through sampling events to make convolution
    sumw = 0.
    sump = 0.
    for ir in rnd:
        rn = _sampling.sdata[ir, jr:jr+3] - _sampling.sctr
        r0 = shape.getLocalPos(rn)
        ki = shape.getLocalDir(_sampling.sdata[ir, jki:jki+3])
        kf = shape.getLocalDir(_sampling.sdata[ir, jkf:jkf+3])
        p = _sampling.sdata[ir, jp]
        sump += p
        r = r0 + dr
        # get event depths and flags for inside/outside events
        [d, d2, ins] = shape.depthLocal(r)
        # calculate extinction factor
        if (iext == 1):
            ex = getExtinction2(r, ki, kf)
        else:
            ex = getExtinction(r, ki, kf)
        # total weight including extinction and scattering probability
        if ifunc is None:
            I0 = np.ones(nx)
        else:
            I0 = ifunc(d)
        w = p*ex*ins*I0
        sumw += w
        wd = w*d
        yd += wd
        yd2 += wd*d
        for j in range(3):
            rc[:,j] += w*r[:,j]
        for j in range(6):
            [j1,j2] = iv[j]
            cov[:,j] += w*r[:,j1]*r[:,j2]
   
    # This will avoid division by zero warnings:
    sg = np.array( (sumw > 0) , dtype=int)
    sumw += (1-sg)*1
    pos = sg*yd/sumw
    epos = sg*yd2/sumw
    width = 2.3548*sg*np.sqrt(np.absolute(epos - pos**2))
    ctr = np.zeros((nx,3))
    covar = np.zeros((nx,6))
    for j in range(3):
        ctr[:,j] = sg*rc[:,j]/sumw
    for j in range(6):
        cc = sg*cov[:,j]/sumw
        [j1,j2] = iv[j]
        covar[:,j] = cc - ctr[:,j1]*ctr[:,j2]
    res = {}
    res['x'] = x.reshape((nx,1))
    res['pos'] = pos.reshape((nx,1))
    res['width'] = width.reshape((nx,1))
    res['ctr'] = ctr
    res['cov'] =  covar
    return res

def convGauge(x, xdir, iext, nev, ifunc=None):
    """Calculate sampling volume center and pseudo strains.

    Parameters
    ----------
    x :  array
        scan positions [mm] along xdir
    xdir : array
        scan direction in local coordinates
    iext :  int
        model for extinction
        0 = small radius approximation (fast)
        1 = exact (slow)
    nev :  int
        number of events from sdata to integrate
    ifunc : function
        A function which returns scattering probability distribution
        as a function of the information depth.
    
    Returns
    -------
    array:
        [x, pos, width, cnts, eps, ctr, err, err_cnt]
    
    Where:
       -  x : scan positions [mm]
       -  pos : information depth
       -  width : information width (FWHM) 
       -  cnts : intensity (=sampling volume)
       -  eps :  pseudo-strains [1e-6] 
       -  ctr : xyz coordinates of the sampling centre of gravity
       -  err : errors of pseudo-strain [1e-6] 
       -  err_cnt : error of intensity
       
    """
    global shape, _sampling
    # initialize local variables
    nx = x.shape[0]
    yc = np.zeros(nx)
    yc2 = np.zeros(nx)
    yd = np.zeros(nx)
    yi = np.zeros(nx)
    yi2 = np.zeros(nx)
    yd2 = np.zeros(nx)
    rc = np.zeros((nx,3)) # centre
    xs = np.matrix(x)
    ex = np.zeros(nx)
    # where to find r, ki, kf, p, dhkl
    if (_sampling.nev != nev):
        _sampling.setRange(nev)
    rnd = range(0, _sampling.nev)
    [jr, jki, jkf, jp, jd] = _sampling.idata[0:5]
    d0 = _sampling.d0
    # scan positions relative to r
    dr = np.array(xs.T*xdir)

    # loop through sampling events to make convolution
    sumw = 0.
    sump = 0.
    sumn = 0
    for ir in rnd:
        rn = _sampling.sdata[ir, jr:jr+3] - _sampling.sctr
        r0 = shape.getLocalPos(rn)
        ki = shape.getLocalDir(_sampling.sdata[ir, jki:jki+3])
        kf = shape.getLocalDir(_sampling.sdata[ir, jkf:jkf+3])
        p = _sampling.sdata[ir, jp]
        sump += p
        # pseudo-strain at the samplig point
        deps = (_sampling.sdata[ir, jd] - d0)/d0
        r = r0 + dr
        # get event depths and flags for inside/outside events
        [d, d2, ins] = shape.depthLocal(r)
        # calculate extinction factor
        if (iext == 1):
            ex = getExtinction2(r, ki, kf)
        else:
            ex = getExtinction(r, ki, kf)
        # total weight including extinction and scattering probability
        eps = deps
        if ifunc is None:
            I0 = np.ones(nx)
        else:
            I0 = ifunc(d)
        w = p*ex*ins*I0
        wd = w*d
        weps = w*eps
        sumw += w
        sumn += ins
        yi += w
        yi2 += w**2
        yc += weps
        yc2 += weps*eps
        yd += wd
        yd2 += wd*d
        
        #if ir/30 - int(ir/30)  == 0:
        #    fmt = 6*'{:g}\t'
        #    print(fmt.format(w[0], eps*1e6, weps[0]*1e6, sumw[0], 
        #                     yc[0]*1e6, yc[0]/sumw[0]*1e6))
        
        for j in range(3):
            rc[:,j] += w*r[:,j]
    # This will avoid division by zero warnings:
        
    sg = np.array((sumn > 0) & (sumw > 0) & (sump > 0), dtype=int)
    sumn += (1-sg)*1
    sump += (1-sg)*1.
    sumw += (1-sg)*1
    cnts = sg*yi/sump
    cnts2 = sg*yi2/sump
    eps = sg*yc/sumw
    err = sg*yc2/sumw
    pos = sg*yd/sumw
    epos = sg*yd2/sumw
    ctr = np.zeros((nx,3))
    for j in range(3):
        ctr[:,j] = sg*rc[:,j]/sumw
    width = 2.3548*sg*np.sqrt(np.absolute(epos - pos**2))
    err = sg*np.sqrt(np.absolute(err - eps**2)/sumn) #+ (1-sg)*np.average(eps)
    err_cnt = sg*np.sqrt(np.absolute(cnts2 - cnts**2)/sumn)
    return np.array([x, pos, width, cnts, 1e6*eps, 
                     ctr[:,0], ctr[:,1], ctr[:,2], 1e6*err, err_cnt]).T

def convIntensity(x, model, iext=0):
    """Convolution of the scattering intensity distribution

    Arguments:
        x -- scan positions [mm] along xdir
        model -- MCCfit model for intensity
        iext -- model for extinction
                0 = small radius approximation (fast)
                1 = exact (slow)
    Returns
        [I, err, pos, epos, width] = convoluted intensity, error, info depth, error, info width
    """
    global shape, _sampling
    xdir = model.xdir
    nev = model.nev
    # initialize local variables
    nx = x.shape[0]
    yc = np.zeros(nx)
    yc2 = np.zeros(nx)
    yd = np.zeros(nx)
    yd2 = np.zeros(nx)
    xs = np.matrix(x)
    ex = np.zeros(nx)
    if (_sampling.nev != nev):
        _sampling.setRange(nev)
    rnd = range(0, _sampling.nev)
    # where to find r, ki, kf, p, dhkl
    [jr, jki, jkf, jp, jd] = _sampling.idata[0:5]
    # scan positions relative to r
    dr = np.array(xs.T*xdir)
 
    # loop through sampling events
    sump = 0.
    sumn = 0.
    for ir in rnd:
        rn = _sampling.sdata[ir, jr:jr+3] - _sampling.sctr
        r0 = shape.getLocalPos(rn)
        ki = shape.getLocalDir(_sampling.sdata[ir, jki:jki+3])
        kf = shape.getLocalDir(_sampling.sdata[ir, jkf:jkf+3])
        p = _sampling.sdata[ir, jp]
        sump += p
        r = r0.T + dr
        # get event depths and flags for inside/outside events
        [d, d2, ins] = shape.depthLocal(r)
        # calculate extinction factor
        if (iext == 1):
            ex = getExtinction2(r, ki, kf)
        else:
            ex = getExtinction(r, ki, kf)
        # total intensity including extinction and scattering probability
        I = ex*ins*model.intFnc(d)
        yc += p*I
        yc2 += p*I**2
        yd += p*I*d
        yd2 += p*I*d**2
        sumn += ins
    # This will avoid division by zero warnings:
    sg = np.array((yc > 0.) & (sumn > 0), dtype=int)
    yc += (1-sg)*1.
    sumn += (1-sg)*1
    # normalize sums
    val = sg*yc/sump
    err = sg*yc2/sump    
    pos = sg*yd/yc
    derr = sg*yd2/yc
    # Calculate sampling position, width, intensity and error
    width = sg*np.sqrt(np.absolute(derr - pos**2))
    err = sg*np.sqrt(np.absolute(err - val**2)/sumn) # + (1-sg)*np.average(val)  
    epos = width/np.sqrt(sumn) # + (1-sg)*np.average(pos)  
    return [val, err, pos, epos, 2.3548*width]

def convStrain(x, model, iext=0):
    """Convolution of the strain distribution

    Arguments:
        x -- scan positions [mm] along xdir
        model -- MCCfit model for intensity
        iext -- model for extinction
                0 = small radius approximation (fast)
                1 = exact (slow)
    Returns
        [I, err, pos, epos, width] = convoluted strain, error, info depth, error, info width
    """
    global shape, _sampling
    xdir = model.xdir
    nev = model.nev
    # initialize local variables
    nx = x.shape[0]
    yc = np.zeros(nx)
    yc2 = np.zeros(nx)
    yd = np.zeros(nx)
    yd2 = np.zeros(nx)
    xs = np.matrix(x)
    ex = np.zeros(nx)
    if (_sampling.nev != nev):
        _sampling.setRange(nev)
    rnd = range(0, _sampling.nev)
    # where to find r, ki, kf, p, dhkl
    [jr, jki, jkf, jp, jd] = _sampling.idata[0:5]
    d0 = _sampling.d0
    # scan positions relative to r
    dr = np.array(xs.T*xdir)

    # loop through sampling events to make convolution
    sumw = 0.
    sumn = 0
    for ir in rnd:
        rn = _sampling.sdata[ir, jr:jr+3] - _sampling.sctr
        r0 = shape.getLocalPos(rn)
        ki = shape.getLocalDir(_sampling.sdata[ir, jki:jki+3])
        kf = shape.getLocalDir(_sampling.sdata[ir, jkf:jkf+3])
        p = _sampling.sdata[ir, jp]
        # pseudo-strain at the samplig point
        deps = (_sampling.sdata[ir, jd] - d0)/d0
        r = r0 + dr
        # get event depths and flags for inside/outside events
        [d, d2, ins] = shape.depthLocal(r)
        # calculate extinction factor
        if (iext == 1):
            ex = getExtinction2(r, ki, kf)
        else:
            ex = getExtinction(r, ki, kf)
        # total weight including extinction and scattering probability
        I0 = model.intFnc(d)
        # eps = intrinsic strain + pseudo-strain
        eps = model.strainFnc(d) + deps
        w = p*ex*ins*I0
        wd = w*d
        weps = w*eps
        sumw += w
        sumn += ins
        yc += weps
        yc2 += weps*eps
        yd += wd
        yd2 += wd*d
    # This will avoid division by zero warnings:
    sg = np.array((sumn > 0) & (sumw > 0), dtype=int)
    yc += (1-sg)*1.
    sumn += (1-sg)*1
    sumw += (1-sg)*1
    eps = sg*yc/sumw
    err = sg*yc2/sumw
    pos = sg*yd/sumw
    epos = sg*yd2/sumw
    width = sg*np.sqrt(np.absolute(epos - pos**2))
    err = sg*np.sqrt(np.absolute(err - eps**2)/sumn) # + (1-sg)*np.average(eps)  
    epos = width/np.sqrt(sumn) # + (1-sg)*np.average(pos)  
    return [1e6*eps, 1e6*err, pos, epos, 2.3548*width]
