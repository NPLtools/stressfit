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

class Sampling():
    def __init__(self, events=None, columns=None, ctr = None):
        """Set array of sampling events.
    
        The sampling events can be obtained by MC ray-tracing simulation
        of the instrument at given setup. Generation of such a file
        is implemented in SIMRES ver. 6.3.5 and above
    
        Arguments
        ---------
            events -- array[:,:]
            columns -- defines column indices for r, ki, kf, p, dhkl
                        r -- event position
                        ki, kf -- initial and final wave vector
                        p -- probability
                        dhkl -- dhkl sampled by this event
        """
        self.sdata = events
        self.idata = columns
        self.src = None
        if (ctr is not None):
            self.sctr = np.array(ctr)
        else:
            self.sctr = np.array([0., 0., 0.])
        if self.sdata is not None:
            self.setRange(self.sdata.shape[0], setd0=True)

    @classmethod
    def fromdict(cls, source, ctr=None):
        """Create class from dictionary.
        
        Requires at least following keys:
            - data
            - columns
            - ctr (optional)
        
        Parameters
        ----------
        source: array[:,:]
            scattering events coordinates
        ctr: list(3)
            gauge centre coordinates
        
        """
        if 'data' in source and 'columns' in source:
            if (ctr is None) and ('ctr' in source):
                ctr = source['ctr']
            s = cls(events=source['data'], columns=source['columns'])
            s.src = source
        else:
            raise Exception('Required fields not found in source.')
        return s
    
    
    def setRange(self, nev, setd0=False):
        [jr, jki, jkf, jp, jd] = self.idata[0:5]
        self.nev = nev
        sumd = 0
        sump = 0
        for ir in range(self.nev):
            p = self.sdata[ir, jp]
            dhkl = self.sdata[ir, jd]
            sump += p
            sumd += p*dhkl
        if (setd0):
            self.d0 = sumd/sump  
        self.sump = sump

def getSampling(nev):
    _sampling.setRange(nev)
    return _sampling

def setSampling(source, ctr=None):
    """Define sampling from dict."""
    global _sampling
    _sampling = Sampling.fromdict(source, ctr=ctr)
    

def rotate(v, axis, angle):
    """Simple vector rotation function

    Arguments:
        v -- 3D vector to be rotated
        axis -- rotation axis index 0 .. 2
        angle -- rotation angle in rad
    Returns
        rotated vector
    """
    s, c = np.sin(angle), np.cos(angle)
    if (axis == 1):
        R = np.array([[c, 0., s], [0., 1., 0.], [-s, 0., c]])
    elif (axis == 0):
        R = np.array([[1., 0., 0.], [0., c, -s], [0., s, c]])
    elif (axis == 2):
        R = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
    else:
        R = np.eye(3)
    vv = np.array(v).reshape((3, 1))
    r = R.dot(vv).reshape(3,)
    return r


def setExtinction(mu = None, table=None):
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
    arg = -(p.T.dot(mu.T))
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
    x = np.array([1., 1.])
    if (ext_mu is not None):
        lam = 2.*np.pi/norm(ki)
        mu = x*np.interp(lam,
                         ext_lambda, ext_mu, left=ext_mu[0], right=ext_mu[-1])
    else:
        mu = ext_coeff*x
    arg = -p.dot(mu.T)
    res = ins*np.exp(arg)
    return res

def shuffleEvents():
    """get random sample of nev events from sdata"""
    idx = np.arange(_sampling.sdata.shape[0])
    np.random.shuffle(idx)
    sel = _sampling.sdata[idx,:]
    _sampling.sdata = sel

def convGauge(x, xdir, iext, nev):
    """Calculate sampling volume center and pseudo strains

    Arguments
    ---------
        x:  array
            scan positions [mm] along xdir
        xdir: array
            scan direction in local coordinates
        iext:  int
            model for extinction
            0 = small radius approximation (fast)
            1 = exact (slow)
        nev:  int
            number of events from sdata to integrate
    Returns
        [x, pos, width, eps] as array
        x: scan positions [mm]
        pos: information depth
        width: information width (FWHM)
        eps:  pseudo-strain (NOT converted to 10^-6 units!)       
    """
    global shape, _sampling
    # initialize local variables
    nx = x.shape[0]
    yc = np.zeros(nx)
    yd = np.zeros(nx)
    yi = np.zeros(nx)
    yd2 = np.zeros(nx)
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
        eps = deps
        w = p*ex*ins
        wd = w*d
        weps = w*eps
        sumw += w
        yi += w
        yc += weps
        yd += wd
        yd2 += wd*d
    # This will avoid division by zero warnings:
    sg = np.array( (sumw > 0) , dtype=int)
    yc += (1-sg)*1.
    sumw += (1-sg)*1
    cnts = sg*yi/_sampling.sump
    eps = sg*yc/sumw
    pos = sg*yd/sumw
    epos = sg*yd2/sumw
    width = 2.3548*sg*np.sqrt(np.absolute(epos - pos**2))
    return np.array([x, pos, width, cnts, eps]).T

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
    # where to find r, ki, kf, p, dhkl
    if (_sampling.nev != nev):
        _sampling.setRange(nev)
    rnd = range(0, nev)
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
    val = yc/sump
    err = yc2/sump    
    dval = yd/yc
    derr = yd2/yc
    # Calculate sampling position, width, intensity and error
    pos = sg*dval
    width = sg*np.sqrt(np.absolute(derr - pos**2))
    err = sg*np.sqrt(np.absolute(err - val**2)/sumn) + (1-sg)*np.average(val)  
    epos = width/np.sqrt(sumn) + (1-sg)*np.average(pos)  
    val = sg*val
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
    rnd = range(0, nev)
    nx = x.shape[0]
    yc = np.zeros(nx)
    yc2 = np.zeros(nx)
    yd = np.zeros(nx)
    yd2 = np.zeros(nx)
    xs = np.matrix(x)
    ex = np.zeros(nx)
    # where to find r, ki, kf, p, dhkl
    if (_sampling.nev != nev):
        _sampling.setRange(nev)
    rnd = range(0, nev)
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
    err = sg*np.sqrt(np.absolute(err - eps**2)/sumn) + (1-sg)*np.average(eps)  
    epos = width/np.sqrt(sumn) + (1-sg)*np.average(pos)  
    return [1e6*eps, 1e6*err, pos, epos, 2.3548*width]
