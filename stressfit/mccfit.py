# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 23:53:16 2017.

@author: Jan Saroun, saroun@ujf.cas.cz
"""
# TODO Extract reporting to a dedicated module
# TODO Model definition data and fitting functions should be in separate modules/classes.
# TODO Move the MC convolution procedures from smaple module to an extra module.
# TODO The sample module should only encapsulate sample properties.

import numpy as np
from lmfit import Minimizer, Parameters
from scipy.interpolate import CubicSpline as Cspline
from scipy.interpolate import PchipInterpolator as PCHspline
from scipy.interpolate import Akima1DInterpolator as Aspline
from scipy.interpolate import splrep, splev

import datetime
#from scipy.interpolate import splrep, splev
#import scipy.integrate as integrate
from abc import ABC, abstractmethod
# imports from stressfit
import stressfit.sample as sam
import stressfit.graphs as gr
import stressfit.dataio as dataio

_chi2 = np.inf
_reg = np.inf
_ispline = None # spline interpolator for intensity, used by SFit
_quiet = False

_log = dataio.logger()

_prog = dataio.FitProgressLogger()

intpmodels = ['natural','clamped', 'PCHIP', 'Akima']

def intClear():
    """Clear intensity fit result.
    
    Scattering probability distribution is assumed uniform afterwards.
    """
    global _ispline
    _ispline = None

def intDefined():
    """Query if intensity fit is defined."""
    return _ispline is not None

def intFnc(x):
    """Return interpolated intensity function value."""
    if (_ispline is None):
        y = np.ones(x.shape)
    else:
        y = splev(x, _ispline, ext=1)
    return np.maximum(0., y)

def quiet(b):
    """Set quiet mode for progress logs."""
    global _quiet
    _quiet=b

def path2win(name):
    """Convert windows-like paths with backslashes.
    
    Obsolete, left for backward compatibility. 
    Use pathlib or os packages instead.
    """
    out = name.replace('\\','/')
    if (out[-1] != '/'):
        out += '/'
    return out

def deriveFilename(file, ext='', sfx=''):
    """Return string derived from given filename.
    
    1)  Remove file extension (if any)
    2)  add suffix 
    3)  add extensio
     
    Parameters
    ----------
         file: string
             base filename
         ext: string
             extension to be added
         sfx: string
             a suffix to be added at the end of filename (before extension)
    
    """
    a = file.split('.')
    if (len(sfx) > 0):
        sf = '_'+sfx
    else:
        sf = ''
    if (len(ext) > 0):
        ex = '.'+ext
    else:
        ex = ''        
    if (len(a)>1):
        s = '.'.join(a[0:len(a)-1])
        s += sf + ex
    else:
        s = file + sf + ex
    return s 

def params2array(params):    
    """Convert params to an array.
    
    Returns
    -------
        array of parameters and std. errors
    """
    npar = len(params)
    res = np.zeros((npar,2))
    i = -1
    for p in params.values():
        i += 1
        res[i,0] = p.value
        res[i,1] = p.stderr
    return res



def params2dist(params):    
    """Convert distribution part of params to array.
    
    Returns
    -------
        array of [x, y, xerr, yerr] 
    """
    nm = int((len(params) - 3)/2)  
    x = np.zeros(nm)
    dx = np.zeros(nm)
    y = np.zeros(nm)
    dy = np.zeros(nm)
    for i in range(nm):
        p =  params['x_'+str(i)]
        x[i] = p.value
        if p.stderr: dx[i] = p.stderr
    for i in range(nm):
        p = params['y_'+str(i)]
        y[i] = p.value 
        if p.stderr: dy[i] = p.stderr
    return np.array([x, y, dx, dy]).T

def array2params(par, err, params):    
    """Convert array to params.
    
    NOTE: we rely that params is an ordered dictionary
    
    Parameters
    ----------
        par, err : array like
            values and std. errors
    
    Returns
    -------
        dict
            Fit parameters in the format used by mcfit.
    """
    dim = np.array([par.size, err.size, len(params)])
    b = dim - dim[0]
    if b.any():
        raise RuntimeError('Incompatible array and Parameters dimensions.')
    n = dim[0]
    nm = int((n - 3)/2)    
    
    params['A'].value = par[2*nm]
    params['A'].stderr = err[2*nm]

    params['B'].value = par[2*nm+1]
    params['B'].stderr = err[2*nm+1]

    params['xc'].value = par[2*nm+2]
    params['xc'].stderr = err[2*nm+2]
    
    
    for i in range(nm):
        params['x_'+str(i)].value = par[i]
        params['x_'+str(i)].stderr = err[i]
        params['y_'+str(i)].value = par[i+nm]
        params['y_'+str(i)].stderr = err[i+nm]



def getSmoothness(x, y):
    """Get smoothness for regularization.
 
    Returns
    -------
    float
        average of |dy/dx|
    """
    nx = x.shape[0]
    sumy = 0.
    # find index range closest to the given interval
    for i in range(nx-1):
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]        
        if (dx <= 0):
            raise RuntimeError('getSmoothness: x is not monotonously increasing')
        
        sumy += abs(dy/dx)
    return sumy/(nx-1)

def getAverage(x, y, rang=[0., 1.]):
    """Get strain integral over given range."""
    nx = x.shape[0]
    ix = [0, nx-1]
    dx = np.zeros(len(x))
    # find index range closest to the given interval
    for i in range(nx-1):
        dx[i] = x[i+1] - x[i]
        if (dx[i] <= 0):
            raise RuntimeError('getAverage: x is not monotonously increasing')
        for k in range(2):
            res1 = (x[i]-rang[k])/dx[i]
            res2 = (x[i+1]-rang[k])/dx[i]
            if (res1 <= 0. and res2 >= 0.):
                if (abs(res1) < abs(res2)):
                    j = 0
                else:
                    j = 1
                ix[k] = i+j
    dx[nx-1]=dx[nx-2]
    # integrate and divide by the interval    
    #av = np.average(y[ix[0]:ix[1]+1])
    av = np.sum(y[ix[0]:ix[1]+1]*dx[ix[0]:ix[1]+1])
    # av = integrate.simps(y[ix[0]:ix[1]+1], x=x[ix[0]:ix[1]+1])/(x[ix[1]] - x[ix[0]])
    return av
    

def getChi2(resid, params):
    """Return reduced chi2 from the previously calculated residuals."""
    ivar = 0
    for key in params:
        p = params[key]
        ivar += p.vary
    chi2 = np.sum(resid**2)/(resid.size - ivar)
    return chi2

    
def fitcallb(params, iteration, resid, *fcn_args, **fcn_kws):
    """Perform actions after each iteration.
    
    Callback function to provide fit progress info.
    """
    global _chi2, _reg
    reg = getChi2(resid, params)
    chi = getChi2(fcn_args[0].resid, params)
    if (reg < _reg):
        _reg = reg
        _chi2 = chi
        it = max(0,iteration)
        args = {'iter': it, 'chi2': chi, 'reg': reg}
        _prog.prog(**args)
        #fmt = 'iter={:d}, chi2={:g}, reg={:g}'
        #if (not _quiet): _log.progress(fmt.format(it, chi, reg))


# define objective function for fitting: returns the array to be minimized
def costFnc(params, model, areg, guess):
    """Residual array to be passed to the minimizer: intensity scan."""
    model.params = params
    res = model.calResid(guess=guess)
    nres = res.shape[0]
    nreg = 0
    # add regularization term
    if (areg > 0):
        # model.updateDistribution()
        der = model.calDeriv()
        reg =np.sqrt(areg)*np.absolute(der)
        nreg = reg.shape[0]       
        w = np.sqrt(1.0*abs(nres + nreg - model.nvar)/abs(nres - model.nvar))
        res = np.concatenate((res, reg), axis=0)*w
    if (model.constraint is not None):
        w = np.sqrt(1.0*abs(nres + nreg - model.nvar))
        y = w*model.constraint(params)
        res = np.concatenate((res, np.array([y])), axis=0)
    return res

def runFit(model, maxiter=200, areg=0., bootstrap=False, loops=3, guess=False):
    """Execute least square fit.
    
    Parameters
    ----------
        model : mccfit class
            class handling model parameters
        maxiter : int
            maximum number of function evaluations
        areg : float
            regularization coefficient
        bootstrap : boolean
            run loops to estimate confidence limits
        loops : int
            number of bootstrap loops
        guess : bool
            if true, run only the guess fit without deconvolution
    Returns
    -------
        results : boolean
                success
    """
    global _chi2, _reg
    if (model.data is None):
        msg = 'Set value of {}.data.'.format(model.__class__.__name__)
        raise Exception('No data defined for fitting. '+msg)

    if (bootstrap):
        res = runFitEx(model, maxiter=maxiter, loops=loops, areg=areg)
        return res
    if not (_quiet or guess) and maxiter>0:
        _prog.start(**{'iter': maxiter})
        #_log.clear(what='prog')
        #_log.progress('Starting fit for < {:d} iterations.'.format(maxiter))
    ndim = model.dist.shape[0]
    ndata = model.data.shape[0]
    xdata = model.data[:,0]
    _chi2 = np.inf
    _reg = np.inf
    # initialize distribution array
    xdis = model.dist[:,0]
    model.fitInit()
    if (maxiter>0):
        if guess:
            iterf = None
        else:
            iterf = fitcallb
        minner = Minimizer(costFnc, model.params, fcn_args=(model, areg, guess), iter_cb=iterf)
        model.result = minner.minimize(method='leastsq', 
                        params=model.params, ftol=model.ftol, epsfcn=model.epsfcn, max_nfev=maxiter)
        model.params = model.result.params
    model.areg = areg
    # par =params as array (par, serr)
    model.par = params2array(model.params)
    # fit = array(x, y, yerr) with fitted curve
    [yfit, fit_err, ypos] = model.getSmearedFnc(model.data[:,0], guess=False)
    model.fit = np.array([xdata, yfit, fit_err]).T
    # pos = array(x, y, yerr), with nominal and information depths
    model.pos = np.array([xdata, ypos, np.zeros(ndata)]).T
    # guess = array(x, y, yerr) with the guess function (no convolution)
    if (guess):
        [yfit, fit_err, ypos] = model.getSmearedFnc(model.data[:,0], guess=True)
        model.guess = np.array([xdata, yfit, fit_err]).T
        model.params_ini = model.params
    # dist = array(x, y, xerr, yerr) with interpolated distribution
    ydis = model.getDistribution(xdis)
    dis_err = np.zeros(ndim)
    model.dist = np.array([xdis, ydis, dis_err]).T
    model.chi  = _chi2 
    model.reg = _reg
    if (model.avgrange is not None):
        av = getAverage(xdis, ydis, rang = model.avgrange)  
        model.avg = [av, 0.]
    #if not (_quiet or guess) and maxiter>0: 
        #_log.progress('runFit finished')
    if (model.result is not None):
        res = model.result.success
    else:
        res = True
    model.fitFinal()
    if not guess and maxiter>0: 
        _prog.finished(**{'completed':res, 'chi2':model.chi, 'reg': model.reg})
    #if not res:
    #    _log.progress('Fit not completed.')
    return res

def runFitEx(model, maxiter=200, loops=5, areg=0, guess=False):
    """Execute a least square fit.
    
    Extended version with error estimates by bootstrap method.

    Parameters
    ----------
        model : mccfit class
            class handling model parameters
        maxiter : int
            maximum number of function evaluations
        loops : int
            number of loops in the boostrap cycle
        areg : float
            regularization coefficient
        guess : bool
            if true, run only the guess fit without deconvolution

    Returns
    -------
        results : boolean
                    success
    """
    global _chi2, _reg
    if (model.data is None):
        msg = 'Set value of {}.data.'.format(model.__class__.__name__)
        raise Exception('No data defined for fitting. '+msg)
    
    if not (_quiet or guess) and maxiter>0:
        #_log.clear(what='prog')
        #fmt = 'Starting fit for < {:d} iterations'
        #fmt += ' and {:d} loops for error estimate.'
        #_log.progress(fmt.format(maxiter, loops))
        _prog.start_loops(**{'iter': maxiter, 'loops': loops})
    ndim = model.dist.shape[0]
    # length of distribution model arrays (number of nodes):  
    # initialize arrays
    pval = np.zeros(loops)
    ndata = model.data.shape[0]
    xdata = model.data[:,0]
    yfit = np.zeros(ndata)
    yfit2 = np.zeros(ndata)
    ypos = np.zeros(ndata)
    ypos2 = np.zeros(ndata)
    ydis = np.zeros(ndim)
    ydis2 = np.zeros(ndim)    
    ypar = np.zeros(len(model.params))
    ypar2 = np.zeros(len(model.params))
    data0 = model.data.copy()
    chi0 = 0
    reg0 = 0
    sump = 0.
    avg = 0.
    avg2 = 0.
    # initialize distribution array
    xdis = model.dist[:,0]
    model.fitInit()
    for it in range(loops):
        _chi2 = np.inf
        _reg = np.inf
        model.resetParams()
        sam.shuffleEvents()
        model.data = model.addNoise(data0)
        minner = Minimizer(costFnc, model.params, fcn_args=(model, areg, guess), iter_cb=fitcallb)
        model.result = minner.minimize(method='leastsq', 
                            params=model.params, ftol=model.ftol, epsfcn=model.epsfcn, max_nfev=maxiter)
        model.params = model.result.params
        # get derived results
        dist = model.getDistribution(xdis)
        if (model.avgrange is not None):
            av = getAverage(xdis, dist, rang = model.avgrange)  
        [y, ey, pos] = model.getSmearedFnc(xdata, guess=guess)
        par = params2array(model.params)
        _reg = model.result.redchi
        _chi2 = getChi2(model.resid, model.params)  
        p = np.exp(-_reg)
        pval[it] = p
        sump += p
        chi0 += p*_chi2
        reg0 += p*_reg
        # accumlate statistics
        yfit += p*y
        yfit2 += p*y**2
        ypos += p*pos
        ypos2 += p*pos**2
        ydis += p*dist
        ydis2 += p*dist**2
        ypar += p*par[:,0]
        ypar2 += p*par[:,0]**2
        if (model.avgrange is not None):
            avg += p*av
            avg2 += p*av**2
        #reg = costFnc(model.params, model, areg)
        #if not  model.result.success:
        #    _log.progress('Fit not completed.')
        if not (_quiet or guess):
            #_log.progress('Loop {:d}: chi2={:g}, reg={:g}\n'.format(it, _chi2, _reg))
            args = {'loop':it, 'chi2':_chi2, 'reg':_reg}
            _prog.finished_loop(completed=model.result.success, **args)
    model.data = data0
    model.chi  = chi0/sump
    model.reg = reg0/sump
    model.areg = areg
    ydis = ydis/sump
    ydis2 = ydis2/sump
    dis_err = np.sqrt(np.abs(ydis2 - ydis**2))
    yfit = yfit/sump
    yfit2 = yfit2/sump
    fit_err = np.sqrt(np.abs(yfit2 - yfit**2))
    ypos = ypos/sump
    ypos2 = ypos2/sump
    pos_err = np.sqrt(np.abs(ypos2 - ypos**2) )   
    ypar = ypar/sump
    ypar2 = ypar2/sump
    ypar_err = np.sqrt(np.abs(ypar2 - ypar**2))
    avg = avg/sump
    avg2 = avg2/sump
    avg_err = np.sqrt(np.abs(avg2 - avg**2))    
    # save results
    # copy parameters
    model.params = model.result.params
    # par =params as array (par, serr)
    model.par = np.array([ypar,ypar_err]).T
    # set params to means and stderr from the bootstrap
    array2params(ypar,ypar_err, model.params)
    # fit = array(x, y, yerr) with fitted curve
    model.fit = np.array([xdata, yfit, fit_err]).T
    # pos = array(x, y, yerr), with nominal and information depths
    model.pos = np.array([xdata, ypos, pos_err]).T
        # guess = array(x, y, yerr) with the guess function (no convolution)
    if (guess):
        [yfit, fit_err, ypos] = model.getSmearedFnc(model.data[:,0], guess=True)
        model.guess = np.array([xdata, yfit, fit_err]).T
        model.params_ini = model.params
    # dist = array(x, y, xerr, yerr) with interpolated distribution
    model.dist = np.array([xdis, ydis, dis_err]).T
    # averaged value
    model.avg = [avg, avg_err]
    if (model.result is not None):
        if not (_quiet or guess):
            #_log.progress('runFitEx weights: '+str(pval)+'\n')
            _prog.log('Loop weights: '+str(pval)+'\n')
        res = model.result.success
    else:
        res = True
    model.fitFinal()
    return res


class MCCfit(ABC):
    """Set parameters needed to run the fit procedure.
        
    Arguments
    ---------                
    nev : int
        number of points from sampling events used in convolution
    xdir : int
        scan direction (local coordinates)  
    ftol : float
        tolerance to be passed to the leastsq minimizer
    epsfcn : function
        Function used by scipy.optimize.leastsq to calculate optimal step length.
        
    """
    
    def __init__(self, nev=1000, xdir=[0., 0., -1.], ftol=1.e-4, epsfcn=None):

        # assign main logger for messages
        self._log = dataio.logger()
        # init parameters
        self.nev = nev
        xd = np.array(xdir)
        nm = np.sqrt(np.sum(xd*xd))
        self.xdir = xd/nm
        self.ftol = ftol
        self.epsfcn = epsfcn
        self.nd = 0 # set by defDistribution
        # current and inital parameters
        self.params = None
        self.params_ini = None
        self.data = None # array(x, y, yerr) with actually fitted data
        # results
        self.resid = None # residuals
        self.par = None  # params as array (par, err)
        self.mod = None  # disttribution model as array (x, y, xerr, yerr)
        self.result = None  # Results class 
        self.fit = None  # array(x, y, yerr) with fitted curve
        self.dist = None # array(x, y, xerr, yerr) with interpolated distribution
        self.tck = None # interpolation polynomial for dist
        self.pos = None # array(x, y, yerr), with nominal and information depths
        self.infodepth = None # array(x, y, width, int, eps) 
        self.resolution = None # dict with spatial resolution parameters
        self.avg = [0., 0.] # [avg, err_avg], distribution average over given interval calculated by RunFit or RunFitEx
        self.avgrange = None # interval for averaging
        self.chi = 0.
        self.reg = 0.
        self.areg = 0.
        self.constraint = None
        self.guess = None
        self.nvar = 0 # number of free parameters
        
        # interpolation data and model
        self.intp = None
        self.intpmod = 1
        self.scaled = False # if true, scaling factor and shift are included in reports
        
    def fitInit(self):
        """Do common tasks before fitting."""
        ivar = 0
        for key in self.params:
            p =self.params[key]
            ivar += p.vary
        self.nvar = ivar
        self.resetParams()
        xdata = self.data[:,0]
        self.calInfoDepth(xdata)
        
    def fitFinal(self):
        """Do common tasks after fitting."""
        # mod = disttribution model as array (x, y, xerr, yerr)
        self.mod = params2dist(self.params)
        self.setModelPoints(self.mod[:,0], self.mod[:,1])
               
    def defDistribution(self, par, vary, ndim=100, scaled=False):
        """Define distribution model parameters.
    
        Creates Parameters object required by Minimizer from lmfit.
        
        Parameters
        ----------
        par : list [x, y]
            node positions (x) and values (y)
        vary : list [fx, fy]
            corresponding flags 0|1 for fixed|free parameters
            interpolation order (1..3)
        ndim: int
            number of points for interpolation
        scaled : bool
            If true, scaling parameters will be included in the result of 
            getDistribution.
        """
        # validate input and plot
        self.scaled = scaled
        x = np.array(par[0])
        y = np.array(par[1])
        fx = np.array(vary[0])
        fy = np.array(vary[1])
        self.nd = x.size  # number of dept/value pairs
        dim = np.array([x.size, y.size, fx.size, fy.size])
        b = dim - dim[0]
        if b.any():
            raise RuntimeError('Check the input arrays. They must have equal length.')
            
        if (self.nd < 2):
            raise RuntimeError('At least 2 points are needed to define a distribution.')

        self.params = Parameters()   
        for i in range(self.nd):
            if (fx[i]):
                if (i < self.nd-1):
                    maxx = 0.5*(x[i] + x[i+1])
                else:
                    maxx = x[i]
                if (i > 0):
                    minx = 0.5*(x[i] + x[i-1])
                else:
                    minx = x[i]
                self.params.add('x_'+str(i), value=x[i], vary=fx[i], min=minx, max=maxx)
            else:
                self.params.add('x_'+str(i), value=x[i], vary=fx[i])
        for i in range(self.nd):
            self.params.add('y_'+str(i), value=y[i], vary=fy[i])
        # define basic scaling as fixed 
        self.params.add('A', value=1., vary=0)
        self.params.add('B', value=0., vary=0)
        self.params.add('xc', value=0., vary=0)
        # store as initial parameters
        self.params_ini = self.params
        # initialize interpolated values
        self.mod = params2dist(self.params)
        xdis = np.linspace(self.mod[0,0], self.mod[-1,0], num=ndim, endpoint=True)
        ydis = self.getDistribution(xdis)
        dis_err = np.zeros(ndim)
        self.dist = np.array([xdis, ydis, dis_err]).T

    def addNoise(self, data):
        """Return fit data + gaussian noise.
        
        Assume data = array(x,y,err)
        """
        n = data.shape[0]
        x = np.random.normal(size=n)
        e = data[:,2]*x
        res = data.copy()
        res[:,1] += e
        res[:,2] = np.sqrt(2.0)*data[:,2]
        return res


    def resetParams(self):
        """Reset params to initial values."""
        self.params = self.params_ini
        
    def defScaling(self, par, vary, minval=None, maxval=None):
        """Define scaling for the distribution model.
        
        Adds items to the Parameters object required by Minimizer from lmfit.
        Call to this routine assumes previous call to defDistribution !
        
        Parameters
        ----------
        par : [A, B, xc]
            scaling for the distribution: y = A*dist(x-xc) + B
        vary : [fA, fB, fxc]
            corresponding flags 0|1 for fixed|free parameters
        """
        # validate input and plot
        dim = np.array([len(par), len(vary), 3])
        b = dim - dim[0]
        if b.any():
            raise RuntimeError('Check the input. Expected par = [A, B, zc].')            
        if (self.params is None):
            raise RuntimeError('Call defDistribution first to initialize Parameters.')
        if (minval is None):
            vmin = -np.inf*np.array([1., 1., 1,])
        else:
            vmin = np.array(minval)
            
        if (maxval is None):
            vmax = np.inf*np.array([1., 1., 1,])
        else:
            vmax = np.array(maxval)  
        self.params.add('A', value=par[0], vary=vary[0], min=vmin[0], max=vmax[0])
        self.params.add('B', value=par[1], vary=vary[1], min=vmin[1], max=vmax[1])
        self.params.add('xc', value=par[2], vary=vary[2], min=vmin[2], max=vmax[2])
        self.params_ini = self.params

    def getScaling(self):
        """Extract scaling parameters.
            
        Returns
        -------
            Model parameters and scaling : [ A, B, xc]

        """
        if ('A' in self.params):
            A = self.params['A'].value
        else:
            A = 1.
        if ('B' in self.params):
            B = self.params['B'].value
        else:
            B = 0.
        if ('xc' in self.params):
            xc = self.params['xc'].value
        else:
            xc = 0.
        return [A, B, xc]
    def getParams(self):
        """Convert params to a list of parameter values.
            
        Returns
        -------
            Model parameters and scaling : [xy, A, B, xc]
        """
        [A, B, xc] = self.getScaling()
        xy = params2dist(self.params)
        #x = np.zeros(self.nd)
        #y = np.zeros(self.nd)
        #for i in range(self.nd):
        #    x[i] = self.params['x_'+str(i)].value
        #for i in range(self.nd):
        #    y[i] = self.params['y_'+str(i)].value
        # return [x, y, A, B, xc]
        return [xy, A, B, xc]


    def getDistribution(self, xval):
        """Return interpolated distribution."""
        [xy, A, B, xc] = self.getParams()
        self.setModelPoints(xy[:,0], xy[:,1])
        if (self.scaled):
            yval = self.intp(xval)
            res = A*yval + B
        else:
            res = self.intp(xval)
        return res

    def updateDistribution(self):
        """Shortcut to calculate internally stored distribution."""
        dis = self.getDistribution(self.dist[:,0])
        self.dist[:,1] = dis[:]
        self.tck = splrep(self.dist[:,0], self.dist[:,1], k=1, s=0)
        return

    def calResid(self, guess=False):
        """Calculate and return residuals for leastsq fiting."""
        [yf, ef, pos] = self.getSmearedFnc(self.data[:,0], guess=guess)
        er = np.sqrt(ef**2 + self.data[:,2]**2)
        self.resid = (yf - self.data[:,1])/er
        return self.resid

    def calDeriv(self):
        """Calculate 1st derivative of the distribution."""
        nx = self.dist.shape[0]
        res = np.zeros(nx-1)
        x = self.dist[:,0]
        y = self.dist[:,1]
        for i in range(nx-1):
            dx = x[i+1] - x[i]
            dy = y[i+1] - y[i] 
            res[i] = dy/dx
        return res

    def formatResultModel(self):
        """Format a table with deconvoluted curves.
        
        Returns
        -------
            String
                tab-delimited table with self.dist and self.mod arrays 
        """
        ss = "# Date: " + str(datetime.datetime.now()) + "\n"
        hdr = ['xi', 'yi', 'yi_err', 'x', 'y', 'x_err', 'y_err']
        ss += '# Header: ' + '\t'.join(f for f in hdr) + '\n'
        # number of rows in the arrays
        hasFit = (self.mod is not None)
        hasFit = hasFit and (self.dist is not None)
        if hasFit:
            nm = self.mod.shape[0]
            nf = self.dist.shape[0]
            nc = self.dist.shape[1]
            for i in range(max(nm, nf)):
                if (i < nm):
                    a = self.mod[i,]
                    sm = '\t'.join('{:g}'.format(f) for f in tuple(a))                
                else:
                    sm = ''
                if (i < nf):
                    a = self.dist[i,]
                    sf = '\t'.join('{:g}'.format(f) for f in tuple(a))
                else:
                    sf = '\t' + '\t'*(nc-1)
                ss += sf + '\t' + sm + '\n'
        return ss
   
    
    def formatResultFit(self):
        """Format a table with fitted data.
        
        Appends fitted curve to the table of original data.
        The resulting table has 6 columns: x, y, err, pos, fit, err.
        Adds header with time and header info.
        
        Returns
        -------
            String
                tab-delimited table with self.data self.pos and self.fit[:,1:3] 
        """
        ss = "# Date: " + str(datetime.datetime.now()) + "\n"
        hasData = (self.data is not None)
        hasFit = ((self.fit is not None) and (self.pos is not None))
        hasGuess = (self.guess is not None)
        if (hasData):
            nm = self.data.shape[0]
            nr = 3
            hdr = ['x', 'y', 'err']
            if hasFit:
                nr += 3
                hdr.extend(['infdepth', 'fit', 'err'])
            if (hasGuess):
                nr += 1
                hdr.extend(['guess'])
            fmt = '{:g}'+'\t{:g}'*(nr-1)
            ss += '# Header: ' + '\t'.join(f for f in hdr) + '\n'   
            for i in range(nm):
                row = []
                x = self.data[i,0]
                y = self.data[i,1]
                err = self.data[i,2] 
                row = [x, y, err]
                if hasFit:
                    fmt = '{:g}'+'\t{:g}'*5
                    pos = self.pos[i,1]
                    yf = self.fit[i,1]
                    ef = self.fit[i,2]
                    row.extend([pos, yf, ef])
                else:
                    sm = '{:g}\t{:g}\t{:g}'.format(x, y, err)
                if (hasGuess):
                    yg = self.guess[i,1]
                    row.extend([yg])
                sm = fmt.format(*row)
                ss += sm + '\n'
        else:
            ss += 'No data, no fun ...\n'
        return ss

    def formatResultDepth(self):
        """Format a table with depth scale.
        
        Returns
        -------
            String
                tab-delimited table with self.infodepth 
        """
        ss = "# Date: " + str(datetime.datetime.now()) + "\n"
        hdr = ['position', 'depth', 'x', 'y', 'z', 'width', 'intensity', 'pseudo_strain']
        ss += '# Header: ' + '\t'.join(f for f in hdr) + '\n'   
        hasData = (self.infodepth is not None)
        if (hasData):
            nm = self.infodepth.shape[0]
            fmt = '{:g}'+7*'\t{:g}'
            for i in range(nm):
                row = self.infodepth[i,[0,1,5,6,7,2,3,4]]
                sm = fmt.format(*row) 
                ss += sm + '\n'
        return ss


    def formatResultResol(self):
        """Format a table with spatial resolution parameters.
                
        Returns
        -------
            String
                tab-delimited table with self.resolution values 
        """
        ss = "# Date: " + str(datetime.datetime.now()) + "\n"
        hdr = ['position', 'depth', 'width']
        hdr.extend(['x0','y0','z0'])
        hdr.extend(['cxx','cyy','czz','cyz','cxz','cxy'])
        
        ss += '# Header: ' + '\t'.join(f for f in hdr) + '\n'   
        
        # shortcut 
        r = self.resolution
        # convert to array
        arr = np.concatenate((r['x'],r['pos'],r['width'],r['ctr'],r['cov']),
                            axis=1)
        hasData = (self.resolution is not None)
        if (hasData):
            nm = self.resolution['x'].shape[0]
            fmt = '{:g}'+11*'\t{:g}'
            for i in range(nm):
                row = arr[i,:]
                sm = fmt.format(*row) 
                ss += sm + '\n'
        return ss

    def formatResultLog(self):
        """Format a table with fitted parameters and chi2.
               
        Returns
        -------
            String
                headers and table with parameters and errors
        """
        try:
            # number of parameters
            na = len(self.params)
            val = np.zeros(na)
            err = np.zeros(na)
            vary = np.zeros(na)
            # cov = self.result.covar
            # hasCov = (cov is not None)
            ss = "# Date: " + str(datetime.datetime.now()) + "\n"
            if (self.result is not None):
                ss += "# Method: " + self.result.method + "\n"
            ss += "# Interpolation: " + intpmodels[self.intpmod]  + "\n"
            ss += "# Chi2: {:g}".format(self.chi)  + "\n"
            ss += "# Reg: {:g}".format(self.reg)  + "\n"
            ss += "# Areg: {:g}".format(self.areg)  + "\n"
            if (self.avgrange is not None):
                sm = '# Integral over x in [{:g}, {:g}]: {:g} +- {:g}\n'
                ss += sm.format(self.avgrange[0],self.avgrange[1], 
                                self.avg[0], self.avg[1])                
            hdr = ['name', 'value', 'stdev', 'vary']
            ss += '\t'.join(f for f in hdr) + '\n'
            i = -1
            for p in self.params.values():
                i += 1
                val = p.value
                vary = p.vary
                err = p.stderr
                if (err is None):
                    err = 0.
                sm = '{:s}\t{:g}\t{:g}\t{:d}'.format(p.name, val, err, vary)
                ss += sm + "\n"
            res = ss
        except Exception as e:
            res = ss
            msg = '{}\n'.format(str(e))
            msg += 'Can''t evaluate fit result: '
            msg += 'the fit was not completed successfully or wrong result structure.'
            self._log.exception(msg)
            # raise RuntimeError(strmsg)
        return res

    def saveInfoDepth(self, outpath, fname):                        
        """Save info depth table - obsolete. Use savePseudoStrain."""
        self.savePseudoStrain(outpath, fname, sfx='depth')
   
    def savePseudoStrain(self, outpath, fname, sfx='deps'):                        
        """Save pseudo-strain, intensity and information depth table."""
        if ((self.infodepth is not None) and fname):
            f = dataio.derive_filename(fname, ext='dat', sfx=sfx)
            fn = str(dataio.get_output_file(f, path=outpath))
            ss = self.formatResultDepth()
            print('Pseudo-strain data saved in '+fn)
            with open(fn, 'w') as f:
                f.write('# Pseudo-strain data: '+fname + "\n")
                f.write(ss)
                f.close             
   
    def saveResolution(self, outpath, fname, sfx='resol'):                        
        """Save table with spatial resolution parameters."""
        if ((self.resolution is not None) and fname):
            f = dataio.derive_filename(fname, ext='dat', sfx=sfx)
            fn = str(dataio.get_output_file(f, path=outpath))
            ss = self.formatResultResol()
            print('Resolution data saved in '+fn)
            with open(fn, 'w') as f:
                f.write('# Resolution data: '+fname + "\n")
                f.write(ss)
                f.close
                
    def saveResults(self, outpath, fname, reglog=None):
        """Save fit results.
        
        Parameters
        ----------
        outpath : str
            Output path.
        fname : str
            File name.
        reglog : dict
            Log data from the regularization loop if any.
        
        """
        hasData = (self.data is not None)
        hasFit = (self.fit is not None)
        if (hasData and hasFit and fname):  
            f = dataio.derive_filename(fname, ext='dat', sfx='model')
            fn = str(dataio.get_output_file(f, path=outpath))
            # fn = deriveFilename(outpath+fname, ext='dat', sfx='model')
            ss = self.formatResultModel()
            print('Model saved in '+fn)
            with open(fn, 'w') as f:
                f.write('# Fitted data: '+fname + "\n")
                f.write(ss)
                f.close
            # Save original data and fit 
            f = dataio.derive_filename(fname, ext='dat', sfx='fit')
            fn = str(dataio.get_output_file(f, path=outpath))
            #fn = deriveFilename(outpath+fname, ext='dat', sfx='fit')
            ss = self.formatResultFit()
            print('Fit saved in '+fn)
            with open(fn, 'w') as f:
                f.write('# Fitted data: '+fname + "\n")
                f.write(ss)
                f.close
                
            # Save log with parameters 
            f = dataio.derive_filename(fname, ext='log', sfx='')
            fn = str(dataio.get_output_file(f, path=outpath))
            # fn = deriveFilename(outpath+fname, ext='log', sfx='')
            ss = self.formatResultLog()
            
            if (reglog is not None):
                ss += '# Regularization log:\n'
                for ia in range(reglog.shape[0]):
                    ss += '{:g}\t{:g}\t{:g}\n'.format(*reglog[ia,:])
            
            print('Log saved in '+fn)
            with open(fn, 'w') as f:
                f.write('# Fit result: '+fname + "\n")
                f.write(ss)
                f.close
    
    def calInfoDepth(self, x, use_int=False):
        """Calculate information depth data.
        
        Performs MC convolution of the sampling distribution with the sample 
        for each scan step. Output is saved as a 2D array in self.infodepth.         
        The output columns are:
            
        -  x : scan positions [mm]
        -  pos : information depth
        -  width : information width (FWHM) 
        -  cnts : intensity (=sampling volume)
        -  eps :  pseudo-strains [1e-6] 
        -  ctr : xyz coordinates of the sampling centre of gravity
        -  err : errors of pseudo-strain [1e-6] 
        
        Parameters
        ----------
        x : array
            scan positions
        use_int : bool
            True of the convolution takes into account the scattering 
            probability distribution (it should be previously determined by 
            fitting the intensities).
        """
        ifunc = None
        if use_int:
            ifunc = intFnc
        [A, B, xc] = self.getScaling()
        self.infodepth = sam.convGauge(x-xc, self.xdir, 0, self.nev, 
                                       ifunc=ifunc)

    def calResolution(self, x, use_int=False):
        """Calculate spatial resolution characteristics.
        
        Performs MC convolution of the sampling distribution with the sample 
        for each scan step. Output is saved as a dictionary with following 
        items:
        
            
        - x : scan positions [mm]
        - pos : information depth (depends on sample shape definition)
        - width : information width (FWHM of depth)
        - ctr :  sampling centre of mass (x,y,z) in local coordinates
        - cov : covariance matrix of sampling distribution 
          (6 elements in Voigt notation)
        
        Parameters
        ----------
        x : array
            scan positions
        use_int : bool
            True of the convolution takes into account the scattering 
            probability distribution (it should be previously determined by 
            fitting the intensities).
        """
        ifunc = None
        if use_int:
            ifunc = intFnc
        [A, B, xc] = self.getScaling()
        self.resolution = sam.convResol(x-xc, self.xdir, 0, self.nev,
                                        ifunc=ifunc)
        
    def setModelPoints(self, x, y):
        """Set x,y points to define the model."""
        if not ((x is None) and (y is None)):
            nx = x.size
            ny = y.size
            if ((nx != ny) or (nx < 2)):
                raise RuntimeError("Incorrect size of x,y arrays [{}, {}".format(nx, ny))
            if (self.intpmod==0):
                self.intp = Cspline(x,y,bc_type='natural', extrapolate=True)
            elif (self.intpmod==1):
                self.intp = Cspline(x,y,bc_type='clamped', extrapolate=True)
            elif (self.intpmod==2):
                self.intp = PCHspline(x,y, extrapolate=True) 
            elif (self.intpmod==3):
                self.intp = Aspline(x,y)

    def setInterpModel(self, model='clamped'):
        """Set interpolation method as an index or string.
        
        Available interpolation models are defined by the list
        mcfit.intpmodels. The choice affects which of the 
        scipy.interpolate package methods is used for inteprolation between 
        mdel distribution nodes.
        
        """
        try:
            if (type(model) is int):
                if (model<0 or model > len(intpmodels)-1):
                    raise Exception
                self.intpmod = model
            else:
                self.intpmod = intpmodels.index(model)
        except:
            self._log.error('Invalid interpolation model.')
            msg = 'Choose string or index to define interpolation:\n'
            fmt = '{}:\t{}\n'
            for i in range(len(intpmodels)):    
                msg += fmt.format(i, intpmodels[i])
            self._log.info(msg)
            
    # Descendants must implement a smearing model
    @abstractmethod
    def getSmearedFnc(self, x, guess=False):
        """Make convolution of given distribution with sampling events.
        
        Arguments
        ---------
            x: ndarray[:]
                scan positions [mm]
            guess: boolean
                if true, get unsmeared function, only rescaled for information 
                depth and pseudostrains
        
        Returns
        -------
        y: ndarray
            values at x
        ey: ndarray
            errors at x
        pos: ndarray
            information depth at x
        """
        pass

# class for fitting smeared intensities
class Ifit(MCCfit):
    """Descendant of MCCfit class for fitting intensity.
    
    Arguments
    ---------                
    nev : int
        number of points from sampling events used in convolution
    xdir : int
        scan direction (local coordinates)  
    ftol : float
        tolerance to be passed to the leastsq minimizer
    epsfcn : function
        Function used by scipy.optimize.leastsq to calculate optimal step length.
    """
    
    def getSmearedFnc(self, x, guess=False):
        """Calculate smeared intensity distribution.
        
        Performs MC convolution of the sampling cloud of events with 
        the sample geometry and its intrinsic properties including
        distribution of intrinsic scattering probability.
        

        Parameters
        ----------
        x : array
            Scan positions [mm].
        guess : bool
            If true, no convolution is done, only returns intrinsic scattering 
            probability distribution as a function of scan positions.

        Returns
        -------
        list
            [y, ey, pos] arrays with intensity values, stddev and information 
            depth values.

        """
        [A, B, xc] = self.getScaling()
        self.updateDistribution()
        if (guess):
            res = self.intFnc(self.infodepth[:,1])*self.infodepth[:,3]
            err = np.zeros(res.shape)
            pos = self.infodepth[:,0]
        else:
            [res, err, pos, epos, width] = sam.convIntensity(x-xc, self)
        y = A*res + B
        ey = A*err
        return [y, ey, pos]

    def getDistribution(self, xval):
        """Return interpolated distribution.
        
        Overrides the default getDistribution, impose >=0 values.
        """
        d = super().getDistribution(xval)
        res = np.maximum(0., d)        
        return res        

    def defDistribution(self, par, vary, ndim=100, scaled=False):
        """Define distribution model parameters.
    
        Overrides the default getDistribution, impose >=0 values.
        
        Creates Parameters object required by Minimizer from lmfit.
        
        Parameters
        ----------
        par : list [x, y]
            node positions (x) and values (y)
        vary : list [fx, fy]
            corresponding flags 0|1 for fixed|free parameters
            interpolation order (1..3)
        ndim: int
            number of points for interpolation
        scaled : bool
            If true, scaling parameters will be included in the result of 
            getDistribution.
        """
        super().defDistribution(par, vary, ndim=ndim, scaled=scaled)
        for i in range(self.nd):
            self.params['y_'+str(i)].set(min=0.)

    def reportFit(self, outpath='', file='', plotSampling=False, **kwargs):
        """Report on the fit results.
        
        Plots and saves the results of the fitting procedure.
        
        Parameters
        ----------
        outpath : str
            Output path.
        file : str
            Output file name.
        plotSampling : bool
            If true, also plot and save the spatial resolution characteristics.
        **kwargs : 
            Allows to pass 'use_int' parameter to calInfoDepth.

        """
        if (file):
            saveit=True
            f = dataio.derive_filename(file, ext='png', sfx='fit')
            outpng = dataio.get_output_file(f, path=outpath)
            f = dataio.derive_filename(file, ext='png', sfx='depth')
            outpngdepth = dataio.get_output_file(f, path=outpath)
        else:
            saveit=False
            outpng=''
            outpngdepth=''
        # Plot result, fit & model
        gr.plot_fit(self, what='int', toplot=['fit','model'], 
                    inline=True, save=saveit, file=outpng)
        # Save results
        self.saveResults(outpath, file)

        # Plot information depth and sampling width        
        if plotSampling:
            use_int = False
            if 'use_int' in kwargs:
                use_int = kwargs['use_int']
            self.calInfoDepth(self.data[:,0], use_int=use_int)
            gr.plot_resolution(self, depth=True, cog=True, inline=True, 
                               save=saveit, file=outpngdepth)
            self.saveInfoDepth(outpath, file)

    def intFnc(self, x):
        """Return interpolated intensity function value."""
        y = splev(x, self.tck, ext=1)
        # disable negative intensities
        return np.maximum(0., y)
    
    def fitFinal(self):
        """Do common tasks after fitting ends."""
        global _ispline
        super().fitFinal()
        _ispline = self.tck


# class for fitting smeared strains
class Sfit(MCCfit):
    """Descendant of MCCfit class for fitting strain distribution.
    
    Arguments
    ---------                
    nev : int
        number of points from sampling events used in convolution
    xdir : int
        scan direction (local coordinates)  
    ftol : float
        tolerance to be passed to the leastsq minimizer
    epsfcn : function
        Function used by scipy.optimize.leastsq to calculate optimal step length.
    """
    
    def __init__(self, nev=1000, xdir=[0., 0., -1.], ftol=1.e-4, epsfcn=None):
        super().__init__(nev=nev, xdir=xdir, ftol=ftol, epsfcn=epsfcn)
        self.imodel = None
        
    def getSmearedFnc(self, x, guess=False):
        """Calculate smeared strain distribution.
        
        Performs MC convolution of the sampling cloud of events with 
        the sample geometry and its intrinsic properties including
        lattice strain distribution.

        Parameters
        ----------
        x : array
            Scan positions [mm].
        guess : bool
            If true, no convolution is done, only returns intrinsic scattering 
            probability distribution as a function of scan positions.

        Returns
        -------
        list
            [y, ey, pos] arrays with fitted strain values, stddev and 
            information depth values.

        """        
        [A, B, xc] = self.getScaling()
        self.updateDistribution()
        if (guess):
            res = 1e6*self.strainFnc(self.infodepth[:,1]) + self.infodepth[:,4]
            err = np.zeros(res.shape)
            pos = self.infodepth[:,0]
        else:
            [res, err, pos, epos, width] = sam.convStrain(x-xc, self)
        y = A*res + B
        ey = A*err
        return [y, ey, pos]

    def reportFit(self, outpath='', file='', reglog=None, **kwargs):
        """Report on the fit results.
        
        Plots and saves the results of the fitting procedure.
        
        Parameters
        ----------
        outpath : str
            Output path.
        file : str
            Output file name.
        reglog : dict
            Log data from the regularization loop if any.
        **kwargs : 
            Not used in this version. 

        """        
        if (file):
            saveit=True
            f = dataio.derive_filename(file, ext='png', sfx='fit')
            outpng = dataio.get_output_file(f, path=outpath)
        else:
            saveit=False
            outpng=''
        
        # Plot result, fit & model
        
        gr.plot_fit(self, what='eps', toplot=['fit','model'], 
                    inline=True, save=saveit, file=outpng)
        # Save results
        self.saveResults(outpath, file, reglog=reglog)
        
        # Report averaged value
        if (self.avgrange is not None):
            fmt = 'Integral over x in [{:g}, {:g}]: {:g} +- {:g}\n' 
            msg = fmt.format(self.avgrange[0], self.avgrange[1], self.avg[0], 
                             self.avg[1])
            if (not _quiet): self._log.info(msg)
        
    def intFnc(self,x):
        """Return interpolated strain function value."""
        if (_ispline is None):
            y = np.ones(x.shape)
        else:
            y = splev(x, _ispline, ext=1)
        return np.maximum(0., y)
    
    def strainFnc(self,x):
        """Return interpolated strain function value."""
        # note: the underlying model defines strain in 10^-6 units
        # but this function should return absolute strain
        y = 1e-6*splev(x, self.tck, ext=3)
        return y
