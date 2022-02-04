# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 23:53:16 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""
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

intpmodels = ['natural','clamped', 'PCHIP', 'Akima']

def intClear():
    global _ispline
    _ispline = None

def intDefined():
    return _ispline is not None

def intFnc(x):
    """Returns interpolated intensity function value """
    if (_ispline is None):
        y = np.ones(x.shape)
    else:
        y = splev(x, _ispline, ext=1)
    return np.maximum(0., y)

def quiet(b):
    global _quiet
    _quiet=b

def path2win(name):
    out = name.replace('\\','/')
    if (out[-1] != '/'):
        out += '/'
    return out

def deriveFilename(file, ext='', sfx=''):
    """Return string derived from given filename.
    
      1)  Remove file extension (if any)
      2)  add suffix 
      3)  add extensio
       
       Arguments
       ---------
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
    """ Convert params to array.
    
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
    """ Convert distribution part of params to array.
    
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
    """ Convert array to params.
    
    NOTE: we rely that params is an ordered dictionary
    
    Arguments
    ---------
        par, err : arrays
            values and std. errors
    
    Returns
    -------
       Parameters
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
    """Get smoothness for regularization:
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
    """Get strain integral over given range. 
    """
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
    """return reduced chi2 from the previously calculated residuals
    """
    ivar = 0
    for key in params:
        p = params[key]
        ivar += p.vary
    chi2 = np.sum(resid**2)/(resid.size - ivar)
    return chi2

    
def fitcallb(params, iter, resid, *fcn_args, **fcn_kws):
    """ Callback to provide fit progress info.
    """
    global _chi2, _reg
    reg = getChi2(resid, params)
    chi = getChi2(fcn_args[0].resid, params)
    if (reg < _reg):
        _reg = reg
        _chi2 = chi
        fmt = 'iter={:d}, chi2={:g}, reg={:g}'
        if (not _quiet): print(fmt.format(iter, chi, reg))


# define objective function for fitting: returns the array to be minimized
def costFnc(params, model, areg, guess):
    """ Residual array to be passed to the minimizer: intensity scan.
    """
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

def runFit_alt(model, maxiter=200, maxc=10, areg=0., bootstrap=False, loops=3):
    """ Execution of least square fit.
    
        Arguments
        ---------
            model : mccfit class
                class handling model parameters
            maxiter : int
                maximum number of function evaluations
            areg: float
                regularization coefficient
            bootstrap: boolean
                run loops to estimate confidence limits
            loops: int
                number of bootstrap loops
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
    ndim = model.dist.shape[0]
    ndata = model.data.shape[0]
    xdata = model.data[:,0]
    _chi2 = np.inf
    _reg = np.inf
    # initialize distribution array
    xdis = model.dist[:,0]
    model.fitInit()
    iterf = None
    if (maxiter>0):
        data_orig = model.data.copy()
        ydata = data_orig[:,1]
        for i in range(maxc):
            minner = Minimizer(costFnc, model.params, 
                               fcn_args=(model, areg, True), 
                               iter_cb=iterf)
            model.result = minner.minimize(method='leastsq', 
                            params=model.params, 
                            ftol=model.ftol, 
                            epsfcn=model.epsfcn, 
                            max_nfev=maxiter)
            model.params = model.result.params
                  
            [yfit, fit_err, ypos] = model.getSmearedFnc(xdata)
            er = np.sqrt(fit_err**2 + model.data[:,2]**2)
            resid = (yfit - data_orig[:,1])/er
            _chi2 = getChi2(resid, model.params)   
            _reg = model.result.redchi
            print('chi2 = {:g}, reg = {:g}'.format(_chi2, _reg))
            delta = ydata - yfit
            model.data[:,1] = ydata + delta
        model.data = data_orig 
    else:
        [yfit, fit_err, ypos] = model.getSmearedFnc(xdata)
    model.areg = areg
    # par =params as array (par, serr)
    model.par = params2array(model.params)
    # fit = array(x, y, yerr) with fitted curve
    model.fit = np.array([xdata, yfit, fit_err]).T
    # pos = array(x, y, yerr), with nominal and information depths
    model.pos = np.array([xdata, ypos, np.zeros(ndata)]).T

    # dist = array(x, y, xerr, yerr) with interpolated distribution
    ydis = model.getDistribution(xdis)
    dis_err = np.zeros(ndim)
    model.dist = np.array([xdis, ydis, dis_err]).T
    model.chi  = _chi2 
    model.reg = _reg
    if (model.avgrange is not None):
        av = getAverage(xdis, ydis, rang = model.avgrange)  
        model.avg = [av, 0.]
    if (not _quiet) : print('runFit finished\n')
    if (model.result is not None):
        res = model.result.success
    else:
        res = True
    model.fitFinal()
    return res

def runFit(model, maxiter=200, areg=0., bootstrap=False, loops=3, guess=False):
    """ Execution of least square fit.
    
        Arguments
        ---------
            model : mccfit class
                class handling model parameters
            maxiter : int
                maximum number of function evaluations
            areg: float
                regularization coefficient
            bootstrap: boolean
                run loops to estimate confidence limits
            loops: int
                number of bootstrap loops
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
    if (not _quiet) and (not guess): print('runFit finished\n')
    if (model.result is not None):
        res = model.result.success
    else:
        res = True
    model.fitFinal()
    return res

def runFitEx(model, maxiter=200, loops=5, areg=0, guess=False):
    """ Execution of least square fit.
    
        Extended version with error estimates by bootstrap method.
    
        Arguments
        ---------
            model : mccfit class
                class handling model parameters
            maxiter : int
                maximum number of function evaluations
            loops : int
                number of loops in the boostrap cycle
            areg : float
                regularization coefficient
        -------
        Returns
        -------
            results : boolean
                    success
    """
    global _chi2, _reg
    if (model.data is None):
        msg = 'Set value of {}.data.'.format(model.__class__.__name__)
        raise Exception('No data defined for fitting. '+msg)
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
        if (not _quiet): print('Loop {:d}: chi2={:g}, reg={:g}\n'.format(it, _chi2, _reg))    
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
        if (not _quiet): print('runFitEx weights: '+str(pval)+'\n')
        res = model.result.success
    else:
        res = True
    model.fitFinal()
    return res


class MCCfit(ABC):
    def __init__(self, nev=1000, xdir=[0., 0., -1.], ftol=1.e-4, epsfcn=None):
        """ Set parameters needed to run the fit procedure.
            
            Arguments
            ----------                
                nev : int
                    number of points from sampling events used in convolution
                xdir : int
                    scan direction (local coordinates)  
                ftol : float
                    tolerance to be passed to the leastsq minimizer
        """
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
        """ Common tasks to do before fitting """
        ivar = 0
        for key in self.params:
            p =self.params[key]
            ivar += p.vary
        self.nvar = ivar
        self.resetParams()
        xdata = self.data[:,0]
        self.calInfoDepth(xdata)
        
    def fitFinal(self):
        """ Common tasks to do after fitting """
        # mod = disttribution model as array (x, y, xerr, yerr)
        self.mod = params2dist(self.params)
        self.setModelPoints(self.mod[:,0], self.mod[:,1])
               
    def defDistribution(self, par, vary, ndim=100, scaled=False):
        """ Define distribution model parameters.
        
            Creates Parameters object required by Minimizer from lmfit.
            
            Arguments
            ---------
            par : list [x, y]
                node positions (x) and values (y)
            vary : list [fx, fy]
                corresponding flags 0|1 for fixed|free parameters
                interpolation order (1..3)
            ndim: int
                number of points for interpolation
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
        """return self.data + gaussian noise 
        
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
        """Reset params to initial values"""
        self.params = self.params_ini
        
    def defScaling(self, par, vary, minval=None, maxval=None):
        """ Define scaling for the distribution model.
        
            Adds items to the Parameters object required by Minimizer from lmfit.
            Call to this routine assumes previous call to defDistribution !
            
            Arguments
            ---------
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
        """ Extract scaling parameters


            
        Returns
        -------
            model parameters and scaling : [ A, B, xc]

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
        """ Convert params to a list of parameter values.
            
        Returns
        -------
            model parameters and scaling : [xy, A, B, xc]
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
        """ Return interpolated distribution.
        """
        [xy, A, B, xc] = self.getParams()
        self.setModelPoints(xy[:,0], xy[:,1])
        if (self.scaled):
            yval = self.intp(xval)
            res = A*yval + B
        else:
            res = self.intp(xval)
        return res

    def updateDistribution(self):
        """ Shortcut to calculate internally stored distribution 
        """
        dis = self.getDistribution(self.dist[:,0])
        self.dist[:,1] = dis[:]
        self.tck = splrep(self.dist[:,0], self.dist[:,1], k=1, s=0)
        return

    def calResid(self, guess=False):
        """ Calculate and return residuals for leastsq fiting.
        """
        [yf, ef, pos] = self.getSmearedFnc(self.data[:,0], guess=guess)
        er = np.sqrt(ef**2 + self.data[:,2]**2)
        self.resid = (yf - self.data[:,1])/er
        return self.resid

    def calDeriv(self):
        """ Calculate 1st derivative of the distribution.
        """
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
        """Format a table with fitted data
        
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
        """Format a table with fitted parameters and chi2
               
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
            print(e)
            msg = 'Can''t evaluate fit result: '
            msg += 'the fit was not completed successfully or wrong result structure.'
            # raise RuntimeError(msg)
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
        ifunc = None
        if use_int:
            ifunc = intFnc
        [A, B, xc] = self.getScaling()
        self.infodepth = sam.convGauge(x-xc, self.xdir, 0, self.nev, 
                                       ifunc=ifunc)

    def calResolution(self, x):
        [A, B, xc] = self.getScaling()
        self.resolution = sam.convResol(x-xc, self.xdir, 0, self.nev)
        
    def setModelPoints(self, x, y):
        """Set x,y points to define the model
        """
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
        """ Set interpolation model as index or string"""
        try:
            if (type(model) is int):
                if (model<0 or model > len(intpmodels)-1):
                    raise Exception
                self.intpmod = model
            else:
                self.intpmod = intpmodels.index(model)
        except:
            print('Invalid interpolation model. Choose string or index:')
            fmt = '{}:\t{}\n'
            for i in range(len(intpmodels)):    
                print(fmt.format(i, intpmodels[i]))
            
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
    def getSmearedFnc(self, x, guess=False):
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
        """ Override the default getDistribution, impose >=0 values
        """
        d = super().getDistribution(xval)
        res = np.maximum(0., d)        
        return res        

    def defDistribution(self, par, vary, ndim=100, scaled=False):
        """ Override the default getDistribution, impose >=0 values
        """
        super().defDistribution(par, vary, ndim=ndim, scaled=scaled)
        for i in range(self.nd):
            self.params['y_'+str(i)].set(min=0.)

    def reportFit(self, outpath='', file='', plotSampling=False, **kwargs):
        if (file):
            saveit=True
            f = dataio.derive_filename(file, ext='png', sfx='fit')
            outpng = dataio.get_output_file(f, path=outpath)
            # outpng = deriveFilename(outpath+file, 'png', 'fit')
            f = dataio.derive_filename(file, ext='png', sfx='model')
            outpngmodel = dataio.get_output_file(f, path=outpath)
            #outpngmodel = deriveFilename(outpath+file, 'png', 'model')
            f = dataio.derive_filename(file, ext='png', sfx='depth')
            outpngdepth = dataio.get_output_file(f, path=outpath)
            # outpngdepth = deriveFilename(outpath+file, 'png', 'depth') 
        else:
            saveit=False
            outpng=''
            outpngmodel=''
            outpngdepth=''
        gr.plotIntensityFit(self, save=saveit, file=outpng)  
        # Plot result, model   
        gr.plotIntensityModel(self, save=saveit, file=outpngmodel)
        # Save results
        self.saveResults(outpath, file)

        # Information depth and sampling width
        if plotSampling:
            self.calInfoDepth(self.data[:,0])
            gr.plotInfoDepth(self, save=saveit, file=outpngdepth)
            self.saveInfoDepth(outpath, file)

    def intFnc(self, x):
        """Returns interpolated intensity function value """
        y = splev(x, self.tck, ext=1)
        # disable negative intensities
        return np.maximum(0., y)
    
    def fitFinal(self):
        """ Common tasks to do after fitting """
        global _ispline
        super().fitFinal()
        _ispline = self.tck


# class for fitting smeared strains
class Sfit(MCCfit):

    def __init__(self, nev=1000, xdir=[0., 0., -1.], ftol=1.e-4, epsfcn=None):
        super().__init__(nev=nev, xdir=xdir, ftol=ftol, epsfcn=epsfcn)
        self.imodel = None
        
    def getSmearedFnc(self, x, guess=False):
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
        # Plot result, fit
        if (file):
            saveit=True
            f = dataio.derive_filename(file, ext='png', sfx='fit')
            outpng = dataio.get_output_file(f, path=outpath)
            # outpng = deriveFilename(outpath+file, 'png', 'fit')
            f = dataio.derive_filename(file, ext='png', sfx='model')
            outpngmodel = dataio.get_output_file(f, path=outpath)
            # outpngmodel = deriveFilename(outpath+file, 'png', 'model')
        else:
            saveit=False
            outpng=''
            outpngmodel=''
        
        gr.plotStrainFit(self, save=saveit, file=outpng)         
        # Plot result, model
        gr.plotStrainModel(self, save=saveit, file=outpngmodel)       
        # Save results
        self.saveResults(outpath, file, reglog=reglog)
        
        # Report averaged value
        if (self.avgrange is not None):
            fmt = 'Integral over x in [{:g}, {:g}]: {:g} +- {:g}\n' 
            if (not _quiet): print(fmt.format(self.avgrange[0], self.avgrange[1], self.avg[0], self.avg[1]))
        
    def intFnc(self,x):
        """Returns interpolated intensity function value """
        if (_ispline is None):
            y = np.ones(x.shape)
        else:
            y = splev(x, _ispline, ext=1)
        return np.maximum(0., y)
    
    def strainFnc(self,x):
        """Returns interpolated intensity function value """
        # note: the underlying model defines strain in 10^-6 units
        # but this function should return absolute strain
        y = 1e-6*splev(x, self.tck, ext=3)
        return y
