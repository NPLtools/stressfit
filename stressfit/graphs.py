# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 08:47:39 2017

@author: Jan Saroun, saroun@ujf.cas.cz
"""
import numpy as np
from matplotlib import pyplot as plt
from numpy.linalg import norm

# set default plot properties
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (6, 4),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large'}
plt.rcParams.update(params)
 
def plotFit(mcfit, title = '', ylabel='', save = False, file = 'fit.png'):
    """Plot intensity data and fits as a function of scan position.
    
    Arguments
    ---------
        mcfit: MCCFit class
            instance of MCCFit with intensity distribution and data assigned
        save: boolean
            if True, save plot as PNG figure
        file: string
            file name for the PNG output
    """
    plt.xlabel('Scan position, mm')
    plt.ylabel(ylabel)
    plt.title(title)
    if (mcfit.data is not None):
        if (mcfit.data.shape[1] > 2):
            err = mcfit.data[:,2]
        else:
            err = np.zeros(mcfit.data.shape[0])
        if (mcfit.fit is None):
            [yfit, fit_err, ypos] = mcfit.getSmearedFnc(mcfit.data[:,0])
            mcfit.fit = np.array([mcfit.data[:,0], yfit, fit_err]).T
            pos_err = np.zeros(ypos.shape[0])
            mcfit.pos = np.array([mcfit.data[:,0], ypos, pos_err]).T
            
        ymin=np.min(np.concatenate((mcfit.fit[:,1],mcfit.data[:,1]), axis=0))
        ymax=np.max(np.concatenate((mcfit.fit[:,1],mcfit.data[:,1]), axis=0))
        dy = ymax - ymin
        if ((ymin>=0) and (ymin<0.1*dy)):
            ymin = 0.
        plt.ylim(ymin-0.1*dy, ymax+0.1*dy)         
        plt.errorbar(mcfit.data[:,0], mcfit.data[:,1], yerr=err, fmt='ko', label='data')            
        plt.errorbar(mcfit.fit[:,0], mcfit.fit[:,1], yerr=mcfit.fit[:,2], fmt='r-', label='fit')
        if (mcfit.guess is not None):
            plt.errorbar(mcfit.guess[:,0], mcfit.guess[:,1], fmt='b--', label='guess')
        plt.legend(loc='best', frameon=False)
        fn = str(file)
        if (save and fn):
            plt.savefig(fn, bbox_inches='tight')
        plt.show()
    
def plotModel(mcfit, title = '', ylabel='', save = False, 
              file = 'model.png'): 
    """Plot distribution model as a function of information depth.
    
    Arguments
    ---------
        mcfit: MCCFit class
            instance of MCCFit with distribution and data assigned
        save: boolean
            if True, save plot as PNG figure
        file: string
            file name for the PNG output
    """
    plt.xlabel('Information depth, mm')
    plt.ylabel(ylabel)
    plt.title(title)
    
    model = mcfit.mod
    intp = mcfit.dist
    if ((model is None) or (intp is None)):
        msg = '{}: distribution not initialized.'
        raise Exception(msg.format(mcfit.__class__.__name__))
    ymin = min(intp[:,1])
    ymax = max(intp[:,1])
    if (ymax == ymin):
        margin = 0.05*ymax
    else:
        margin = 0.05*(ymax-ymin)
    if (margin==0):
        margin = 0.1
    ymin -= margin
    ymax += margin
    plt.xlim(min(model[:,0]), max(model[:,0]))
    plt.ylim(ymin, ymax)
    if (model.shape[1] > 3):
        xerr = model[:,2]
        yerr = model[:,3]
    else:
        xerr = np.zeros(model.shape[0])
        yerr = np.zeros(model.shape[0])
    plt.errorbar(model[:,0], model[:,1], xerr=xerr, yerr=yerr, 
                 fmt='ro', label='model points')
    if (intp is not None):
        if (intp.shape[1] > 2):
            ferr = intp[:,2]
        else:
            ferr = np.zeros(intp.shape[0])
        plt.errorbar(intp[:,0], intp[:,1], yerr=ferr, 
                     fmt='r-', label='interpolation')
    plt.legend(loc='best', frameon=False)
    fn = str(file)
    if (save and fn):
        plt.savefig(fn, bbox_inches='tight')
    plt.show()    
  
    
def plotIntensityFit(mcfit, save = False, file = 'strainFit.png'):
    """  Descendand of plotFit, only sets title and ylabel ...
    """
    plotFit(mcfit,title='Measured vs. model intensity', 
            ylabel='Intensity', save=save, file=str(file))

    
def plotIntensityModel(mcfit, save = False, file = 'strainModel.png'): 
    """  Calls ploModel with correct title and ylabel ...
    """
    plotModel(mcfit,title='Scattering probability distribution', 
              ylabel='Scattering probability', save=save, file=str(file))
    
def plotStrainFit(mcfit, save = False, file = 'strainFit.png'):
    """  Descendand of plotFit, only sets title and ylabel ...
    """
    plotFit(mcfit,title='Measured vs. model strain', 
            ylabel='Strain, $\mu\epsilon$', save=save, file=str(file))

    
def plotStrainModel(mcfit, save = False, file = 'strainModel.png'): 
    """  Calls ploModel with correct title and ylabel ...
    """
    plotModel(mcfit,title='Intrinsic strain distribution', 
              ylabel='Strain, $\mu\epsilon$', save=save, file=str(file))

    
def plotInfoDepth(model, save = False, file = 'infodepth.png'):
    """Plot information depth and width as a function of scan position.
    
    Arguments
    ---------
        model: MCCfit
            fit model, must be initialized before
        xpos: ndarray
            encoder positions
        save: boolean
            if True, save plot as PNG figure
        file: string
            file name for the PNG output
    """
    data = model.infodepth
    if (data is not None):
        fig, ax1 = plt.subplots()
        fig.suptitle('Information depth')
        ax2 = ax1.twinx()   
        ax1.set_xlabel('Scan position, mm')
        ax1.set_ylabel('Information depth,  mm', color='k')
        ax2.set_ylabel('Sampling width,  mm', color='k')
        ln1 = ax1.errorbar(data[:,0], data[:,1], fmt='bo-', label='depth')
        ln2 = ax2.errorbar(data[:,0], data[:,2], fmt='ro-', label='width') 
        lns = (ln1, ln2)
        labs = ('depth', 'width')
        ax1.legend(lns, labs, loc='lower right', frameon=False)
        fn = str(file)
        if (save and fn):
            plt.savefig(fn, bbox_inches='tight')
        plt.show()
        


def plotScene(rang, proj, shape, ki, kf, sdir, sampling, save = False, 
              file='scene.png', arrows=True):

    """ Plot figure with scattering geometry 
    (events, scattering vectors, scan direction, sample contours).
    
    Arguments
    ---------
        rang: touple[2]
            plot area in mm
        proj: int
            projection plane: 0: (z, y); 1: (x, y); 2: (x, z)
        shape: Shape
            Shape object with sample geometry etc.
        ki, kf: ndarray[3]
            incident and exit k-vector
        sdir: ndarray[3]
            scan direction - shows motoin of events in stationary sample
        sampling: :class:`Sampling
            Instance of the Sampling object
        save: boolean
            if True, save plot as PNG figure
        file: string
            file name for the PNG output        
    
    """
    nd = sampling.nev
    [jr, jki, jkf, jp, jd] = sampling.idata[0:5]
    r = sampling.sdata[0:nd, jr:jr+3]-sampling.sctr
    p = sampling.sdata[0:nd, jp]
    dhkl = sampling.sdata[0:nd, jd]
# gr.plotScene(rang, proj, shape, ki, kf, sdir, r-sam.sctr, p, dhkl, save = True, file = outpng)

    nd = r.shape[0]
    meanp = np.mean(p)
    w = 1*p/meanp
    r0 = np.zeros((nd, 3))
    rc = shape.getLocalPos([0., 0., 0.])
    for j in range(nd):
        rloc = shape.getLocalPos(r[j, :])
        r0[j, :] = rloc[:]

    # centre in local coord.
    if (proj == 0):
        ix = 2
        iy = 1
    elif (proj == 2):
        ix = 0
        iy = 1
    else:
        ix = 0
        iy = 2
    xmin = -0.5*rang[0] + rc[ix]
    xmax = +0.5*rang[0] + rc[ix]
    ymin = -0.5*rang[1] + rc[iy]
    ymax = +0.5*rang[1] + rc[iy]


    # make color scale for dhkl 
    d0 = np.sum(p*dhkl)/np.sum(p)  # mean dhkl
    dd = 1e6*(dhkl - d0)/d0  # convert to strain
    dmin = np.min(dd)  
    dmax = np.max(dd)
    # dirty trick to set color scale symmetric around zero
    # moving 1st two events outside range and set dd to limits
    out = np.max(abs(np.array(rang)))
    r0[0:2,] = 2*np.array([out, out, out])
    w[0:2] = 0.
    # reduce scale to improve contrast
    dlim = 0.5*max(abs(dmin), abs(dmax))
    dd = np.minimum(dd, dlim)
    dd = np.maximum(dd, -dlim)
    dd[0] = -dlim
    dd[1] = dlim

    # make space for color legend
    wleg = 1.3
    wfig = 5  # size of the figure
    if (rang[0] >= rang[1]):
        size = (wleg + wfig, wfig*rang[1]/rang[0])
    else:
        size = (wleg+ wfig*rang[0]/rang[1], wfig)
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=size)
  # plt.figure(figsize=size)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    lbl = ["x", "y", "z"]
    ax.set_xlabel(lbl[ix]+", mm")
    ax.set_ylabel(lbl[iy]+", mm")
    
    # plot shape
    shape.plotContours(ax, proj, "#444444", "dashed")
    # plot events                       
    cm = plt.cm.get_cmap('RdYlBu')
    sc = ax.scatter(r0[:, ix], r0[:, iy], s=w, cmap=cm, c=dd)
    plt.colorbar(sc)
    # plot ki, kf arrows and scan direction
    if (arrows):
        arrow_len = 0.25*rang[0]
        arrw = 0.02*rang[0]
        
        #/ax = plt.gca()
        ki0 = arrow_len*shape.getLocalDir(ki)/norm(ki)
        kf0 = arrow_len*shape.getLocalDir(kf)/norm(kf)
        p1 = rc-ki0  # ki origin
        # scan arrow shows motion of the events w.r.t. sample shape !!!
        xd = -arrow_len*np.array(sdir)/norm(sdir)
        # ki
        ax.arrow(p1[ix], p1[iy], ki0[ix], ki0[iy],
                 head_width=0.01*arrw, head_length=0.01*arrw, fc='k', ec='k')
        # kf
        ax.arrow(rc[ix], rc[iy], kf0[ix], kf0[iy],
                 head_width=arrw, head_length=arrw, fc='k', ec='k')
        # scan direction
        ax.arrow(rc[ix], rc[iy], xd[ix], xd[iy],
                 head_width=arrw, head_length=arrw, fc='r', ec='r')
    plt.title("Scattering geometry")
    ax.invert_xaxis()
    fn = str(file)
    if (save and fn):
        plt.savefig(fn, bbox_inches='tight')
    plt.show()

