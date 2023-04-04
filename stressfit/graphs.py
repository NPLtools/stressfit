# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 08:47:39 2017.

@author: Jan Saroun, saroun@ujf.cas.cz
"""
# TODO rewrite 2D scene plot to allow Axes as output and use in plot_collection.
# TODO write function for plotting multiple 2D scenes (for a set of scan geometries)
import numpy as np
from matplotlib import pyplot as plt
from numpy.linalg import norm
import copy

# set default plot properties
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (6, 4),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'large',
          'ytick.labelsize':'large'}
plt.rcParams.update(params)

#%% Creating plot data structures

def _plot_data_template():
    toplot = {}
    toplot['type'] = 'graph'
    toplot['title'] = 'Title'
    toplot['xlabel'] = 'x'
    toplot['ylabel'] = 'y'
    toplot['legend'] = True
    toplot['grid'] = True
    toplot['items'] = []
    return toplot
    
def get_plot_pseudostrain(mcfit):
    """Collect pseudo-strain data for plotting.
    
    Uses `mcfit.infodepth` data, which must be initialized by calling 
    :meth:`mcfit.calInfoDepth`.
    
    Parameters
    ----------
        mcfit : :class:`MCCfit`
            instance of MCCFit with model distribution and data assigned
    Return
    ------
    dict
        Data structure prepared for plotting by :func:`plot_collection`
    
    """
    data = mcfit.infodepth
    result = {}
    sres = _plot_data_template()
    ires = _plot_data_template()
    # strain
    sres['title'] = 'Pseudo-strain'
    sres['xlabel'] = 'Scan position, mm'
    sres['ylabel'] = 'Strain,  1e-6'
    sres['legend'] = False
    item = {}
    item['x'] = data[:,0]
    item['y'] = data[:,4]
    item['yerr'] = data[:,8]
    item['label'] = 'simulation'
    item['fmt'] = 'b-'
    sres['items'] = [item]
    result['strain'] = sres
    # intensity
    ires['title'] = 'Sampling intensity'
    ires['xlabel'] = 'Scan position, mm'
    ires['ylabel'] = 'Intensity, rel. units'
    ires['legend'] = False
    item = {}
    item['x'] = data[:,0]
    item['y'] = data[:,3]
    item['yerr'] = data[:,9]
    item['label'] = 'simulation'
    item['fmt'] = 'b-'
    ires['items'] = [item]
    result['intensity'] = ires
    return result

def get_plot_fit(mcfit, title='', ylabel=''):
    """Collect fit results for plotting.
    
    The type of fit results (intensity or strain) depends
    on the actual type of mcfit (IFit or Sfit). Corresponding
    title and y-axis label must be provided as arguments.
    The x-axis is always nominal scan position. 
    
    Parameters
    ----------
        mcfit : :class:`MCCfit`
            Fit model object with results.
        title : str
            plot tite 
        ylabel : str
            y-axis label
    Return
    ------
    dict
        Data structure prepared for plotting by :func:`plot_results`
    
    """
    toplot = _plot_data_template()
    toplot['title'] = title
    toplot['xlabel'] = 'Scan position, mm'
    toplot['ylabel'] = ylabel        
    if (mcfit.data is not None):
    # exp. data
        item = {}
        item['label'] = 'data'
        item['fmt'] = 'ko'    
        item['x'] = mcfit.data[:,0]
        item['y'] = mcfit.data[:,1]        
        if (mcfit.data.shape[1] > 2):
            item['yerr'] = mcfit.data[:,2]
        else:
            item['yerr'] = np.zeros(mcfit.data.shape[0])
        toplot['items'].append(item)
        ymin=np.min(mcfit.data[:,1])
        ymax=np.max(mcfit.data[:,1])
    # fitted function
        if (mcfit.fit is None):
            [yfit, fit_err, ypos] = mcfit.getSmearedFnc(mcfit.data[:,0])
            xfit = mcfit.data[:,0]
        else:
            xfit = mcfit.fit[:,0]
            yfit = mcfit.fit[:,1]
            fit_err = mcfit.fit[:,2]            
        item = {}
        item['label'] = 'fit'
        item['fmt'] = 'r-'
        item['x'] = xfit
        item['y'] = yfit
        item['yerr'] = fit_err
        toplot['items'].append(item)
        ymin=np.min(np.concatenate((yfit, mcfit.data[:,1]), axis=0))
        ymax=np.max(np.concatenate((yfit, mcfit.data[:,1]), axis=0))
        # set ymin=0 if near zero
        dy = ymax - ymin
        if ((ymin>=0) and (ymin<0.1*dy)):
            ymin = 0.
            toplot['ylim'] = [ymin-0.1*dy, ymax+0.1*dy]
    # guess
        if (mcfit.guess is not None):
            item = {}
            item['label'] = 'guess'
            item['fmt'] = 'b--'
            item['x'] = mcfit.guess[:,0]
            item['y'] = mcfit.guess[:,1]
            toplot['items'].append(item)
    return toplot


def get_plot_model(mcfit, title='', ylabel=''):
    """Collect fitted model distribution for plotting.
    
    Parameters
    ----------
         mcfit : :class:`MCCfit`
            Fit model object with model distribution assigned.
        title : str
            plot tite 
        ylabel : str
            y-axis label

    Return
    ------
    dict
        Data structure prepared for plotting by :func:`plot_results`
        
    """
    modp = mcfit.mod
    intp = mcfit.dist
    if ((mcfit is None) or (intp is None)):
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
    xlim = [min(modp[:,0]), max(modp[:,0])]
    ylim = [ymin, ymax]
    
    if (modp.shape[1] > 3):
        xerr = modp[:,2]
        yerr = modp[:,3]
    else:
        xerr = np.zeros(modp.shape[0])
        yerr = np.zeros(modp.shape[0])
            
    toplot = _plot_data_template()
    toplot['title'] = title
    toplot['xlabel'] = 'Information depth, mm'
    toplot['ylabel'] = ylabel
    toplot['xlim'] = xlim
    toplot['ylim'] = ylim
    item = {}
    item['label'] = 'points'
    item['fmt'] = 'ro'                
    item['x'] = modp[:,0]
    item['y'] = modp[:,1]
    item['xerr'] = xerr
    item['yerr'] = yerr
    toplot['items'].append(item)    
    if (intp is not None):
        if (intp.shape[1] > 2):
            ferr = intp[:,2]
        else:
            ferr = np.zeros(intp.shape[0])
        item = {}
        item['label'] = 'interpolation'
        item['fmt'] = 'r-'
        item['x'] = intp[:,0]
        item['y'] = intp[:,1]
        item['yerr'] = ferr
        toplot['items'].append(item)

    return toplot

def get_plot_depth(mcfit):
    """Collect plot data for information depth.
    
    Parameters
    ----------
        mcfit : :class:`MCCfit`
            Model object initialized with :meth:`calResolution`.
    Return
    ------
    dict
        Data structure prepared for plotting by :func:`plot_twins`
        
    """ 
    data = None
    toplot = _plot_data_template()
    toplot['type'] = 'twinx'
    if mcfit.resolution is not None:
        data = mcfit.resolution
        x = data['x']
        yd = data['pos']
        yw = data['width']
    elif mcfit.infodepth is not None:
        data = mcfit.infodepth
        x = data[:,0]
        yd = data[:,1]
        yw = data[:,2]
    if data:
        toplot['title'] = 'Depth scale'
        toplot['xlabel'] = 'Scan position, mm'
        toplot['ylabel'] = ['Information depth,  mm']
        toplot['ylabel'].append('Sampling width,  mm')
        item = {}
        item['iax'] = 0
        item['fmt'] = 'bo-'
        item['x'] = x
        item['y'] = yd
        item['label'] = 'depth'
        toplot['items'].append(item)
        item = {}
        item['iax'] = 1
        item['fmt'] = 'ro-'
        item['x'] = x
        item['y'] = yw
        item['label'] = 'width'
        toplot['items'].append(item)
    return toplot
    

def get_plot_cog(mcfit):
    """Collect plot data for sampling center of gravity.
    
    Parameters
    ----------
        mcfit : :class:`MCCfit`
            Model object initialized with :meth:`calResolution`.
    Return
    ------
    dict
        Data structure prepared for plotting by :func:`plot_graph`
        
    """ 
    data = None
    toplot = _plot_data_template()
    if mcfit.resolution is None:
        return toplot
    x = mcfit.resolution['x']
    data = mcfit.resolution['ctr']
    toplot['title'] = 'Centre of gravity'
    toplot['xlabel'] = 'Scan position, mm'
    toplot['ylabel'] = 'Centre of gravity, mm'
    
    fmts = ['ro-','bo-','go-']
    labels = ['x','y','z']
    for i in range(3):
        item = {}
        item['fmt'] = fmts[i]
        item['x'] = x
        item['y'] = data[:,i]
        item['label'] = labels[i]
        toplot['items'].append(item)

    toplot['grid'] = {}
    toplot['grid']['minor'] = {'color':'0.7','linestyle': '-'}
    toplot['grid']['major'] = {'color':'0.7','linestyle': '-'}
   
    return toplot


def get_plot_comparison(simdata, expdata):
    """Collect plot data for comparing experimental and simulated data.

    Lists of experimental and simulated data sets are provided as  
    arguments (dict). The keys of simdata and expdata must match in 
    order to be plotted on the same plot. simdata items without corresponding
    expdata item are ignored.

    Parameters
    ----------
    simdata : dict
        Named list of data structures used by :func:`plot_graph`
        Each plot plot_graph input is addressed as `simdata[key][what]`,
        where `what` is either `strain` or `intensity`. The key must match
        corresponding item in expdata.
    expdata : dict
        Experimental data. Named list of data structures returned by 
        :func:`stressfit.ui.config.uiconfig().get_scan()`.

    Returns
    -------
    dict
        Data structure prepared for plotting by :func:`plot_collection`.

    """
    def exp_to_dict(expdata, what='int'):
        # convert exp. data do dict for plotting
        out = _plot_data_template()
        edata = {}
        edata['label'] = 'experiment'
        edata['fmt'] = 'ko'
        out['xlabel'] = 'Scan position, mm'
        if what=='int':
            out['title'] = expdata['intfile']
            out['ylabel'] = 'Intensity, rel. units'
            edata['x'] = expdata['int'][:,0]
            edata['y'] = expdata['int'][:,1]
            edata['yerr'] = expdata['int'][:,2]
        else:
            out['title'] = expdata['epsfile']
            out['ylabel'] = 'Strain,  1e-6'
            edata['x'] = expdata['eps'][:,0]
            edata['y'] = expdata['eps'][:,1]
            edata['yerr'] = expdata['eps'][:,2]
        out['items'].append(edata)
        return out
    
    def scale(data, A=1, B=0, x0=0):
        """Scale intensity data."""
        out = {}
        out['x'] = data['x']-x0
        out['y'] = A*data['y'] + B
        if 'yerr' in data:
            out['yerr'] = A*data['yerr']
        else:
            out['yerr'] = np.zeros(len(data['x']))
        return out
    
    coll = {}
    for key in expdata:
        dexp = expdata[key]
        if key in simdata and simdata[key] is not None:
            dsim = simdata[key]
        else:
            dsim = None
        # strain
        toplot = exp_to_dict(dexp, what='eps')
        if dsim and 'strain' in dsim:
            toplot['items'].extend(dsim['strain']['items'])
        coll[key+'_s'] = toplot            
        # intensity
        if 'int' in dexp and dexp['int'] is not None:
            toplot  = exp_to_dict(dexp, what='int')
            if dsim and 'intensity' in dsim:
                # try to scale intensity to match the data
                for item in dsim['intensity']['items']:
                    A = dexp['int'][:,1].mean()/item['y'].mean()
                    scaled_data = scale(item, A=A, B=0, x0=0)
                    item_scaled = copy.deepcopy(item)
                    item_scaled['x'] = scaled_data['x']
                    item_scaled['y'] = scaled_data['y']
                    item_scaled['yerr'] = scaled_data['yerr']
                    toplot['items'].append(item_scaled)
            coll[key+'_i'] = toplot      
    return coll    

def get_plot_reglog(reglog):
    """Collect data for regularization plot.
    
    Parameters
    ----------
    reglog : list
            List of records returned by :meth:`get_fit_record` of mccfit.   
    Return
    ------
    dict
        Data structure prepared for plotting by :func:`plot_twins`
        
    """ 
    toplot = _plot_data_template()
    toplot['type'] = 'twinx'
    toplot['title'] = 'Regularization scan'
    toplot['xlabel'] = 'a_r = log10(areg)+10'
    toplot['ylabel'] = ['chi2']
    toplot['ylabel'].append('Reg. term value')    
    # collect data
    nr = len(reglog)
    ar = nr*[0]
    chi = nr*[0]
    reg = nr*[0]
    for i in range(nr):
        areg = reglog[i]['areg']
        ar[i] = np.log10(areg)+10
        chi[i] = reglog[i]['chi2']
        reg[i] = reglog[i]['reg']    
    # chi2
    item = {}
    item['iax'] = 0
    item['fmt'] = 'ro-'
    item['x'] = ar
    item['y'] = chi
    item['label'] = 'chi2'
    toplot['items'].append(item)
    # regularization term
    item = {}
    item['iax'] = 1
    item['fmt'] = 'b^-'
    item['x'] = ar
    item['y'] = reg
    item['label'] = 'reg'
    toplot['items'].append(item)    

    toplot['yscale'] = 'log'
    toplot['grid'] = {}
    toplot['grid']['minor'] = {'color':'0.7','linestyle': '-'}
    toplot['grid']['major'] = {'color':'0.7','linestyle': '-'}
   
    return toplot

    
#%% Top level plot functions        

def plot_resolution(mcfit, depth=True, cog=False, inline=True, save=False, 
                    file=''):
    """Plot information about sampling centre vs. scan position.
    
    Parameters
    ----------
        mcfit : :class:`MCCfit`
            Fit model, must be initialized with :meth:`mcfit.calResolution` 
        depth : bool
            Plot information depth vs. scan position
        cog : bool
            Plot gauge centre-of-gravity vs. scan position
        inline : bool
            Plot all in one row (else arrange in a column)
        save : boolean
            Save plot as PNG file
        file : string
            File name for the PNG output
    """
    if mcfit.resolution is None: return
    coll = []
    if depth:
        coll.append(get_plot_depth(mcfit))
    if cog:
        coll.append(get_plot_cog(mcfit))
        
    # nothing to plot ... 
    if len(coll)==0: return
    
    plot_collection(coll, dim=len(coll), title='Sampling positions', 
                    save=save, file=file)    
        
    
def plot_fit(mcfit, what='int', toplot=['fit','model'],
             inline=True, save=False, file=''):
    """Plot fit results for intensity and strain distributions.
    
    Parameters
    ----------
    mcfit : :class:`MCCfit`
        Fit model object.
    what : str
        Either 'int' or 'eps' for intensity and strain, respectively.
    toplot : list
        Which graphs should be plotted. Allowed elements are 'fit' and 'model'
    inline : bool
            Plot all in one row (else arrange in a column)
    save : boolean
       Save plot as PNG.
    file : str
        File name for saved plot.

    """
    if what=='int':
        ylabel1 = 'Intensity, rel. units'
        ylabel2 = 'Scattering probability, rel. units'
        title1 = 'Measured vs. model intensity'
        title2 = 'Scattering probability distribution'
    elif what=='eps':
        ylabel1 = 'Strain,  1e-6'
        ylabel2 = 'Strain,  1e-6'
        title1 = 'Measured vs. model strain'
        title2 = 'Intrinsic strain distribution'
    else:
        raise Exception('Undefined fit type: {}'.format(what))
    
    coll = {}
    dim = 2
    inline = inline
    if 'fit' in toplot:
        coll['fit'] = get_plot_fit(mcfit, title=title1, ylabel=ylabel1)
    if 'model' in toplot:
        coll['model'] = get_plot_model(mcfit, title=title2, ylabel=ylabel2)
    plot_collection(coll, dim=dim, inline=inline, title='',
                    save=save, file=file)      

    
def plot_comparison(simdata, expdata, title='', inline=True, save=False, file=''):
    """Plot comparison of experimental and simulated data.

    Lists of experimental and simulated data sets are provided as  
    arguments (dict). The keys of simdata and expdata must match in 
    order to be plotted on the same plot. simdata items without corresponding
    expdata item are ignored.

    Parameters
    ----------
    simdata : dict
        Named list of data structures used by :func:`plot_graph`
        Each plot plot_graph input is addressed as `simdata[key][what]`,
        where `what` is either `strain` or `intensity`. The key must match
        corresponding item in expdata.
    expdata : dict
        Experimental data. Named list of data structures returned by 
        :func:`stressfit.ui.config.uiconfig().get_scan()`.

    """
    coll = get_plot_comparison(simdata, expdata)
    if len(coll)>0:    
        plot_collection(coll, dim=2, inline=inline, title=title,
                        save=save, file=file)      
    

    
def plot_reglog(reglog, save=False, file=''):
    """Plot results of regularization scan.
    
    Parameters
    ----------
    reglog : list
            List of records returned by :meth:`get_fit_record` of mccfit.
    save : boolean
       Save plot as PNG.
    file : str
        File name for saved plot.
    """
    toplot = get_plot_reglog(reglog)
    plot_twins(toplot, save=save, file=file)
    


#%% General plot functions using dict structures on input

def plot_collection(data, dim = None, inline=True, title='', save=False, file=''):
    """Plot multiple results in a grid of plots.
    
    Parameters
    ----------
        data : dict or list
            A collection of datasets to be passed to one of the plotting 
            functions: :func:`plot_graph`, :func:`plot_twins` 
        dim : list or int
            Optional dimensions of the grid, [nx, ny] 
        inline : bool
            Plot all in one row (else plot one per row).
        title : str
            Optional main title.
        save : boolean
            if True, save plot as PNG figure
        file : string
            file name for the PNG output
    """
    # convert data sets to a list
    if isinstance(data, dict):
        src = list(data.values())
    else:
        src = data
    # define plot grid
    if isinstance(dim, list):
        [nx, ny] = dim
    elif isinstance(dim, int):
        nx = dim
        ny = int(0.9999*len(src)/nx) + 1
    else:
        nf = len(src)
        if nf==0:
            return
        if inline:
            nx = nf
            ny = 1
        else:
            nx = 1
            ny = nf
    xwidth = params['figure.figsize'][0]
    yheight = params['figure.figsize'][1]
    fig, axs = plt.subplots(nrows=ny, ncols=nx, 
                            figsize=(nx*xwidth,ny*yheight))
#    plt.subplots_adjust(wspace=0.4, hspace=0.3)
    if title:
        fig.suptitle(title, fontsize='x-large', weight ='bold')
    if np.size(axs)==1:
        axs = np.array([axs],dtype=object)
    axx = axs.reshape((axs.size,))
    for i in range(len(src)):
        d = src[i]
        
        if d is not None:
            if 'type' in d:
                if d['type'] == 'twinx':
                    plot_twins(d, ax=axx[i])
                elif d['type'] == 'graph':
                    plot_graph(d, ax=axx[i])
            else:
                plot_results(d, ax=axx[i])
    fig.tight_layout(w_pad=2.0)
    fn = str(file)
    if (save and file):
        plt.savefig(fn, bbox_inches='tight')
    plt.show()
    # print info below the figure ...
    if (save and file):
        print('Figure saved in {}'.format(fn))

def plot_results(data, ax=None, save=False, file=''):
    """Plot data on given axis.
    
    Obsolete, provided for backward compatibility only. 
    Use :func:`plot_graph` instead.
    
    Parameters
    ----------
        data : dict
            Includes at least: `title, xlabel, ylabel, x, y`.
            It can include also: `xerr, yerr, args`, where args is a dict 
            with other keyword arguments to be passed to :func:`errorbar`.
             
            It can also include a key `other`, which defines additional 
            data to be shown on the plot. It should include at least: 
            `x, y, fmt, label`. 
            Optionally, It can include also: `xerr, yerr` .
              
        ax : Axis
            If defined, plot on given Axis object, otherwise create its own.
        save : boolean
            If True, save plot as PNG figure.
        file : string
            File name for the PNG output.
    """
    def get_args(d):
        args = {'label':None, 'fmt':'ko-', 'xerr': None, 'yerr':None}
        for key in ['label','fmt','xerr','yerr', 'valid']:
            if key in d:
                args[key] = d[key]
        if 'args' in d:
            args.update(d['args'])
        return args
    items = []
    item = {}
    item['x'] = data['x']
    item['y'] = data['y']
    item['args'] = get_args(data)
    items.append(item)
    if 'other' in data:
        other = data['other']
        if isinstance(other, list):
            dlist = other
        else:
            dlist = [other]
        for d in dlist:
            items.append(d)
    data['items'] = items
    data['grid'] = True
    # redirect to the newer plot_graph
    plot_graph(data, ax=ax, save=save, file=file)
    

def plot_graph(data, ax=None, save=False, file=''):
    """Plot data on a simple graph with x,y axes.

    The data keys of data must include: **xlabel, ylabel, items**. Optional
    keys are: **title, legend, grid**. 
    
    The `items` include data for indivitual data sets to plot. Each
    item must include: **x, y**. 
    
    Optional keys **xerr, yerr, valid, label, fmt** are passed to errobar function as 
    arguments. Additional arguments can be added in a list with the key='args'.
    The `valid` array is used to mark the masked data (valid=0) in gray. 
    
    Pass legend location as the value of data['legend'] (optional).
    
    Parameters
    ----------
        data : dict
            Data to plot.
        ax : Axis
            If defined, plot on given Axis object, otherwise create its own.
        save : boolean
            If True, save plot as PNG figure
        file : string
            File name for the PNG output.
    """
    def get_args(d):
        args = {'label':None, 'fmt':'ko-', 'xerr': None, 'yerr':None}
        for key in ['label','fmt','xerr','yerr', 'valid']:
            if key in d:
                args[key] = d[key]
        if 'args' in d:
            args.update(d['args'])
        return args
    
    def get_slice(d,lim=[0,10]):
        out = copy.deepcopy(d)
        rang = range(*lim)
        for key in ['x','y','xerr','yerr']:
            out[key] = d[key][rang]
        return out
    
    def split_to_segments(d):
        """Split plot to masked and unmasked segments."""
        mask = d['valid']
        # mask margin indices
        idx = [i for i in range(len(mask)-1) if mask[i] != mask[i+1]]
        # number of segments
        out = []
        if len(idx) == 0:
            d['masked'] = False
            out.append(d)
        else:
            i0 = 0
            idx.append(len(mask)-1)
            for iseg in idx:
                slc = get_slice(d, lim=[i0, iseg+1])
                if mask[i0]:
                    slc['masked'] = False                    
                else:
                    slc['masked'] = True
                out.append(slc)
                i0 = iseg+1
        return out
    
    def plot(ax, d):
        """Plot data on given axis.
        
        If data includes data mask, plot separately valid data and then
        the masked data in shadow collor.
        """
        mcolor = 'grey'
        mstyle = '-'
        args = get_args(d)
        if 'valid' in d:
            # data mask is provided
            # split to slices
            slcs = split_to_segments(d)
            for i in range(len(slcs)):
                s = slcs[i]
                args = get_args(s)
                if s['masked']:
                    del args['label']
                    del args['fmt']
                    ax.errorbar(s['x'], s['y'], color=mcolor, 
                        linestyle=mstyle, **args)
                else:
                    if i != 0:
                        del args['label']
                    ax.errorbar(s['x'], s['y'], **args)
        else:
            ax.errorbar(d['x'], d['y'], **args)
    
    if data is None:
        return
    if not 'items' in data:
        return
    if ax is None:
        fig, ax1 = plt.subplots() 
    else:
        ax1 = ax
    ax1.set_title(data['title'])
    ax1.set_xlabel(data['xlabel'])
    ax1.set_ylabel(data['ylabel'])
    for d in data['items']:
        plot(ax1, d)
    
    if 'grid' in data:
        if isinstance(data['grid'], dict):
            if 'major' in data['grid']:
                if 'minor' in data['grid']:
                    ax1.minorticks_on()
                    ax1.grid(True, which='minor', **data['grid']['minor'])
                ax1.grid(True, which='major', **data['grid']['major'])
            else:
                ax1.grid(**data['grid'])
        else:
            ax1.grid(True)
    if 'legend' in data:
        if isinstance(data['legend'], str):
            loc = data['legend']
        else:
            loc = 'best'
        ax1.legend(loc=loc, frameon=True)
    if 'xlim' in data:
        ax1.set_xlim(data['xlim'])
    if 'ylim' in data:
        ax1.set_ylim(data['ylim'])
        
    if 'xscale' in data:    
        ax1.set_xscale(data['xscale']) 
    if 'yscale' in data:    
        ax1.set_yscale(data['yscale']) 
    
    # don't show/save if it makes part of a multi-cell plot 
    if ax is None:
        fn = str(file)
        if (save and fn):
            plt.savefig(fn, bbox_inches='tight')
        plt.show()

def plot_twins(data, ax=None, save=False, file=''):
    """Plot data on a graph with twin y-axes.
    
    The data keys of data must include: **xlabel, ylabel, items**. Optional
    keys are: **title, legend**. ylabel is a list with the 2 strings 
    for the y-axes labels.
    
    The `items` include data for indivitual data sets to plot. Each
    item must include: **x, y, iax**, where `iax` is the index (0,1) 
    to the corresponding y-axis.  
    
    Optional keys **xerr, yerr, label, fmt** are passed to errobar function as 
    arguments. Additional arguments can be added in a list with the key='args'.
    
    Pass legend location as the value of data['legend'] (optional).

    Parameters
    ----------
        data : dict
            Data to plot.
        ax : Axis
            If defined, plot on given Axis object, otherwise create its own.
        save : boolean
            If True, save plot as PNG figure
        file : string
            File name for the PNG output.
    """
    def get_args(d):
        """Extract arguments for the errorbar function."""
        args = {'label':None, 'fmt':'ko-', 'xerr': None, 'yerr':None}
        for key in ['label','fmt','xerr','yerr']:
            if key in d:
                args[key] = d[key]
        if 'args' in d:
            args.update(d['args'])
        return args

    if ax is None:
        fig, ax1 = plt.subplots()
    else:
        ax1 = ax
    
    ax2 = ax1.twinx()   
    axes = (ax1, ax2)
    if 'title' in data:
        ax1.set_title(data['title'])
    ax1.set_xlabel(data['xlabel'])
    ax1.set_ylabel(data['ylabel'][0], color='k')
    ax2.set_ylabel(data['ylabel'][1], color='k')
    
    lns = []
    labs = []
    for item in data['items']:
        iax = item['iax']
        if iax in [0,1]:
            args = get_args(item)
            ln = axes[iax].errorbar(item['x'], item['y'], **args)
            lns.append(ln)
            if 'label' in args:
                labs.append(args['label'])
            else:
                labs.append(args[''])
    if 'xlim' in data:
        ax1.set_xlim(data['xlim'])
    if 'xscale' in data:    
        ax1.set_xscale(data['xscale'])         
    if 'ylim' in data:
        for iax in [0,1]:
            ax[iax].set_ylim(data['ylim'])
    if 'yscale' in data: 
        for iax in [0,1]:
            axes[iax].set_yscale(data['yscale']) 
    if 'legend' in data:
        if isinstance(data['legend'], str):
            loc = data['legend']
        else:
            loc = 'best'
        ax1.legend(lns, labs, loc=loc, frameon=True)

    if 'grid' in data:
        if isinstance(data['grid'], dict):
            if 'major' in data['grid']:
                if 'minor' in data['grid']:
                    ax1.minorticks_on()
                    ax1.grid(True, which='minor', **data['grid']['minor'])
                ax1.grid(True, which='major', **data['grid']['major'])
            else:
                ax1.grid(**data['grid'])
        else:
            ax1.grid(True)
    
    if ax is None:
        fn = str(file)
        if (save and fn):
            plt.savefig(fn, bbox_inches='tight')
        plt.show()


#%% 2D scene plots

def plotShape(rang, proj, shape, arrows=True, file=''): 
    """Plot sample in 2D projection.
    
    Parameters
    ----------
        rang: touple[2]
            plot area in mm
        proj: int
            projection plane: 0: (z, y); 1: (x, z); 2: (x, y)
        shape: Shape
            Shape object with sample geometry etc.
        file: string
            file name for the PNG output (leave empty to suppress file saving)
    
    """    
    projections = [[2,1], [0,2], [0,1]]
    [ix,iy] = projections[proj]
    # centre in local coord.
    rc = shape.getLocalPos([0., 0., 0.])
    # plot range
    xmin = -0.5*rang[0] + rc[ix]
    xmax = +0.5*rang[0] + rc[ix]
    ymin = -0.5*rang[1] + rc[iy]
    ymax = +0.5*rang[1] + rc[iy]
    
    # make space for color legend
    wleg = 0.0
    wfig = 5  # size of the figure
    if (rang[0] >= rang[1]):
        size = (wleg + wfig, wfig*rang[1]/rang[0])
    else:
        size = (wleg+ wfig*rang[0]/rang[1], wfig)
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=size)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    lbl = ["x", "y", "z"]
    ax.set_xlabel(lbl[ix]+", mm")
    ax.set_ylabel(lbl[iy]+", mm")
    
    # plot shape
    shape.plotContours(ax, proj, "#444444", "dashed")
    
    # plot ki, kf arrows and scan direction
    if (arrows):
        arrow_len = 0.25*rang[0]
        arrw = 0.02*rang[0]
        # scan arrow shows motion of the events w.r.t. sample shape !!!
        xd = arrow_len*np.array(shape._sdir)/norm(shape._sdir)
        # scan direction
        ax.arrow(rc[ix], rc[iy], xd[ix], xd[iy],
                 head_width=arrw, head_length=arrw, fc='r', ec='r')
    plt.title("Sample shape {}".format(shape.shape_type))
    if proj<2:
        ax.invert_xaxis()
    fn = str(file)
    if file:
        plt.savefig(fn, bbox_inches='tight')
    plt.show()
    if file:
        print('Figure saved in {}.'.format(fn))
        
def plotScene(rang, proj, shape, ki, kf, sdir, sampling, save = False, 
              file='scene.png', arrows=True):
    """Plot 2D figure with scattering geometry.
    
    Plot sampling events with scattering vectors, scan direction, 
    sample contours, ...
    
    Parameters
    ----------
        rang: touple[2]
            plot area in mm
        proj: int
            projection plane: 0: (z, y); 1: (x, z); 2: (x, y)
        shape: Shape
            Shape object with sample geometry etc.
        ki, kf: ndarray[3]
            incident and exit k-vector
        sdir: ndarray[3]
            scan direction - shows motion of events in stationary sample
        sampling: :class:`Sampling
            Instance of the Sampling object
        save: boolean
            if True, save plot as PNG figure
        file: string
            file name for the PNG output        
    
    """
    has_sampling = sampling is not None
    
    # centre in local coord.
    rc = shape.getLocalPos([0., 0., 0.])
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
    
    if has_sampling:
        nd = sampling.nev
        [jr, jki, jkf, jp, jd] = sampling.idata[0:5]
        r = sampling.sdata[0:nd, jr:jr+3]-sampling.sctr
        p = sampling.sdata[0:nd, jp]
        dhkl = sampling.sdata[0:nd, jd]
        nd = r.shape[0]
        meanp = np.mean(p)
        w = 1*p/meanp
        r0 = np.zeros((nd, 3))
        for j in range(nd):
            rloc = shape.getLocalPos(r[j, :])
            r0[j, :] = rloc[:]
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
    
    if has_sampling:
        wleg = 1.3
    else:
        wleg = 0.0
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
    if has_sampling:
        # plot events                       
        cm = plt.cm.get_cmap('RdYlBu')
        sc = ax.scatter(r0[:, ix], r0[:, iy], s=w, cmap=cm, c=dd)
        plt.colorbar(sc, label='pseudo-strain, 1e-6')
    
    
    # plot ki, kf arrows and scan direction
    if (arrows):
        arrow_len = 0.25*rang[0]
        arrw = 0.02*rang[0]
        
        #/ax = plt.gca()
        ki0 = arrow_len*shape.getLocalDir(ki)/norm(ki)
        kf0 = arrow_len*shape.getLocalDir(kf)/norm(kf)
        p1 = rc-ki0  # ki origin
        # scan arrow shows motion of the events w.r.t. sample shape !!!
        xd = arrow_len*np.array(sdir)/norm(sdir)
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
    if proj<2:
        ax.invert_xaxis()
    fn = str(file)
    if (save and file):
        plt.savefig(fn, bbox_inches='tight')
    plt.show()
    if (save and file):
        print('Figure saved in {}.'.format(fn))

