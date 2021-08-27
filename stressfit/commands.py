# -*- coding: utf-8 -*-
# Written by: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
"""Top-level functions useful for running stressfit scripts."""

import numpy as np
import stressfit.sample as sam
import stressfit.mccfit as mc
import stressfit.graphs as gr
from stressfit.geometry import Geometry
import stressfit.dataio as dataio
_deg = np.pi/180.
_geom = Geometry()


def load_sampling(filename, columns=[1, 4, 7, 10, 11], maxn=0, verbose=False):
    """Load Monte Carlo events representing the sampling distribution.
    
    The event table contains neutron coordinates, weights and dhkl values.  
    You need to specify column numbers for position, ki and kf vectors, 
    weights and dhkl.
    
    Imported MC events are defined in the laboratory frame, with the origin 
    at the centre of the instrumental gauge volume. Sample and orientation 
    thus defines zero scan position and scan direction in the sample.
    

    Parameters
    ----------
    filename : str
        Input file name.
    columns : list, optional
        Column indices of r[0], ki[0], kf[0], weight and dhkl (starts from 0)
    maxn : int, optional
        Maximum number or records. If zero, take all records in the file.
    verbose: boolean
        If true, print calculated gauge parameters.

    Returns
    -------
    dict
        - data: event data
        - columns: given input parameter
        - nrec: number of records
        - ctr: sampling centre
        - dmean: mean dhkl
    """
    out = {}
    data = np.loadtxt(filename)
    if maxn:
        nrec = min(maxn, data.shape[0])
    else:
        nrec = data.shape[0]
    # Calculate centre of mass of the distribution 
    P = data[:nrec,columns[3]]/np.sum(data[:nrec,columns[3]])
    ctr = np.zeros(3)
    ki = np.zeros(3)
    kf = np.zeros(3)
    for i in range(3):
        ctr[i] = data[:nrec, columns[0] + i].dot(P)
        ki[i] = data[:nrec, columns[1] + i].dot(P)
        kf[i] = data[:nrec, columns[2] + i].dot(P)
    
    wav = 2*np.pi/np.sqrt(ki.dot(ki))
    dmean = data[:nrec, columns[4]].dot(P)
    tth = np.arcsin(wav/2/dmean)*360/np.pi
    if verbose:
        print('Loaded event list with {:d} records'.format(nrec))    
        print('Gauge centre: [{:g}, {:g}, {:g}] '.format(*ctr))
        print('Mean wavelength: {:g}'.format(wav))
        print('2 theta: {:g}'.format(tth))
        print('d0: {:g}\n'.format(dmean))    
    out['data'] = data[:nrec,:]
    out['columns'] = columns
    out['nrec'] = nrec
    out['ctr'] = ctr
    out['dmean'] = dmean
    out['wav'] = wav
    out['tth'] = tth
    out['ki'] = ki
    out['kf'] = kf
    return out



def plot_scene(nev, filename, rang=[30, 30], proj=1):
    """Plot 2D scene with experiment geometry.

    Parameters
    ----------
    nev : int
        Number of sampling points to plot.
    filename : str
        Output file name (*_scene.png will be added)
    rang : array(2)
        with, height of plotted scene in [mm]
    proj : int
        projection plane: 0: (z, y); 1: (x, y); 2: (x, z)

    """
    # retrieve sampling for given number of events
    nevp = min(nev,10000)
    sampling = sam.getSampling(nevp)
    # format output filename
    f = dataio.derive_filename(filename, ext='png', sfx='scene')
    outpng = str(dataio.get_output_file(f))
    # define ki and kf vectors in laboratory frame
    take_off = sampling.src['tth']*_deg # Scattering angle
    ki = np.array([0., 0., 1.])  # default for SIMRES simulations
    kf = Geometry.rotate(ki, 1, take_off) 
    # do plot
    gr.plotScene(rang, proj, sam.shape, ki, kf, _geom.scandir, sampling,  
                 save=True, file=outpng)

def report_pseudo_strains(scan_range, file, nev=3000):
    """Calculate, plot and save calculated pseuostrains and related data.
    
    Parameters
    ----------
    scan_range : [min, max, n]
        Scan range [mm] given as min, max and number of positions. 
        Positions are relative to the scan centre provided in sample geometry.
    file : str
         Output file name (_gauge.dat will be added).
    nev: int, optional
        Number of events to be used for convolution.
    """
    # Initialize model
    model = mc.Sfit(nev=nev, xdir=_geom.scandir)
    
    # choose randomly a subset of sampling events
    sam.shuffleEvents()
    
    # define strain distribution model
    x = np.linspace(scan_range[0],  scan_range[1], num=scan_range[2])
    y = np.zeros(len(x))
    fx = len(x)*[1]
    fy = len(x)*[1]
    model.defDistribution(par=[x, y], vary=[fx, fy], ndim=100, scaled=True)
    
    f = dataio.derive_filename(file, ext='png', sfx='depth')
    filepng = dataio.get_output_file(f)
    model.calInfoDepth(x)
    gr.plotInfoDepth(model, save=True, file=filepng)
    model.saveInfoDepth('', file)


def set_sampling(filename, path=None, nev=3000):
    """Load sampling points and assign them for use in the current script.
    
    Sampling points are a list of events simulated by Monte Carlo ray-tracing
    and saved as a table in a text file. Each row contains neutron coordinates, 
    weights and associated dhkl values.
    You may need to specify column positions for position, ki and kf vectors, 
    weights and dhk (the default matches the format exported by SIMRES ver. 6).
     
    Imported MC events are defined in the laboratory frame, with the origin at 
    the centre of the instrumental gauge volume. Sample and orientation thus 
    defines zero scan position and scan direction in the sample.
    """
    fs = dataio.get_input_file(filename, path=path)
    sampling = load_sampling(fs, maxn=nev, verbose=True)
    sam.setSampling(sampling)
    
def set_geometry(geometry):
    """Set experiment geometry data.
    
    The geometry is provided as a dictionnary with following keys and values:
        
        scandir : array(3), optional
            Scan direction in sample coordinates
        scanorig : array(3), optional
            Sscan origin (encoder = 0) in sample coordinates, in [mm]
        angles : array(3), optional
            Sample orientation (Euler angles YXY), in [deg]
        rotctr : array(3), optional
            Sample rotation centre (sample coordinates), in [mm]
        
    Parameters
    ----------
    geometry : dict
        Experiment geometry data.

    """
    _geom.define(**geometry)
    if sam.shape is not None:
        sam.shape.reset()
        sam.shape.rotate(*list(_geom.angles))
        sam.shape.moveTo(_geom.scanorig)


def set_attenuation(att, path=None):
    """Set beam attenuation parameter.
    
    This is an attenuation coefficient in [1/cm]. It can be defined 
    by a single number, or as a lookup table as a function of wavelength.
    
    Material attenuation: provide either of
    # File name: A table with 2 columns: wavelength [A], attenuation [1/cm]
    # Float number: attenuation [1/cm]:
    
    **NOTE**: the lookup table is searched for in the package resources,
    or in the directory provided as the optional parameter.

    Parameters
    ----------
    att : str or float
        Attenuation coefficients. Float number is interpretted as a 
        wavelength independent value. String is interpreted as a file name
        with lookup table.
    path: str, optional
        Search directory for the lookup table file. 

    Returns
    -------
    None.

    """
    if isinstance(att, float):
        sam.setExtinction(mu=att)
    else:
        sam.setExtinction(table=dataio.load_data(att, kind='tables', 
                                                 path=path))

def set_shape(shape):
    """Set sample shape properties.
    
    Parameters
    ----------
    shape : obj
        Instance of any descendant of the shapeAbstract class.

    To define sample shape, use classes from stressfit.shaps.

    Examples
    --------    
    Infinite plate:
        `stressfit.shapes.ShapePlate(thickness)`
        
        - z is along thickness
    
    Cylinder:
        `stressfit.shapes.ShapeCyl(radius, height)`
        
        - axis || y    
    Sphere:
        `stressfit.shapes.Sph(radius)`
    
    Hollow cylinder:
        `stressfit.shapes.ShapeShellCyl(radius1, radius2, height)`
        
         - axis || y and outer/inner radii = radius1/radius2. 
              
    Hollow sphere:
        `stressfit.shapes.ShapeShell(Rin, Rout)`      
        
    Curved plate:
        `stressfit.shapes.ShapePlateCurved(thickness, length, height, rho1, rho2)`
        
        - z is along thickness, y is along height
        - rho1 = [hor, ver] are curvatures on the front surface
        - rho2 = [hor, ver] are curvatures on the rear surface
   
    """
    sam.shape = shape
    sam.shape.reset()
    sam.shape.rotate(*list(_geom.angles))
    sam.shape.moveTo(_geom.scanorig)


def set_environment(data=None, tables=None, output=None, instruments=None):
    """Set paths for data input and outpout folders.
    
    Calls dataio.set_path() with the same arguments.
    
    By default, the input paths are the package resource directories, 
    the output path is the current directory. 
    Paths with `None` value remain unchanged.
    
    Parameters
    ----------
    data: str
        Input of strain and intensity data.
    tables: str
        Other auxilliary input such as lookup tables or material data.
    instruments: str
        Instrument configuration files.
    output: str
        Output folder.
    
    """
    dataio.set_path(data=data, output=output, tables=tables, instruments=instruments)
    
def load_input(strain, intensity=None,
              path=None,
              scandir=[0, 0, 1], 
              scanorig=[0, 0, 0], 
              rotctr=[0, 0, 0],
              angles=[0, 0, 0],
              sampling=None):
    """Load experimental data and metadata.
    
    Parameters
    ----------
    strain : str
        File name for strain data: position (encoder value), strain, error
    intensity : str, optional
        File name for intensity data: position (encoder value), strain, error.
    path: str, optional
        Optional search path for the input file. If non, the setting
        from set_environment command is used. 
    scandir : list or array, optional
        Scan direction in sample coordinates
    scanorig : list or array, optional
        Sscan origin (encoder = 0) in sample coordinates
    rotctr : list or array, optional
        Sample rotation centre (sample coordinates)
    angles : list or array, optional
        Sample orientation (Euler angles YXY)
    sampling: dict
        Sampling events loaded by the function load_gauge.

    Returns
    -------
    dict
        Input data: keys are derived from argument names.

    """
    scan = {}
    scan['eps'] = dataio.load_data(strain, path=path)
    scan['epsfile'] = strain
    if intensity:
        scan['int'] = dataio.load_data(intensity, path=path)
        scan['intfile'] = intensity
    else:
        scan['int'] = None
        scan['intfile'] = ''
    scan['scandir'] = np.array(scandir)
    scan['scanorig'] = np.array(scanorig)
    scan['rotctr'] = np.array(rotctr)
    scan['angles'] = np.array(angles)
    scan['sampling'] = sampling

    keys =  ['scandir','scanorig','angles','rotctr']
    scan_geometry = {key: scan[key] for key in keys}
    set_geometry(scan_geometry)
    return scan
  

def report_fit(model, filename, **kwargs):
    """Output all fit results for given model object.

    Parameters
    ----------
    model : obj
        Instance of mcfit.MCCFit
    filename : str
        Output filename 
        (output files wil be mangled with appropriate extensions and suffixe).
    **kwargs : 
        Other arguments passed to model.report_fit()
    """
    if filename:
        fn = filename
    else:
        fn = ''
    model.reportFit(file=str(fn), **kwargs)
    
def define_ifit(scan, nodes, nev, **kwargs):
    """Define model for fitting intensity scan.

    Parameters
    ----------
    scan : dict
        Scan properties, as returned by :func:`~stressfit.commands.load_input`.
    nodes : list(4) or array
        [x,y,fx,fy], where x,y are the node coordinates and fx,fy 
        are flags (0|1) marking corresponding fixed parameters. 
    nev : int
        Number of sampling events to be used for convolution.
    **kwargs : 
        Other arguments pased to Ifit constructor.

    Returns
    -------
    ifit : obj
        Instance of :class:`~stressfit.mccfit.Ifit`

    """
    [x,y,fx,fy] = nodes
    ifit = mc.Ifit(nev=nev, xdir=scan['scandir'], **kwargs)
    ifit.data = scan['int']
    # define the intensity distribution, dim=number of points for interpolation 
    ifit.defDistribution([x, y], [fx, fy], ndim=200)
    return ifit


def estimate_eps0(scan, z0, nev, nmin, nmax):
    """Estimate eps0 from the given range of strain data points.
    
    It is assumed that only pseudostrain is observed in the given range:
        - Calculate eps1 as the average of **pseudostrain** in the given range.
        - Calculate eps2 as the average of **strain data** in the given range.
        - eps0 = eps2 - eps1
 
    **Note**

    We must subtract pseudo-strain since it is taken into account 
    by the convolution procedure. eps0 includes only an intrinsic d0 shift
    or an instrumental effect other than the pseudo-strain (e.g. misalignment)  
 
    Parameters
    ----------
    scan : dict
        Scan properties, as returned by :func:`~stressfit.commands.load_input`.
    z0: float
        Depth missfit of between data (encoder) and scan centre.
    nev : int
        Number of sampling events to be used for convolution.
    nmin,nmax : int
        Range of data to calculate eps0 from.
    
    Returns
    -------
    eps0: float
        Strain value in units of [1e-6] to be subtracted from measured data.

    """
    epsdata = scan['eps']
    # calculate pseudo-strain
    tmp = sam.convGauge(epsdata[:,0]-z0, scan['scandir'], 0, nev)[:,4]
    # subtract pseudo-strain from measured data
    dif = epsdata[nmin:nmax+1,1] - 1e6*tmp[nmin:nmax+1]
    # calculate average from the difference
    eps0 = np.average(dif)
    fmt = 'eps0 calculated form points {:d} to {:d}: eps0 = {:g}'
    print(fmt.format(nmin, nmax, eps0))
    return eps0
    

def define_sfit(scan, nodes, nev, z0=0.0, constFnc=None, eps0=False, 
                avgrange=None, **kwargs):
    """Define model for fitting strain scan.

    Parameters
    ----------
    scan : dict
        Scan properties, as returned by :func:`~stressfit.commands.load_input`.
    nodes : list(4) or array
        [x,y,fx,fy], where x,y are the node coordinates and fx,fy 
    nev : int
        Number of sampling events to be used for convolution.
    z0: float
        Depth missfit of between data (encoder) and scan centre.
    constFnc : function
        Optional constraint function. Accepts dict with fit parameters.
    eps0: boolean or float, optional
        Subtract zero strain eps0 from the data. 
        If the value is boolean and True, then eps0 is calculated from the 
        specified range using the function :func:`estimate_eps0`. The function
        parameters nmin, nmax then must be passed as other named arguments.
        They define the data range from which the average intrinsic eps0 
        is determined. 
    avgrange : list(2) of int, optional
        If defined, an average strain is evaluated in the given range of depths.
    **kwargs : 
        Other arguments pased to Sfit constructor.

    Returns
    -------
    ifit : obj
        Instance of :class:`~stressfit.mccfit.Sfit`
    """
    epsdata = scan['eps']
    [x,y,fx,fy] = nodes
    
    e0 = 0.0
    if eps0:
        if isinstance(eps0, bool) and eps0:
            if 'nmin' in kwargs and 'nmax' in kwargs:
                nmin = kwargs.pop('nmin')
                nmax = kwargs.pop('nmax')
                e0 = estimate_eps0(scan, z0, nev, nmin=nmin, nmax=nmax)
        elif isinstance(eps0, float):
            e0 = eps0
    sfit = mc.Sfit(nev=nev, xdir=scan['scandir'], **kwargs)
    sfit.data = epsdata
    sfit.constraint = constFnc
    if avgrange is not None:
        sfit.avgrange = avgrange
    # define strain distribution, dim=number of points for interpolation 
    # par = nodes [x,y] values
    # vary = corresponding flags for fixed (0) or free(1) variables.
    sfit.defDistribution(par=[x, y], vary=[fx, fy], ndim=100, scaled=True)

    # define function scaling (amplitude, strain offset, depth-shift) 
    sfit.defScaling(par=[1., e0, z0], vary=[0, 0, 0])
    return sfit


def run_fit(model, maxiter=100, areg=1e-3, outname=None, guess=False, 
            bootstrap=False, loops=3):
    """Run fit on given model.

    Parameters
    ----------
    model : obj
        Instance of a model class from ``stressfit.mccfit`` (Ifit or Sfit)
    maxiter : int, optional
        Maximum number of iterations. 
    areg : float
        Regularization parameter.
    outname: str, optional
        Filename for output data. A prefix 'a(n)_' will be aded for yeach loop.
        Set to empty string to suppress output except of intermediate plots. 
        Set to None to suppress any output.
    guess : bool, optional
        If True, run first a guess fit. It is fast, but neglects smearing, 
        only subtracts pseudo-strains. 
    bootstrap : bool, optional
        If True, then estimate errors by bootstrap method. The value of `loops`
        defines the number of cycles.
    loops : int, optional
        Number of cycles for the bootstrap method.

    """
    # choose randomly a subset of sampling events
    sam.shuffleEvents()
    
    if guess:
        # Set maxiter=0 if you don't want to fit, just plot the initial model.
        mc.runFit(model, maxiter=maxiter, areg = areg, guess=True)
    
    res = mc.runFit(model, maxiter=maxiter, areg=areg, bootstrap=bootstrap, 
                        loops=loops)
    if not res:
        print('Fit not finished.')
    
    if outname is not None:
        report_fit(model, outname) 
    
def run_fit_guess(model, maxiter=100, areg=1e-3):
    """Run guess fit.
     
    Fast fit which neglects smearing, only subtracts pseudo-strains. 

    Parameters
    ----------
    model : obj
        Instance of a model class from ``stressfit.mccfit`` (Ifit or Sfit)
    maxiter : int, optional
        Maximum number of iterations. 
    areg : array of float
        Regularization parametesr.The length of the aray defines the number of
        fits.
    """
    mc.runFit(model, maxiter=maxiter, areg = areg, guess=True)
    report_fit(model, '')
        
    
def run_fit_reg(model, maxiter=100, areg=[1e-4, 1e-3], outname=None, guess=True):
    """Run regularization loop with fitting of given model.

    Parameters
    ----------
    model : obj
        Instance of a model class from ``stressfit.mccfit`` (Ifit or Sfit)
    maxiter : int, optional
        Maximum number of iterations. 
    areg : array of float
        Regularization parametesr.The length of the aray defines the number of
        fits. 
    outname: str, optional
        Filename for output data. A prefix 'a(n)_' will be aded for yeach loop.
        Set to empty string to suppress output except of intermediate plots. 
        Set to None to suppress any output.
    guess : bool, optional
        If True, run first a guess fit. It is fast, but neglects smearing, 
        only subtracts pseudo-strains. 
    """
    na = len(areg)
    if guess:
        # Set maxiter=0 if you don't want to fit, just plot the initial model.
        ia = min(max(0,int(na/2)),na-1)  
        mc.runFit(model, maxiter=maxiter, areg=areg[ia], guess=True)
        # report_fit(model, '')

    reglog = np.zeros((na,3))
    for ia in range(na):    
        res = mc.runFit(model, maxiter=maxiter, areg=areg[ia])
        reglog[ia] = [areg[ia], model.chi, model.reg]
        if not res:
            print('Fit not finished.')
        ss = '[areg, chi2, reg] = {:g}\t{:g}\t{:g}\n'.format(*reglog[ia,:])
        print(ss)
        sfx = 'a{:g}_'.format(areg[ia])
        if outname:
            report_fit(model, sfx+outname, reglog=reglog) 
        elif outname=='':
            report_fit(model, '') 
    
    # report the table with regularization progress:
    ss = 'areg\tchi2\treg\n'
    for ia in range(na):
        ss += '{:g}\t{:g}\t{:g}\n'.format(*reglog[ia,:])
    print(ss)
