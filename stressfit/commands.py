# -*- coding: utf-8 -*-
# Written by: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
"""Top-level functions useful for running stressfit scripts."""

import numpy as np
import stressfit.sample as sam
import stressfit.mccfit as mc
import stressfit.graphs as gr
from stressfit.geometry import Geometry
import stressfit.dataio as dataio
from stressfit.compliance import Compliance 
_deg = np.pi/180.
_geom = Geometry()
_compl = None
_setup = {}

def validate_workspace(verbose=True):
    try:
        dataio.workspace().validate_paths()
        if verbose:
            print('Workspace setting:')
            dataio.workspace().print_info(absolute=True)
    except Exception as e:
        print('ERROR:')
        print(e)

def set_user_input(setup):
    """Set user input data.
    
    Parameters
    ----------
    setup: dict
        User input data as dict.
    """
    global _setup
    _setup.clear()
    _setup.update(setup) 
    

def load_sampling(file='', path=None, nev=None, columns=[1, 4, 7, 10, 11], 
                  **kwargs):
    """Load Monte Carlo events representing the sampling distribution.
    
    Sampling points are a list of events simulated by Monte Carlo ray-tracing
    and saved as a table in a text file. Each row contains neutron coordinates, 
    weights and associated dhkl values.
    You may need to specify column positions for position, ki and kf vectors, 
    weights and dhk (the default matches the format exported by SIMRES ver. 6).
     
    Imported MC events are defined in the laboratory frame, with the origin at 
    the centre of the instrumental gauge volume. Sample and orientation thus 
    defines zero scan position and scan direction in the sample.  

    Parameters
    ----------
    file : str
        Input file name.
    path: str
        Input directory
    columns : list, optional
        Column positions of r[0], ki[0], kf[0], weight and dhkl
        (indexing from 0)
    nev : int, optional
        Maximum number or records. If zero, take all records in the file.

    Returns
    -------
    :obj:`stressfit.sample.Sampling`
    """
    src = {}
    if isinstance(file, dict):
        fname = file['file']
        if 'path' in file:
            fpath = file['path'] 
        else:
            fpath = path
    else:
        fname = file
        fpath = path       
    fs = dataio.get_input_file(fname, path=fpath)
    data = np.loadtxt(fs)
    #print('load_sampling: {}'.format(fs))
    # backwars compatibility, maxn = nev
    if 'maxn' in kwargs:
        nev = kwargs['maxn']
    if nev:
        nrec = min(nev, data.shape[0])
    else:
        nrec = data.shape[0]
    src['file'] = fs.name
    src['path'] = fs.parent.as_posix()
    src['data'] = data[:nrec,:]
    src['columns'] = columns
    sampling = sam.Sampling(src=src)
    return sampling



def plot_scene(nev, scan=None, filename='', rang=[30, 30], proj=1, save=True):
    """Plot 2D scene with experiment geometry.

    Parameters
    ----------
    nev : int
        Number of sampling points to plot.
    scan: dict
        Meta-data for the scan
    filename : str
        Output file name (*_scene.png will be added)
        If not defined, tries to derive irt from scan['epsfile']
    rang : array(2)
        with, height of plotted scene in [mm]
    proj : int
        projection plane: 0: (z, y); 1: (x, z); 2: (x, y)

    """
    # retrieve sampling for given number of events
    nevp = min(nev,10000)
    sampling = sam.getSampling(nevp)
    # format output filename
    if filename:
        f = dataio.derive_filename(filename, ext='png', sfx='scene')
        outpng = str(dataio.get_output_file(f))
    elif scan is not None and 'epsfile' in scan:
        f = dataio.derive_filename(scan['epsfile'], ext='png', sfx='scene')
        outpng = str(dataio.get_output_file(f))
    else:
        outpng = ''
    # define ki and kf vectors in laboratory frame
    take_off = sampling.src['tth']*_deg # Scattering angle
    ki = np.array([0., 0., 1.])  # default for SIMRES simulations
    kf = Geometry.rotate(ki, 1, take_off) 
    # do plot
    gr.plotScene(rang, proj, sam.shape, ki, kf, _geom.scandir, sampling,  
                 save=save, file=outpng)



def cal_pseudo_strains(scan_range, nev=3000, use_int=False):
    """Calculate pseudo-strains and pseudo-intensities.
    
    Parameters
    ----------
    scan_range : [min, max, n]
        Scan range [mm] given as min, max and number of positions. 
        Positions are relative to the scan centre provided in sample geometry.
    nev : int, optional
        Number of events to be used for convolution.
    use_int : bool
        Use previously fitted intensity distribution in calculation 
        of pseudo-strain.
        
    Returns
    -------
    dict
        Results as dict with keys `strain`, `intensity`.
        Each of them contains a dict with results including:
            title, xlabel, ylabel, x, y
        In addition, it returns `model`, a pointer to the Sfit object.
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
    model.calInfoDepth(x, use_int=use_int)
    data = model.infodepth
    result = {}
    sres = {}
    ires = {}
    # strain
    sres['title'] = 'Pseudo strain'
    sres['xlabel'] = 'Scan position, mm'
    sres['ylabel'] = 'Strain,  1e-6'
    sres['x'] = data[:,0]
    sres['y'] = data[:,4]
    result['strain'] = sres
    # intensity
    ires['title'] = 'Pseudo intensity'
    ires['xlabel'] = 'Scan position, mm'
    ires['ylabel'] = 'Intensity, rel. units'
    ires['x'] = data[:,0]
    ires['y'] = data[:,3]
    result['intensity'] = ires
    return result, model
        

def report_pseudo_strains(scan_range, file, 
                          nev=3000, 
                          intensity=False,
                          inline=True, 
                          plot=True, 
                          save=True,
                          use_int=False):
    """Calculate, plot and save calculated pseuostrains and related data.
    
    Parameters
    ----------
    scan_range : [min, max, n]
        Scan range [mm] given as min, max and number of positions. 
        Positions are relative to the scan centre provided in sample geometry.
    file : str
         Output file name. Just base name is used, `_depth.dat` will be added.
    nev : int, optional
        Number of events to be used for convolution.
    intensity : bool
        Plot also intensity vs. position
    inline : bool
        Plot all in one row (else plot intensity below the strains)
    plot : bool
        Show plot
    save: bool
        Save figures and table with results.
    use_int : bool
        Use previously fitted intensity distribution in calculation 
        of pseudo-strain.
    """
    
    data, model = cal_pseudo_strains(scan_range, nev=nev, use_int=use_int)
    if not intensity and 'intensity' in data:
        del data['intensity']                          
    if plot:
        f = dataio.derive_filename(file, ext='png', sfx='deps')
        filepng = dataio.get_output_file(f)
        gr.plot_collection(data, inline=inline, save=save, file=filepng)

    if file and save:
        model.savePseudoStrain('', file, sfx='deps')
    return data

def report_resolution(scan_range, file, 
                      nev=3000, 
                      depths=False,
                      cog=False,
                      inline=True, 
                      plot=True, 
                      save=True):
    """Calculate, plot and save spatial resolution data.
    
    Parameters
    ----------
    scan_range : [min, max, n]
        Scan range [mm] given as min, max and number of positions. 
        Positions are relative to the scan centre provided in sample geometry.
    file : str
         Output file name (_depth.dat will be added).
    nev : int, optional
        Number of events to be used for convolution.
    cog : bool
        Plot also centre of gravity of the sampling.
    inline : bool
        Plot all in one row (else plot intensity below the strains)
    plot : bool
        Show plot
    save: bool
        Save figures and table with results.
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
    model.calResolution(x)
    if plot:
        f = dataio.derive_filename(file, ext='png', sfx='dpos')
        filepng = dataio.get_output_file(f)
        #gr.plotInfoDepth(model, save=save, file=filepng)
        #gr.plotPseudoStrain(model, save=save, file=filepng2)
        gr.plot_resolution(model, depth=True, 
                    cog=cog, 
                    inline=inline,
                    save=save, 
                    file=filepng)
    if file and save:
        model.saveResolution('', file, sfx='dpos')
    
def set_sampling(sampling):
    """Assign sampling events for use by convolution models.
    
    Parameters
    ----------
    sampling: :obj:`stressfit.sample.Sampling`
        Object returend by :func:`load_sampling`
    
    """
    sam.setSampling(sampling)
    
    
def get_sampling():
    """Return currently set sampling object."""
    return sam.getSampling()
    
def set_geometry(geometry):
    """Set experiment geometry data.
    
    The geometry is provided as a dictionary with following keys and values:
        
        scandir : array(3), optional
            Scan direction in sample coordinates
        scanorig : array(3), optional
            Scan origin (encoder = 0) in sample coordinates, in [mm]
        angles : array(3), optional
            Sample orientation (Euler angles YXY), in [deg]
        rotctr : array(3), optional
            Sample rotation centre (sample coordinates), in [mm]
        
    Parameters
    ----------
    geometry : dict
        Experiment geometry data.

    """
    global _geom, _compl
    _geom.define(**geometry)
    if sam.shape is None:
        raise Exception('Sample shape not defined. Use set_shape().')
    sam.shape.reset()
    sam.shape.rotate(*list(_geom.angles))
    sam.shape.moveTo(-_geom.scanorig)

def set_scan(scan):
    """Set scan parameters."""
    global _compl
    if sam.shape is None:
        raise Exception('Sample shape not defined. Use set_shape().')
    
    keys =  ['scandir','scanorig','angles','rotctr']
    scan_geometry = {key: scan[key] for key in keys}
    set_geometry(scan_geometry)
    
    # Make sure the scan has sampling points defined.
    # Sampling assigned to the scan data has preference.
    if scan['sampling'] is not None:
        sam.setSampling(scan['sampling'])
    sampling = sam.getSampling()
    if sampling is None:    
        raise Exception('Sampling points are not defined. Use set_sampling().')    
    
    # Assign sampling to the scan if not yet done.
    if (scan['sampling'] is None):
        scan['sampling'] = sampling
    
    # Set mean Q direction to the compliance object
    if _compl is not None:
        ki = sampling.src['ki']
        kf = sampling.src['kf']
        ki0 = sam.shape.getLocalDir(ki)
        kf0 = sam.shape.getLocalDir(kf)
        Q = kf0-ki0
        _compl.set_q(Q)
    

def set_attenuation(att, path=None):
    """Set beam attenuation parameter.
    
    This is an attenuation coefficient in [1/cm]. It can be defined 
    by a single number, or as a lookup table as a function of wavelength.
    
    Material attenuation: provide either of
    
    array_like : 
        A table with 2 columns: wavelength [A], attenuation [1/cm]
    str : 
        A name of file with the lookup table.
    float: 
        Attenuation coefficient [1/cm]:
    
    **NOTE**: the lookup table file is searched for in the package resources,
    or in the directory provided as the optional parameter.

    Parameters
    ----------
    att : str, float or array_like
        Lookup table, file name for a lookup table, or value for the 
        attenuation coefficient.
    path: str, optional
        Search directory for the lookup table file. 

    Returns
    -------
    None.

    """
    if isinstance(att, float):
        sam.setExtinction(mu=att)
    elif isinstance(att, str):
        sam.setExtinction(table=dataio.load_data(att, kind='tables', 
                                                 path=path))
    else:
        sam.setExtinction(table=att)

def set_shape(shape, **kwargs):
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
    if shape==None:
        sam.shape.update(**kwargs)
    else:
        sam.shape = shape
        sam.shape.reset()
        sam.shape.rotate(*list(_geom.angles))
        sam.shape.moveTo(_geom.scanorig)


def set_workspace(workdir):
    """Set new workspace directory.
    
    Other workspace paths can be set by :func:`set_environment`.
    For example set_environment(output=some_direcotry) will define a new output 
    directory. Relative paths are then interpreted as relative to the 
    workspace directory.

    Parameters
    ----------
    workdir : str
        New workspace directory name.

    """
    dataio.workspace().change_workspace(workdir);

def set_environment(data=None, tables=None, output=None, instruments=None):
    """Set paths for data input and outpout folders.
    
    Calls dataio..workspace().set_path() with the same arguments.
    
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
    wks = dataio.workspace()
    wks.set_paths(data=data, output=output, tables=tables, 
                    instruments=instruments)
    
def set_compliance(**kwargs):
    """Set material compliance as stressfit.compliance.Compliance object.
    
    The keyward arguments should define either of
        - [E, nu, hkl] to define isotropic material
        - [resource, phase, hkl] to define compliance from resource file
    
    Parameters
    ----------
    E: 
        Young modulus in MPa
    nu: 
        Poisson ratio
    hkl:
        hkl indices for relfecting plane
    resource: str
        resource name (file in JSON format)
        See resources/compliance in the package data for examples.
    phase: str
        phase name as defined in the resource file
    
    Return
    ------
    Compliance
        Instance of the :class:`~stressfit.compliance.Compliance` class.
        
    """
    global _compl
    # create isotropic Compliance 
    if all (k in kwargs for k in ("E", "nu", "hkl")):
        S = Compliance.create_isotropic(kwargs['E'],
                                             kwargs['nu'],
                                             kwargs['hkl'])
    # or create Compliance from resources    
    elif all (k in kwargs for k in ("resource", "phase", "hkl")):
        S = Compliance.from_resource(kwargs['resource'])
        S.set_hkl(kwargs['phase'], kwargs['hkl'])
    else:
        raise Exception('Missing or undefined arguments.')
    _compl = S
    return S

def load_input(strain, intensity=None,
              path=None,
              scandir=[0, 0, 1], 
              scanorig=[0, 0, 0], 
              rotctr=[0, 0, 0],
              angles=[0, 0, 0],
              sampling=None,
              verbose=True):
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
        Scan direction in sample coordinates. It defines motion of the sample
        across stationary gauge. 
    scanorig : list or array, optional
        Sscan origin (encoder = 0) in sample coordinates
    rotctr : list or array, optional
        Sample rotation centre (sample coordinates)
    angles : list or array, optional
        Sample orientation (Euler angles YXY)
    sampling: dict
        Sampling events loaded by the function :func:`load_sampling`.
    verbose : bool
        If True, print info on loaded file.

    Returns
    -------
    dict
        Input data: keys are derived from argument names.

    """
    scan = {}
    scan['eps'] = dataio.load_data(strain, path=path, verbose=verbose)
    scan['epsfile'] = strain
    if intensity:
        try:
            scan['int'] = dataio.load_data(intensity, path=path, verbose=verbose)
            scan['intfile'] = intensity
        except:
            scan['int'] = None
            scan['intfile'] = ''
    else:
        scan['int'] = None
        scan['intfile'] = ''
    scan['scandir'] = np.array(scandir)
    scan['scanorig'] = np.array(scanorig)
    scan['rotctr'] = np.array(rotctr)
    scan['angles'] = np.array(angles)
    scan['sampling'] = sampling
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
        are flags (0|1) marking corresponding free parameters. 
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
    set_geometry(scan)
    ifit = mc.Ifit(nev=nev, xdir=_geom.scandir, **kwargs)
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
    set_geometry(scan)
    # calculate pseudo-strain
    tmp = sam.convGauge(epsdata[:,0]-z0,_geom.scandir, 0, nev)[:,4]
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
    set_geometry(scan)
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
    sfit = mc.Sfit(nev=nev, xdir=_geom.scandir, **kwargs)
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
        Regularization parameters.The length of the array defines the number of
        fits.
    """
    mc.runFit(model, maxiter=maxiter, areg = areg, guess=True)
    report_fit(model, '')
        
def run_fit_alt(model, maxiter=100, areg=1e-3, maxc=10):
    """Run guess fit.
     
    Fast fit which neglects smearing, only subtracts pseudo-strains. 

    Parameters
    ----------
    model : obj
        Instance of a model class from ``stressfit.mccfit`` (Ifit or Sfit)
    maxiter : int, optional
        Maximum number of iterations. 
    areg : array of float
        Regularization parameters.The length of the array defines the number of
        fits.
    """
    mc.runFit_alt(model, maxiter=maxiter, maxc=maxc, areg = areg)
    #report_fit(model, '')
            
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
