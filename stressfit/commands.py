# -*- coding: utf-8 -*-
# Written by: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
"""Top-level functions useful for running stressfit scripts."""

import numpy as np
from contextlib import nullcontext
import stressfit.sample as sam
import stressfit.shapes as shapes
import stressfit.mccfit as mc
import stressfit.graphs as gr
from stressfit.geometry import Geometry
import stressfit.dataio as dataio
from stressfit.compliance import Compliance 
_deg = np.pi/180.
_geom = Geometry()
_compl = None

def validate_workspace(verbose=True):
    """Check that all wokrspace paths are valid."""
    try:
        dataio.workspace().validate_paths()
        if verbose:
            dataio.workspace().info(absolute=True)
    except Exception as e:
        print('ERROR:')
        print(e)

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
            path = file['path'] 
    else:
        fname = file
    if 'verbose' in kwargs:
        verb = kwargs['verbose']
    else:
        verb = False
    fs = dataio.get_input_file(fname, path=path)
    #data = np.loadtxt(fs)
    data = dataio.load_data(fname, path=path, verbose = verb)
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
    pnames = {'zy':0, 'yz':0, 'xz': 1, 'zx': 1, 'xy': 2, 'yx': 2}
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
    if isinstance(proj, str):
        if proj in pnames:
            p = pnames[proj]
        else:
            p = pnames['xz']
    else:
        p = proj
    gr.plotScene(rang, p, sam.shape, ki, kf, _geom.scandir, sampling,  
                 save=save, file=outpng)

def cal_smeared_strain(scan_range, eps=None, intensity=None, nev=3000):
    """Calculate smeared strain and intensity distributions.
    
    Parameters
    ----------
    scan_range : [min, max, n]
        Scan range [mm] given as min, max and number of positions. 
        Positions are relative to the scan centre provided in sample geometry.
    eps : dict
        Intrinsic strain distribution [1e-6].
        The dictionary should contain two keys: 'x', 'y' with corresponding
        values to initiate strain distribution in :class:`stressfit.mcfit.Sfit`.
        If None, assume zero intrinsic strain. 
    intensity : dict
        Intrinsic intensity distribution.
        The dictionary should contain two keys: 'x', 'y' with corresponding
        values to initiate strain distribution in :class:`stressfit.mcfit.Ifit`.
        If None, assume uniform intrinsic intensity distribution.
    nev : int, optional
        Number of events to be used for convolution.

        
    Returns
    -------
    dict
        Results as dict with keys `strain`, `intensity`.
        Each of them contains a dict with results, including:
            title, xlabel, ylabel, x, y, pos, ctr
        Strain data also include yerr.
        
        The meaning is:
            - x, y, yerr - strain or intensity function with error
            - pos - centre of gravity position along scan
            - ctr - centre of gravity xyz coordinates
    """
    # Initialize model
    sfit = mc.Sfit(nev=nev, xdir=_geom.scandir)
    ifit = mc.Ifit(nev=nev, xdir=_geom.scandir)
    
    
    # define strain distribution model
    xmin = min(scan_range[0],  scan_range[1]) - 20
    xmax = max(scan_range[0],  scan_range[1]) + 20
    
    if eps is None:
        sx = np.linspace(xmin, xmax, num=3)
        sy = np.zeros(len(sx))
    else:
        sx = eps['x']
        sy = eps['y']
    fsx = len(sx)*[0]
    fsy = len(sx)*[0]
    
    if intensity is not None:
        ix = intensity['x']
        iy = intensity['y']
        fix = len(ix)*[0]
        fiy = len(ix)*[0]
        ifit.defDistribution(par=[ix, iy], vary=[fix, fiy], ndim=100)
        ifit.updateDistribution()
        ifit.fitFinal()
    else:
        mc.intClear()
    
    sfit.defDistribution(par=[sx, sy], vary=[fsx, fsy], ndim=100)
    sfit.updateDistribution()
    sfit.fitFinal()
    
    # choose randomly a subset of sampling events
    sam.shuffleEvents()
    
    # define steps
    x = np.linspace(scan_range[0],  scan_range[1], num=scan_range[2])
    # calculate pseudo-strains and pseudo-intensities
    sfit.calInfoDepth(x, use_int=True)
    data = sfit.infodepth
    
    # strain was provided?
    if eps is not None:
        # yes, do extra convolution for strain function
        [eps, eeps, pos] = sfit.getSmearedFnc(x)
    else:
        # no, only pseudostrain is reported
        [eps, eeps, pos] = [data[:,4], data[:,8], data[:,1]]
        
    
    result = {}
    sres = {}
    ires = {}
    # strain
    sres['title'] = 'Smeared strain'
    sres['xlabel'] = 'Scan position, mm'
    sres['ylabel'] = 'Strain,  1e-6'
    sres['pos'] = pos
    ires['ctr'] = data[:,5:8]
    sres['x'] = x
    sres['y'] = eps
    sres['yerr'] = eeps
    result['strain'] = sres
    # intensity
    ires['title'] = 'Smeared intensity'
    ires['xlabel'] = 'Scan position, mm'
    ires['ylabel'] = 'Intensity, rel. units'
    sres['pos'] = data[:,1]
    ires['ctr'] = data[:,5:8]
    ires['x'] = x
    ires['y'] = data[:,3]
    result['intensity'] = ires
    return result
        

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
    :class:`MCCfit.Sfit`
        Model with initialized with :meth:`calInfoDepth`.
    """
    # Initialize model
    model = mc.Sfit(nev=nev, xdir=_geom.scandir)
    
    # choose randomly a subset of sampling events
    sam.shuffleEvents()
    
    # define strain distribution model
    xd = np.linspace(scan_range[0],  scan_range[1], num=scan_range[2])
    
    xmin = min(scan_range[0],  scan_range[1]) - 20
    xmax = max(scan_range[0],  scan_range[1]) + 20
    x = np.linspace(xmin,  xmax, num=3)
    y = np.zeros(len(x))
    fx = len(x)*[1]
    fy = len(x)*[1]
    model.defDistribution(par=[x, y], vary=[fx, fy], ndim=100, scaled=True)
    model.calInfoDepth(xd, use_int=use_int)
    return model
        
def cal_resolution(scan_range, nev=3000):
    """Calculate spatial resolution data.
    
    Parameters
    ----------
    scan_range : [min, max, n]
        Scan range [mm] given as min, max and number of positions. 
        Positions are relative to the scan centre provided in sample geometry.
    nev : int, optional
        Number of events to be used for convolution.
        
    Returns
    -------
    :class:`MCCfit.Sfit`
        Model object initialized with :meth:`calResolution`.
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
    return model
    

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
        
    Return
    ------
    dict
        Plot data.
    """   
    model = cal_pseudo_strains(scan_range, nev=nev, use_int=use_int)
    data = gr.get_plot_pseudostrain(model)
    
    if not intensity and 'intensity' in data:
        del data['intensity']                          
    if plot:
        if file:
            f = dataio.derive_filename(file, ext='png', sfx='deps')
            filepng = dataio.get_output_file(f)
        else:
            filepng  = ''
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
    sfx='dpos'
    model = cal_resolution(scan_range, nev=nev)
    if plot:
        if file:
            f = dataio.derive_filename(file, ext='png', sfx=sfx)
            filepng = dataio.get_output_file(f)
        else:
            filepng  = None
        gr.plot_resolution(model, depth=True, 
                    cog=cog, 
                    inline=inline,
                    save=save, 
                    file=filepng)
    if file and save:
        model.saveResolution('', file, sfx=sfx)
 

def report_data_comparison(scan, nev=3000, file='', save=False, out=None):
    """Compare scan data with pseudo strains and intensities.
    
    Parameters
    ----------
    scan : dict
        Either a single scan data, or a dict of such scans. 
        The scan must include at least 'eps' and 'sampling' items.
    file : str
         Output file name (_depth.dat will be added).
    nev : int, optional
        Number of events to be used for convolution.
    save: bool
        Save figures and table with results.
    out : obj
        Optional plot context (for use with ipywidgets.Output)
    
    """
    expdata = {}
    simdata = {}
    if 'eps' in scan:
        dlist = {'scan':scan}
    else:
        dlist = scan
    for name in dlist:
        scan = dlist[name]
        x = scan['eps'][:,0]        
        scan_range = [min(x), max(x), 2*len(x)+1]      
        # set geometry and sampling for given scan
        geometry = {k:scan[k] for k in Geometry.input_keys}
        set_geometry(geometry)
        set_sampling(scan['sampling'])
        res = report_pseudo_strains(scan_range, '', 
                                         nev=nev,
                                         intensity=True,
                                         inline=True,
                                         plot=False, 
                                         save=False)
        expdata[name] = scan
        simdata[name] = res
    # do plot
    title = 'Experimental data vs. pseudo-stran'
    if out is not None:
        with out:
            gr.plot_comparison(simdata, expdata, save=save, file=file, 
                               title=title)
    else:
        gr.plot_comparison(simdata, expdata, save=save, file=file,
                           title=title)
    
def set_sampling(sampling=None, filename=None, **kwargs):
    """Assign sampling events for use by convolution models.
    
    If sampling is not provided, attempt to load sampling from a file.
    
    Parameters
    ----------
    sampling : :obj:`stressfit.sample.Sampling`
        Object returend by :func:`load_sampling`
    filename : str
        Filename to load sampling from.
    kwargs : dict
        Arguments passed to load_sampling when loaded from file.
        Search is by default in the workspace input data directory. 
    
    """
    if sampling is None:
        # assume filename
        if filename:
            sampling = load_sampling(file=filename, **kwargs)
    if sampling is not None:
        sam.setSampling(sampling)
    else:
        logger = dataio.logger()
        msg = 'Cannot set sampling, no valid object or filename provided'
        logger.error(msg)
        
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
    sam.shape.set_geometry(_geom.to_dict())

def set_scan(scan):
    """Set scan parameters."""
    global _compl
    if sam.shape is None:
        raise Exception('Sample shape not defined. Use set_shape().')
    
    scan_geometry = {key: scan[key] for key in Geometry.input_keys}
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

def set_shape(shape=None, **kwargs):
    """Define sample shape properties.
    
    If shape is provided, creates the shape object by calling 
    the shape constructor :func:`stressfit.shapes.create`. Otherwise updates 
    the current shape parameters. Then applies on it the current 
    geometry - see :meth:`set_geometry`.
    
    For more information on creation of the shape objects, run 
    :func:`stressfit.shapes.help`.
    
    Parameters
    ----------
    shape : obj
        Name of the shape or the instance of the shape class. 
        If None, only updates the parameters of the sahpe already defined before. 
    
    kwargs : dict
        Other named arguments passed to the shape costructor or parameter 
        setter.
    

    Example
    -------    
    `from stressfit.shapes import Plate, Cylinder, Sphere, Tube # ... `
    
    Infinite plate:
        `set_shape(Plate, thickness=10.0)`
        
        An infinitely large flat plate of given thickness, 
        z is along thickness.
    
    Cylinder:
        `set_shape(Cylinder, radius=4.0, height=30.0)`
        
        A cylindrical shape with axis along y.
        
    Sphere:
        `set_shape(Sphere, radius=8.0)`
        
        A spherical sample.
    
    Hollow cylinder:
        `set_shape(Tube, Rin=4.0, Rout=8.0, height=30.0, ctr=[0,0], sref=1)`
        
         A hollow cylinder with axis along y-axis. 
         
         - Rin, Rout are the inner and outer radii.
         - ctr is the x,z position of the hole centre
         - sref defines the reference surface for depth calculation 
         (0/1 for the inner/outer surface)
              
    Hollow sphere:
        `set_shape(Shell, Rin=4.0, Rout=8.0)` 
                
        Rin, Rout are the inner and outer radii.
        
    Curved plate:
        `set_shape(PlateCurved, thickness, length, height, rho1, rho2)`
        
        - z is along thickness, y is along height
        - rho1 = [hor, ver] are curvatures on the front surface
        - rho2 = [hor, ver] are curvatures on the rear surface
   
    Multichannel tube:
        `set_shape(Tubes, **kwargs)`
        
        A cylinder with multiple coaxial holes. See the stressfit.shapes.Tubes 
        docstring for the constructor arguments.
   
    Polygonal bar:
        `set_shape(PolygonBar, **kwargs)`
        
        A bar with polygonal basis. See the stressfit.shapes.PolygonBar 
        docstring for the constructor arguments. 
   
    File:
        `set_shape(File, filename='filename.json')`
        
        Defines the shape from a file. The file can be created by calling
        the :meth:`save` method of the given shape object. Only the
        abbove shape types can be defined, but this method allows 
        to define more complex shape configurations like multichannel
        tubes from a previously written json file.
   
    """
    sh = None
    if sam.shape == None:
        # create default shape
        sh = shapes.create('Tube')
    if shape==None:
        if sam.shape is not None:
            sam.shape.update(**kwargs)
        else:
            logger = dataio.logger()
            logger.error('Sample shape not defined yet.')
    elif isinstance(shape,str):
        try:
            sh = shapes.create(shape, **kwargs)
        except:
            logger = dataio.logger()
            logger.exception('Unknown shape type required: {}'.format(shape))
    elif isinstance(shape, shapes.ShapeAbstract):
        sh = shape
    else:
        logger = dataio.logger()
        logger.error('Unknown shape type required: {}'.format(shape))
    if sh is not None:
        sam.shape = sh
        sam.shape.reset()
        if _geom is not None:
            sam.shape.set_geometry(_geom.to_dict())

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
    dataio.workspace().change_workspace(root=workdir);

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
              sampling=None,
              verbose=True,
              **kwargs):
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
    sampling: dict
        Sampling events loaded by the function :func:`load_sampling`.
        If None, use the previously set one.
    verbose : bool
        If True, print info on loaded file.
    kwargs : dict
        Optional information about scan geometry.
        
    kwargs
    ------        
    scandir : list or array, optional
        Scan direction in sample coordinates. It defines motion of the sample
        across stationary gauge. 
    scanorig : list or array, optional
        Scan origin (encoder = 0) in sample coordinates
    rotctr : list or array, optional
        Sample rotation centre (sample coordinates)
    angles : list or array, optional
        Sample orientation (Euler angles YXY)
        

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
    
    # update geometry from gloal setting
    if _geom:
        g = _geom.to_dict()
        scan.update(g)
    # then update from arguments
    for g in kwargs:
        if g in Geometry.input_keys:
            scan[g] = kwargs[g]
    # set sampling from argument or use the previously set one
    if sampling:
        scan['sampling'] = sampling
    else:
        scan['sampling'] = get_sampling()
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
    
def define_ifit(scan, nodes, nev, ndim=200, **kwargs):
    """Define model for fitting intensity scan.

    Parameters
    ----------
    scan : dict
        Scan properties, as returned by :func:`~stressfit.commands.load_input`.
    nodes : list(4) or array or dict
        [x,y,fitx,fity], where x,y are the node coordinates and fitx,fity 
        are flags (0|1) marking corresponding free parameters. 
    nev : int
        Number of sampling events to be used for convolution.
    ndim : int
        Number of points for interpolated distribution function.
    **kwargs : 
        Other arguments passed to Ifit constructor.

    Returns
    -------
    ifit : obj
        Instance of :class:`~stressfit.mccfit.Ifit`

    """
    if isinstance(nodes, dict):
        [x,y,fx,fy] = [nodes[k] for k in ['x', 'y', 'fitx', 'fity']]
    else:
        [x,y,fx,fy] = nodes
    set_geometry(scan)
    ifit = mc.Ifit(nev=nev, xdir=_geom.scandir, **kwargs)
    ifit.data = scan['int']
    # define the intensity distribution, dim=number of points for interpolation     
    ifit.defDistribution([x, y], [fx, fy], ndim=ndim)
    return ifit

def define_scaling(model, scale, flags, minval=[0., 0., -np.inf], maxval=None):
    """Define scaling parameters for given model.
    
    The scaling formula of y(z) distribution is
    
    y_scaled = A*y(z-zc) + B

    Parameters
    ----------
    scale : dict ot list
        'A', 'B', 'zc' scaling values.
    flags : dict ot list
        'fA', 'fB', 'fzc' flags, for free (1) or fixed (0) variables. 
    minval : list
        Corresponding lower limits.

    """
    if isinstance(scale, dict):
        [A, B, zc] = [scale[k] for k in ['A', 'B', 'zc']]
    else:
        [A, B, zc] = scale

    if isinstance(flags, dict):
        [fA, fB, fzc] = [flags[k] for k in ['fA', 'fB', 'fzc']]
    else:
        [fA, fB, fzc] = flags
        
    model.defScaling([A, B, zc], [fA, fB, fzc], minval=minval, maxval=maxval)


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
    

def define_sfit(scan, nodes, nev, ndim=200, z0=0.0, eps0=False, constFnc=None, 
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
    ndim : int
        Number of points for interpolated distribution function.        
    z0: float
        Depth missfit between data (encoder) and scan centre.
    eps0: boolean or float, optional
        Subtract zero strain eps0 from the data. 
        If the value is boolean and True, then eps0 is calculated from the 
        specified range using the function :func:`estimate_eps0`. The function
        parameters nmin, nmax then must be passed as other named arguments.
        They define the data range from which the average intrinsic eps0 
        is determined. 
    constFnc : function
        Optional constraint function. Accepts dict with fit parameters.
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
    if isinstance(nodes, dict):
        [x,y,fx,fy] = [nodes[k] for k in ['x', 'y', 'fitx', 'fity']]
    else:
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
    sfit.defDistribution(par=[x, y], vary=[fx, fy], ndim=ndim, scaled=True)

    # define function scaling (amplitude, strain offset, depth-shift) 
    sfit.defScaling(par=[1., e0, z0], vary=[0, 0, 0])
    return sfit

def run_fit_guess(model, maxiter=100, ar=5.0, outname=None, **kwargs):
    """Run guess fit.
     
    Fast fit which neglects smearing, only subtracts pseudo-strains. 

    Parameters
    ----------
    model : obj
        Instance of a model class from ``stressfit.mccfit`` (Ifit or Sfit)
    maxiter : int, optional
        Maximum number of iterations. 
    ar : float
        Regularization parameter. Since version 1.1, the values are defined 
        as log10(a)+10, where a is actual regularization coefficient.
    outname: str, optional
        Filename for output data.
        Set to empty string to suppress output except of plots. 
        Set to None to suppress any output.
    """
    # handle backward compatibility, obsolete areg
    if 'areg' in kwargs:
        a = kwargs['areg']
        ar = np.log10(a)+10
    else:
        a = 10**(ar-10)
    mc.runFit(model, maxiter=maxiter, areg=a, guess=True);
    if outname is not None:
        report_fit(model, outname)
        
def run_fit(model, maxiter=100, ar=5.0, outname=None, guess=False, 
            bootstrap=False, loops=3, **kwargs):
    """Run fit on given model.

    Parameters
    ----------
    model : obj
        Instance of a model class from ``stressfit.mccfit`` (Ifit or Sfit)
    maxiter : int, optional
        Maximum number of iterations. 
    ar : float
        Regularization parameter. Since version 1.1, the values are defined 
        as log10(a)+10, where a is actual regularization coefficient.
    outname: str, optional 
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
    # handle backward compatibility, obsolete areg
    if 'areg' in kwargs:
        a = kwargs['areg']
        ar = np.log10(a)+10
    else:
        a = 10**(ar-10)
    # choose randomly a subset of sampling events
    sam.shuffleEvents()
    
    if guess:
        run_fit_guess(model, maxiter=maxiter, areg=ar, outname=None)    
    mc.runFit(model, maxiter=maxiter, areg=a, bootstrap=bootstrap, 
                        loops=loops);
 
    if outname is not None:
        report_fit(model, outname, use_int=True) 
            
def run_fit_reg(model, maxiter=100, ar=[5.0, 6.0], outname=None, 
                guess=True, **kwargs):
    """Run regularization loop with fitting of given model.

    Parameters
    ----------
    model : obj
        Instance of a model class from ``stressfit.mccfit`` (Ifit or Sfit)
    maxiter : int, optional
        Maximum number of iterations. 
    ar : array of float
        Regularization parameters. The length of the array defines the number of
        fits. Since version 1.1, the values are defined as log10(a)+10, 
        where a is actual regularization coefficient.
    outname: str, optional
        Filename for output data. A prefix 'a(n)_' will be added for each loop.
        Set to empty string to suppress output except of intermediate plots. 
        Set to None to suppress any graphic output.
    guess : bool, optional
        If True, run first a guess fit. It is fast, but neglects smearing, 
        only subtracts pseudo-strains. 
    """
    # handle backward compatibility, obsolete areg
    if 'areg' in kwargs:
        a = kwargs['areg']
        na = len(a)
        ar = na*[0]
        for ia in range(na): 
            ar[ia] = np.log10(a[ia])+10
            
    na = len(ar)
    if guess:
        ia = min(max(0,int(na/2)),na-1)
        run_fit_guess(model, maxiter=maxiter, areg=ar[ia], outname=None)

    logger = dataio.logger()
    logger.clear(what='prog')
    reglog = []   
    lout = logger.output_prog
    
    outp = (outname is not None)
    save = outp and (len(str(outname))>0)
    if lout is None:
        lout = nullcontext()
    with lout:
        for ia in range(na):    
            ss = 'Loop {}/{}, ar={:g}'.format(ia+1, na, ar[ia])
            logger.progress(ss)
            a = 10**(ar[ia]-10)
            mc.runFit(model, maxiter=maxiter, areg=a, clearlog=False);
            rec = model.get_fit_record()
            reglog.append(rec)
            sfx = 'a{:g}_'.format(ar[ia])
            if outp:
                if save:
                    report_fit(model, sfx+outname, reglog=reglog)
                else:
                    report_fit(model, '')
            print(' ')          
        # plot reglog
        if outp:
            if save:
                f = dataio.derive_filename(outname, ext='png', sfx='reglog')
                outpng = dataio.get_output_file(f)
            else:
                outpng = None
            gr.plot_reglog(reglog, save=save, file=outpng)
    # report the table with regularization progress:
    ss = model.format_reglog(reglog)
    logger.progress(ss)
