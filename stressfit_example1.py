# coding: utf-8
import numpy as np
import stressfit.shapes as S
import stressfit.mccfit as mc
import stressfit.commands as comm
from IPython.display import HTML
# import warnings
deg = np.pi/180.
HTML('''<script>
code_show=false; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
<b>To hide/show the code blocks, click <a href="javascript:code_toggle()">here</a>.</b>''')

# warnings.simplefilter('error', UserWarning)


# # STRESSFIT
# <p>
# <i>Written by:</i> Jan Saroun, Nuclear Physics Institute CAS, Rez, saroun@ujf.cas.cz<br/>
# <i>Date:</i> 22/05/2018<br/>
# <i>Version:</i> 0.9.2<br/>
# <i>Source:</i> https://github.com/NPLtools/stressfit
# </p>
# <p>
# This script permits to treat pseudo strains in neutron residual strain measurements. Pseudo strains due to partial immersion of sampling volume in the sample, steep intrinsic strain gradients and heterogeneous distribution of scattering probability can be treated. The script employs Monte Carlo convolution method, which uses sampling points produced by Monte Carlo ray-tracing simulation of the diffractometer to model smearing of intrinsic sample properties by instrumental response. Each point has all available information about a scattering event: position, initial and final neutron wave vector, probability and dhkl value associated with this sampling point. Such a sampling distribution can be simulated rather quickly for any instrument setup, for example by the ray-tracing program SIMRES (http://neutron.ujf.cas.cz/restrax). The sampling distribution is then reused to cary out the convolution with any sample shape, position, orientation and strain distribution. This decoupling of MC ray-tracing simulation from the convolution procedure permits to write a code which is rather fast (typically about 1s per one scan at 1 CPU). It is therefore possible to use the MC integration as a part of cost function for least squares fitting.
# </p><p>
# Development of a Python library STRESSFIT aims at providing tools for analysis of residual stress measurements, where pseudo strains can't be avoided and must be taken into account in data processing. Currently, STRESSFIT enables to model pseudo strains for several basic sample shapes (curved plates, full and hollow cylinders and spheres). It also enables least squares fitting of measured scans for a single strain component. Simultaneous fitting of multiple components with a single stress distribution model is envisaged in future versions. 
# </p>
# 
# ## Jupyter viewer
# Examples with output of STRESSFIT are also available via Jupyter viewer server:
# <p>
# <a href='http://nbviewer.jupyter.org/url/neutron.ujf.cas.cz/restrax/download/stressfit/stressfit_example1.ipynb'>
# Example 1</a>: ECNS2019, Fitting of strain gradient under the inner surface of a tube, test on synthetic data for STRESS-SPEC.
# </p>
# 
# ## Documentation
# <p>
# For more information, see: <br/>
# <a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/saroun_MECASENS_2017.pdf'>MECASENS 2017 poster</a><br/>
# <a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/ECRS2018_stressfit.pdf'>ECRS10, 2018, slides</a><br/>
# <a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/saroun_ECNS2019_poster.pdf'>ECNS 2019, poster</a> <br/>
# </p>
# 

#%% USER INPUT

# Define environment: input/output folders (can be absolute or relative)
# Set to non for searching in package resources
input_path = None # input data
output_path = './output' # output data
tables_path = None # lookup tables etc.
sampling_path = None # directory with simulated sampling points 


# ## Sample definition
# Following block defines the <b>sample shape, position and orientation.</b>
# ### Coordinates
# The laboratory frame is defined by y-axis vertical and z-axis pointing along the incident beam.
# The positioning includes rotation by YXY Euler angles ($\omega$, $\chi$, $\phi$) and a linear shift (multiple commands <code>rotate</code> and <code>moveTo</code> are possible to achieve required position).
# 
# Imported MC events are defined in the laboratory frame, with the origin at the centre of the instrumental gauge volume. Sample position and orientation thus define zero scan position and scan direction in the sample.
# 
# ### Sample shape
# Sample dimensions depend on the selected shape. Choose one of the classes defined in the Shapes folder of the package:
# 
# <code> S.ShapePlate(thickness) </code> <br />
# An infinitely large flat plate of given thickness.
# 
# <code> S.ShapePlateCurved(thickness, length, height, [r1x, r1y],  [r2x, r2y]) </code><br />
# A curved plate of given thickness (z), length (x) and height (y). <code>[r1x, r1y]</code> are curvature radii along x and y of the front surcae (z>0). <code>[r2x, r2y]</code> are the radii for the rear surface.
# 
# <code> S.ShapeCyl(radius, height) </code><br />
# A cylindrical shape with axis along y-axis.
# 
# <code> S.ShapeShellCyl(Rin, Rout, height) </code><br /> 
# A hollow cylinder with axis along y-axis. Rin, Rout are the inner and outer radii.
# 
# <code> S.ShapeSph(radius) </code> <br />
# A spherical sample.
# 
# The next block in this template defines a tube: a hollow cylinder with inner radius 4 mm, outer radius 8 mm and height 50 mm. Zero scan position corresponds to the instrumental gauge volume centered at the surface. Measured strain direction is defined by the angle <code>omega</code>.
# 


# Set sample shape (see comment above).
# Dimensions [mm]
radius1 = 4.
radius2 = 8
height = 50.0
shape = S.ShapeShellCyl(radius1, radius2, height)

# Define input data with info about experimental geometry:
# give filename with strain data as the first argument
# intensity = filename with intensity scan
# scandir = scan direction in sample coordinates
# scanorig = scan origin (encoder = 0) in sample coordinates
# rotctr = sample rotation centre (sample coordinates)
# angles = sample orientation (Euler angles YXY) in deg
scan = comm.load_input('eps_SS_hoop.dat', 
                       intensity='int_SS_hoop.dat', 
                       scandir=[0., 0., -1.],
                       scanorig=[0, 0, 0],
                       rotctr=[0,0,0], 
                       angles=[180+45, 0, 0])

# Material attenuation, provide either of
# File name: A table with 2 columns: wavelength [A], attenuation [1/cm]
# Float number: attenuation [1/cm]:
att = 'Fe_mu.dat'

# file with the sampling points  
sampling_file = 'events_S_1mm.dat' 
# number of sampling points to load from the file
nev_load = 3000
# number of sampling points to plot
nev_plot = 3000
# number of sampling points to use in concolution
nev_use = 3000

# 2D plot of experiment geometry:
# scene width,height in [mm]
scene_range = [16, 16]  
# projection plane (zy=0, xz=1, xy=2)
scene_projection = 1  

# ## Initialization commands

# Set environment
comm.set_environment(data=input_path, output=output_path, tables=tables_path)

# Set sampling distribution
comm.set_sampling(sampling_file, path=sampling_path, nev=nev_load)

# Set beam attenuation
comm.set_attenuation(att)    

# Set sample shape
comm.set_shape(shape)

# Set experiment geometry 
comm.set_geometry(scan)   

# Plot experiment geometry
comm.plot_scene(nev_plot, scan['epsfile'], rang=scene_range, proj=scene_projection)


#%% INTENSITY FIT
# Fitting of intensities allows to determine the variation of scattering probability and extinction with scan depth. It can also help to correct for any missfit between the encoder positions (stored in the data file) and true surface position. Note that sometimes these effects can't be distinguished from each other. With a strong variation of scattering probability (e.g. due to a greadient in texture or composition near the surface), it is not possible to reliably determine the surface position and extinction just from the intensity variation. Then some of the parameters must be determined independently and fixed for fitting. On the other hand, it is the product of extinction and scattering probability distributions which affects pseudo-strains, therefore they do not need to be exactly distinguished.
# 
# NOTE:<br/>
# Scattering probability and strain distributions are defined on depth scale, where depth is the shortest distance to the front sample surface. Definition of the <i>front surface</i> depends on the sample shape. For plates, it is the distance to the upper surface (in local coordinates of the sample). For sphers and cylinders, it is the distance from the curved surface (the outer one in the case of hollow shapes).
# 
# ### Define initial values
# 
# The depth distributions are modelled as a set of points interpolated by splines of selected order (1 to 3). Define below a minimum number of depth and intensity values which gives a satisfactory estimate of the intensity variation. Obviously, the intensity values should be kept constant for homogeneous materials.


# MODEL PARAMETERS
#--------------------------------------------------------------------
# Give initial values, followed by flags (0|1), fx=1 means a free variable, 0 for fixed. 

# depth values [mm]
x =  [0., 1., 2.5, 3.5, 4]
fx = [0, 0, 0, 0, 0]

# intrinsic scattering intensity [rel. units]
y = [1., 1.0, 1.0, 1., 1.]
fy = [0, 0, 0, 0, 0]

# Set the method for interpolation between the nodes.
# Use one of 'natural','clamped', 'PCHIP', 'Akima'
interpolation = 'PCHIP' # PCHIP = Piecewise cubic Hermitian interpolation polynomial

# Background and amplitude
A = 33000
B = 0
# fixed (0) or free (1):
fA = 1
fB = 0

# zero encoder position = position of the front surface in encoder scale
# = where is the surcae position in the data table ...
zc = 0.05
# fixed (0) or free (1):
fzc = 1

# Define fit options
# Maximum number of iterations
maxiter = 100
# Use bootstrap method for estimation of confidence limits?
bootstrap = False
# Set loops for the number of bootstrap cycles.
loops = 3
# regularization
areg = 1e-3
# Set False to skip intensity fit
runIFit = True
#--------------------------------------------------------------------

ifit = comm.define_ifit(scan, [x,y,fx,fy], nev_use)

# define scaling (background, amplitude and depth shift). minval,maxval=limits of these parameters.
ifit.defScaling([A, B, zc], [fA, fB, fzc], minval=[0., 0., -np.inf])
# set the interpolation method
ifit.setInterpModel(interpolation)

# Run guess fit with given parameters (see docs for run_fit_guess)
if runIFit:
    comm.run_fit_guess(ifit, maxiter=100, areg=areg)
    
#%% Is the above guess fit OK? Then continue ...

if runIFit:
    comm.run_fit(ifit, maxiter=maxiter, areg=areg, bootstrap=bootstrap, 
                 loops=loops)
    comm.report_fit(ifit, scan['intfile'], plotSampling=True)    


#%% STRAIN FIT

# Fitting of strain depth distribution is similar to the above procedure for fitting intensities. The scattering probability distribution determined above will be automatically taken into account in modelling of pseudo-strains below.
# 
# ### Define initial values
# 
# The depth distributions are modelled as a set of points [depth, $\epsilon$(depth)] interpolated by splines of selected order (1 to 3). Define below a minimum number of depth and strain values which gives a satisfactory estimate of the strain distribution. 

# Define strain depth distribution  
# Give initial values, followed by flags (0|1), fx=1 means a free variable, 0 for fixed.  

# MODEL PARAMETERS
#--------------------------------------------------------------------
# initial depth values [mm]
x =  [0., 1., 2.0 , 2.5, 3.0 , 3.5 , 3.7 ,  4.]
fx = len(x)*[1]
fx[0] = 0
fx[-1] = 0

# initial strain values [1e-6]
y =  len(x)*[0.]
fy = len(y)*[1]
fy[0]=0

# Set the method for interpolation between the nodes.
# Use one of 'natural','clamped', 'PCHIP', 'Akima'
# PCHIP = Piecewise cubic Hermitian interpolation polynomial
interpolation = 'natural'  

# Define a constraint function (optional)
def constraint(params):
    """Constraint function."""
    # constraint example: keeps surface strain at 0 with 50ue tolerance
    dist = mc.params2dist(params)
    y=dist[1,1]
    return (y-0.)/50. 
# Assign constraint to constFnc in order to use it:
# constFnc=constraint
constFnc=None

# Define fit options
# Use bootstrap method for estimation of confidence limits?
bootstrap = True
# Set loops for the number of bootstrap cycles.
loops = 5
# Define a list of regularization factors:
aregs = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
aregs = [1e-7, 1e-6, ]
# maximum iterations for guess fit
maxguess = 100
# maximum iterations for fit
maxiter=200
# Run regularization loop?
runReg = False
# Run strain fit?
runSFit = True
#--------------------------------------------------------------------

# use surface position from intensity fit
zc = ifit.params['xc'].value
print('Using surface position: {:g}\n'.format(zc))

sfit = comm.define_sfit(scan, [x,y,fx,fy], nev_use, z0=zc, constFnc=constFnc)

# define interpolation method
sfit.setInterpModel(interpolation)

# Run guess fit with given parameters (see docs for run_fit_guess)
if runSFit:
    comm.run_fit_guess(sfit, maxiter=maxguess, areg=areg)
    
#%% Is the above guess fit OK? Then continue ...

# Run fit with regularization
if runSFit and runReg:
    comm.run_fit_reg(sfit, maxiter=maxiter, areg=aregs, outname='')

#%% Choose the best areg value and run the final fit
 
areg = 1e-7
if runSFit:
    comm.run_fit(sfit, maxiter=maxiter, areg=areg, outname=scan['epsfile'], 
                 bootstrap=bootstrap, loops=loops)


