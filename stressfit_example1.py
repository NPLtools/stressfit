#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import stressfit.commands as comm
import stressfit.shapes as shapes
from IPython.display import HTML
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


# # STRESSFIT
# <p>
# <i>Written by:</i> Jan Saroun, Nuclear Physics Institute CAS, Rez, saroun@ujf.cas.cz<br/>
# <i>Date:</i> 03/20/2023<br/>
# <i>Source:</i> <a href='https://github.com/NPLtools/stressfit'>https://github.com/NPLtools/stressfit</a>
# </p>
# <p>
# This script implements common workflow for the data treatment with StressFit. On the input, STRESFIT uses a list of neutron scattering events with associated data (position, wave vectors, weight factor and "as-measured" lattice spacing - $d_{hkl}$}) accessible at given instrument setup. This list describes the instrumental sampling distribution independently on the sample. It can be obtained by ray-tracing simulation of the instrument using appropriate software, such as McStas (<a href='http://mcstas.org'>http://mcstas.org</a>) or SIMRES (<a href='https://github.com/saroun/simres'>https://github.com/saroun/simres</a>). 
# 
# STRESSFIT provides tools for 3D convolution of such a sampling list with the sample model and permits to calculate: 
# 
# - “centre of gravity” and size of the neutron sampling volume as a function of sample position (in 3D),
# - variation of intensity and position of diffraction peaks due to the perturbation of sampling distribution (material boundaries, absorption, composition and texture gradients),
# - “as measured” (smeared) intensity and strain distributions including the pseudo-strain effects,
# - least-squares fit of intrinsic strain and intensity distributions.
# 
# </p><p>
# STRESSFIT enables to model pseudo-strains for several sample shapes such as curved plates, cylinders, spheres, tubes (both single- and multi-channel) and polygonal rods. Neutron attenuation tables for several common materials generated with the help of the NCrystal library (<a href='https://github.com/mctools/ncrystal'>https://github.com/mctools/ncrystal</a>) are provided as a part of the package resources. Currently, STRESSFIT allows for least squares fitting of measured scans for a single strain component. Simultaneous fitting of multiple strain or stress tensor components is envisaged in future versions.
# </p>    
# 
# ## Jupyter viewer
# 
# User interface based on Jupyter notebook widgets can be launched by executing the code:
# <code>  
# import stressfit.ui.notebook as nb  
# ui = nb.UI()  
# ui.display()  
# </code> 
#     
# Script examples with output of STRESSFIT are also available via Jupyter viewer server:
# <p>
# <a href='http://nbviewer.jupyter.org/url/neutron.ujf.cas.cz/restrax/download/stressfit/stressfit_example1.ipynb'>
# Example 1</a>: ECNS2019, Fitting of strain gradient under the inner surface of a tube, test on synthetic data for STRESS-SPEC.
# </p>
# 
# ## Documentation
# <p>
# For more information and use examples, see: <br/>
# <a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/ECRS2018_stressfit.pdf'>ECRS10, 2018, slides</a><br/>
# <a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/saroun_ECNS2019_poster.pdf'>ECNS 2019, poster</a> <br/>
# <a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/stressfit_ECNS2023_poster.pdf'>ECNS 2023, poster</a> <br/>
# </p>

# ## Workspace definition
# 
# Set the directories for your workspace:
# 
# `workspace`: root directory, should be an absolute path  
# `data`: input path for sampling and experimental data  
# `tables`: input path for lookup tables etc.  
# `output`: output path for all results
# 
# Relative paths should exist under the workspace root directory.
# Set input path to `None` for searching in package resources.

# In[ ]:


"""Set the root workspace directory. Use None for the current directory."""
workspace = None

"""Set the other input/output folders (can be absolute or relative).
Relative paths should exist under the workspace directory.
Set to None for searching in package resources.
"""
env = {'output': './output', # output path (./output is the default if None)
       'data': None,         # input path for sampling and experimental data 
       'tables': None        # input path for lookup tables etc.
      }

# Set workspace and validate
comm.set_workspace(workspace)
comm.set_environment(**env)
comm.validate_workspace()


# ## Sample shape
# 
# Sample shape can be created by the function `shape.create(ID, **kwargs)`, which take the shape ID string and other named parameters (kwargs) as the arguments. The ID's are defined as constants in the shapes module. The named parameters depend on the shape type. For example,
# <code> S.create(shapes.Plate, thickness=10.0) </code> defines an infinitely large flat plate of given thickness. Dimensions are in mm. For documentation on other shapes, execute the command `shapes.help()`.

# In[ ]:


shapes.help()


# In[ ]:


# Define a tube-like sample
# sref=0 means that the depth is measured from the outer surface ...
shape_param = {'Rin':4.0, 'Rout':8.0, 'height':50.0, 'sref':1}
# execute
comm.set_shape(shapes.Tube, **shape_param)


# ## Geometry
# 
# The geometry data include four vectors:
# 
# `angles`: Defines sample rotation with respect to the laboratory frame. The laboratory frame is defined by y-axis vertical and z-axis pointing along the incident beam. The orientation is described by YXY Euler angles in deg.  
# `rotctr`: Sample rotation centre (in sample coordinates).  
# `scandir`: Scan direction in sample coordinates (where the gauge moves relative to sample)  
# `scanorig`: Scan origin in sample coordinates. It is the point corresponding to the zero encoder value (x-value in the input data).

# In[ ]:


"""Define scan geometry."""
geom = {
# Scan direction in sample coordinates (where the gauge moves)
'scandir': [0., 0., -1.],
# Sample orientation (Euler angles YXY) in deg
'angles': [135, 0, 0],
# Scan origin (encoder = 0) in sample coordinates
'scanorig':[0, 0, 0],
# Sample rotation centre (sample coordinates)
'rotctr': [0, 0, 0]
}
# execute
comm.set_geometry(geom)


# ## Material
# 
# Define material attenuation, either of:
#     
# - File name: A table with 2 columns (wavelength [A], attenuation [1/cm])
# - Float number: attenuation [1/cm]
# 
# Example data files can be found in the package resources. Other files should be placed in the `tables` directory of the current workspace.

# In[ ]:


comm.set_attenuation('mu_Fe_gamma.dat')


# ## Sampling list
# 
# Provide filename and number of events to load.  
# The file should be placed in the data input directory of the current workspace. Example files are provided in the package resources. The format is a text of 12 columns, including:
# $id, \mathbf{r}, \mathbf{k}_i, \mathbf{k}_f, p, d_{hkl}$, which denote the line number, position [mm], incident and final wave vectors in 1/Ang , weight, and "as measured" lattice spacing in Ang.

# In[ ]:


# Load file with the sampling event list.  
comm.set_sampling(filename='events_S_1mm.dat', nev=10000)
# print properties of the loaded sampling.
comm.get_sampling().print_properties()


# ## Show the experiment geometry
# 
# This command plots the 2D scene with the sample shape and sampling events projection.
# Provide the number of events to plot as the 1st argument.
# Other optional arguments are:
# - `filename` : png file to save the plot
# - `rang` : range of the plotted area in mm
# - `proj` : projetion plane, eg. xz, zy, etc.
#     
# You may need to adapt sample/geometry above to match the required setup.

# In[ ]:


comm.plot_scene(3000, filename='scene.png', rang=[16, 16], proj='xz')


# ## Calculate resolution and pseudostrain
# 
# Provide the scan range and calculate pseudo-strains for it, using the command `report_pseudo_strains()`. This function calculates the "as-measured" strain and intensity by making convolution of the sampling events with the sample, assuming zero intrinsic strain and uniform scattering intensity.
# 
# Provide the `scan_range` and output `filename` as the 1st and 2nd arguments.
# Other optional arguments are:
# - `nev` : number of events to use
# - `intensity` : if true, show also the intensity profile.

# In[ ]:


# Define scan range in mm: minimum, maximum and number of steps.
scan_range = [-9, -3, 31]
# Calculate and plot (use ; to suppress the returned value print)
comm.report_pseudo_strains(scan_range, '', nev=3000, intensity=True);


# Calculate and plot spatial resolution characteristics, using the command `report_resolution()`.
# 
# This function calculates the size and centre of gravity (CoG) of the actual sampling distribution at given scan positions.
# 
# Provide the `scan_range` and output `filename` as the 1st and 2nd arguments.
# Other optional arguments are:
# - `nev` : number of events to use
# - `cog` : if true, plot also the xyz positions of the sampling CoG.

# In[ ]:


comm.report_resolution(scan_range, '', nev=3000, cog=True)


# ## Load input data
# 
# Below, define file names for the input data: measured strain and intensity. They should be in 3-column text format with scan position, value and error. 
# 
# By default, the previously defined geometry and sampling are used.You can associate the data with other geometry/sampling by providing the values as optional keyward arguments, or afterwards by setting corresponding 
# values to the scan variable. 
# 
# For example: <code>scan['scandir'] = [1,0,0]</code>
# 
# At the end, the pseudo-strains and intensities are compared with the loaded data.

# In[ ]:


# Load and set the input data.
strain = 'eps_SS_rad.dat'
intensity = 'int_SS_rad.dat'
scan = comm.load_input(strain, intensity=intensity)
comm.set_scan(scan)   

# Calculate pseudo-strains and intensities and compare with the data.
comm.report_data_comparison(scan, nev=3000)


# ## Fit intensities  - setup model
# 
# Fitting of intensities allows to determine the variation of scattering 
# probability and extinction with scan depth. It can also help to correct 
# for any missfit between the encoder positions (stored in the data file) 
# and true surface position. 
# 
# Note that sometimes these effects can't be distinguished from each other. 
# With a strong variation of scattering probability 
# (e.g. due to a greadient in texture or composition near the surface), 
# it is not possible to reliably determine the surface position and extinction 
# just from the intensity variation. Then some of the parameters must be 
# determined independently and fixed for fitting. 
# On the other hand, it is the product of extinction and scattering probability 
# distributions which affects pseudo-strains, therefore they do not need to be 
# exactly distinguished.
# 
# *NOTE*:
#     
# Scattering probability and strain distributions are defined on the depth scale.
# Definition of the <i>depth</i> depends on the sample shape. For plates, it is 
# the distance to the upper surface (in local coordinates of the sample). 
# For sphers and cylinders, it is the distance from the curved surface 
# (the outer one in the case of hollow shapes). 
# For complex samples like ETubes, the 'depth' is a position inside the sample 
# projected on the scan direction.
# 
# ### Distribution model
# 
# The depth distributions are modelled as a set of points interpolated by 
# splines of selected order (1 to 3). Define below a minimum number of depth 
# and intensity values which gives a satisfactory estimate of the intensity 
# variation. Obviously, the intensity values should be kept constant for 
# homogeneous materials.
# 
# Define the `x`, `y` distribution values and associated fit-flags. For example,
# `fitx=1` means a free x-variable, 0 means fixed.
# 
# `x` = depth values in mm  
# `y` = intrinsic scattering intensity in rel. units  
# 
# In addition, define the method for `interpolation` between the nodes.
# Use one of `natural`,`clamped`, `PCHIP`, `Akima`. See documentation in the lmfit package. 
# 
# ### Scaling
# 
# Set scaling parameters A, B and zc, and corresponding fit-flags. 
# The scaling formula is `y_scaled = A*y(z-zc) + B`.
# 
# ### Fit configuration
# 
# Following parameters control the fitting process:
# 
# `maxiter` = Maximum number of iterations  
# `bootstrap` = Use bootstrap method for estimation of confidence limits?  
# `loops` = Set loops for the number of bootstrap cycles.  
# `areg` = Regularization parameter  
# `runIFit` = Set False to skip the intensity fit
# 
# ### Initial estimate
# The command `run_fit_guess` below executes the simple fit estimate (without the slower convolution process).
# If necessary, adjust the above parameters to improve the initial model, and execute the block below again.

# In[ ]:


# Distribution model
distribution = {
    'x':    [0., 1., 2.5, 3.5, 4],
    'fitx': [0, 0, 0, 0, 0],
    'y':    [1., 1., 1., 1., 1.],
    'fity': [0, 0, 0, 0, 0]
    }
interpolation = 'natural'

# Scaling
scaling = {
    'A':  20000,
    'B':  0,
    'zc': 0.05
    }
f_scaling = {
    'fA':  1,
    'fB':  0,
    'fzc': 1
    }   

# Maximum number of iterations
maxiter = 100
# Use bootstrap method for estimation of confidence limits?
bootstrap = False
# Set loops for the number of bootstrap cycles.
loops = 3
# regularization parameter
areg = 3
# Set False to skip intensity fit
runIFit = True

"""Commands to define the fitting model."""
ifit = comm.define_ifit(scan, distribution, 3000)
comm.define_scaling(ifit, scaling, f_scaling)
ifit.setInterpModel(interpolation)

"""Guess fit."""
if runIFit:
    comm.run_fit_guess(ifit, maxiter=maxiter, ar=areg, outname='')


# ## Fit intensities - run and report results
# 

# In[ ]:


if runIFit:
    comm.run_fit(ifit, maxiter=maxiter, ar=areg, bootstrap=bootstrap, 
                 loops=loops, outname=scan['intfile'])


# ## Fit strain - setup model
# 
# Fitting of strain depth distribution is similar to the above procedure for 
# fitting intensities. The scattering probability distribution determined 
# above will be automatically taken into account in modelling of 
# pseudo-strains below.
# 
# The depth distributions are modelled as a set of points 
# (depth, strain) interpolated by splines of selected order (1 to 3).
# Define below a minimum number of depth and strain values which is required
# to gives a satisfactory estimate of the strain distribution. 
#  
# Define the `x`, `y` distribution values and associated fit-flags. For example,
# `fitx=1` means a free x-variable, 0 means fixed.
# 
# `x` = depth values in mm  
# `y` = strain in 1e-6 units   
# 
# The 1st and last x-value should always be fixed ...
# 
# In addition, define the method for `interpolation` between the nodes.
# Use one of `natural`,`clamped`, `PCHIP`, `Akima`. See documentation in the lmfit package. 
# 
# ### Fit configuration
# 
# Following parameters control the fitting process:
# 
# `maxiter` = Maximum number of iterations  
# `maxguess` = Maximum iterations for guess fit.
# `bootstrap` = Use bootstrap method for estimation of confidence limits?  
# `loops` = Set loops for the number of bootstrap cycles.  
# `areg` = Regularization parameter  
# `aregs` = A list of regularization factors to scan through
# `runSFit` = Set False to skip the strain fit
# `runReg` = Run regularization loop?
# 
# ### Initial estimate
# The command `run_fit_guess` below executes the simple fit estimate (without the slower convolution process).
# If necessary, adjust the above parameters to improve the initial model, and execute the block below again.
# 

# In[ ]:


distribution = {
    'x':   [0., 1., 2.0 , 2.5, 3.0 , 3.5 , 3.7 ,  4.],
    'fitx': [0] + 6*[1] + [0],
    'y':  8*[0],
    'fity': 8*[1]
    }
interpolation = 'natural'

# You can use surface position from the intensity fit:
zc = ifit.params['xc'].value
print('Using surface position: {:g}\n'.format(zc))

"""Commands to define the fitting model."""
sfit = comm.define_sfit(scan, distribution, 3000, z0=zc)
sfit.setInterpModel(interpolation)

# Use bootstrap method for estimation of confidence limits?
bootstrap = True
# Set loops for the number of bootstrap cycles.
loops = 3
# Define a list of regularization factors:
aregs = [1, 2, 3, 4, 5]
# maximum iterations for guess fit
maxguess = 100
# maximum iterations for fit
maxiter = 250
# Run regularization loop?
runReg = True
# Run strain fit?
runSFit = True

# Run guess fit with given parameters (see docs for run_fit_guess)
if runSFit:
    comm.run_fit_guess(sfit, maxiter=maxguess, ar=areg, outname='')


# ### Run regularization loop
# 
# Regularization is necessary to avoid overfitting. It introduces additional penalty for non-smooth solutions (minimizes 2nd derivative of the distribution). Its weight is defined by the `areg` parameter. Usually, the optimum value of `areg` is somewhere bewteen 0 and 5. To acieve an optimal result, a scan with reguilarization parameter values should be done (takes long). For this set above  `runReg=True` and define the list of `values in `aregs`. 
# 

# In[ ]:


if runSFit and runReg:
    comm.run_fit_reg(sfit, maxiter=maxiter, ar=aregs, outname='')


# ##  Fit strain - run final fit
# 
# Choose below the optimum regularization value (areg) and run the final fit.

# In[ ]:


areg = 3
if runSFit:
    comm.run_fit(sfit, maxiter=maxiter, ar=areg, outname=scan['epsfile'], 
                 bootstrap=bootstrap, loops=loops)


# In[ ]:




