#!/usr/bin/env python
# coding: utf-8
"""Template script for standard workflow when using StressFit.

This script executes in sequence the tasks usually needed when fitting the 
neutron strain mapping data. The data are assumed to be the 3-column
text files with scan position [mm], strain [1e-6] or intensity, and errors 
(stdev). 

The script should work as is with the input data provided by the package 
resources. However, it is intended for use as a template for real work. 
Individual steps are explained in the comments below. 

It is expected that the script is executed interactively, section by section, 
from suitable environment such as Spyder or Jupyter notebook. With the default 
values, it should run as a single command

`python -m stressfit_template.py`

Note that the script includes fiting in several loops for regularization and
error estimate, therefore the execution time may be long.

@author: Jan Saroun, saroun@ujf.cas.cz
"""
import stressfit.commands as comm
import stressfit.shapes as shapes

#%% Workspace definition
"""Set the root workspace directory. Use None for the current directory."""
workspace = None
#workspace = r'your workspace root directory'


"""
Set the other input/output folders (can be absolute or relative).
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


#%% Sample definition

"""
Define sample shape. Dimensions are in mm.
See docstring for comm.set_shape. 
"""
shape = shapes.Tube
# sref=0 means that the depth is measured from the outer surface ...
shape_param = {'Rin':4.0, 'Rout':8.0, 'height':50.0, 'sref':1}
# execute
comm.set_shape(shape, **shape_param)

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


"""
Define material attenuation, either of:
    
- File name: A table with 2 columns: wavelength [A], attenuation [1/cm]
- Float number: attenuation [1/cm]

Example data files can be found in the package resources.
"""
comm.set_attenuation('mu_Fe_gamma.dat')


"""
Load file with the sampling event list.

Provide filename and number of events to loaad.
The file should be placed in the data input directory of the current workspace.
"""
comm.set_sampling(filename='events_S_1mm.dat', nev=10000)
comm.get_sampling().print_properties()


#%% Show the experiment geometry
"""
Plot the 2D scene with the sample shape and sampling events projection.
Provide the number of events to plot as the 1st argument.
Other optional arguments are:
    filename : png file to save the plot
    rang : range of the plotted area in mm
    proj : projetion plane, eg. xz, zy, etc.
    
You may need to adapt sample/geometry above to match the required setup.
"""
comm.plot_scene(3000, filename='scene.png', rang=[16, 16], proj='xz')


#%% Calculate resolution and pseudostrain 

"""Define scan range in mm: minimum, maximum and number of steps."""
scan_range = [-9, -3, 31]

"""
Calculate and plot pseudo-strains for given scan range.

This function calculates the "as-measured" strain and intensity
by making convolution of the sampling events with the sample, assuming zero
intrinsic strain and uniform scattering intensity.

Provide the range and output filename as the 1st and 2nd arguments.
Other optional arguments are:
    nev : number of events to use
    intensity : if true, show also the intensity profile.
"""
comm.report_pseudo_strains(scan_range, '', nev=3000, intensity=True)

"""
Calculate and plot spatial resolution characteristics.

This function calculates the size and centre of gravity (CoG)
of the actual sampling distribution at given scan positions.

Provide the range and output filename as the 1st and 2nd arguments.
Other optional arguments are:
    nev : number of events to use
    cog : if true, plot also the xyz positions of the sampling CoG.
"""
comm.report_resolution(scan_range, '', nev=3000, cog=True)
                    

#%% Load input data
"""
Define input data for measured strain and intensity.

3-column text format with scan position, value and error
"""
strain = 'eps_SS_rad.dat'
intensity = 'int_SS_rad.dat'

"""
Load the data and associated meta-data (scan geometry, sampling).

By default, the previously defined geometry and sampling are used.
You can associate the data with other geometry/sampling by providing the 
values as optional keyward arguments, or afterwards by setting corresponding 
values to the scan variable. 

For example: scan['scandir'] = [1,0,0]
"""
scan = comm.load_input(strain, intensity=intensity)
comm.set_scan(scan)   

"""Calculate pseudo-strains and intensities and compare with the data."""
comm.report_data_comparison(scan, nev=3000)
 

#%% Fit intensities - setup model

"""
Fitting of intensities allows to determine the variation of scattering 
probability and extinction with scan depth. It can also help to correct 
for any missfit between the encoder positions (stored in the data file) 
and true surface position. 

Note that sometimes these effects can't be distinguished from each other. 
With a strong variation of scattering probability 
(e.g. due to a greadient in texture or composition near the surface), 
it is not possible to reliably determine the surface position and extinction 
just from the intensity variation. Then some of the parameters must be 
determined independently and fixed for fitting. 
On the other hand, it is the product of extinction and scattering probability 
distributions which affects pseudo-strains, therefore they do not need to be 
exactly distinguished.


NOTE:
    
Scattering probability and strain distributions are defined on the depth scale.
Definition of the <i>depth</i> depends on the sample shape. For plates, it is 
the distance to the upper surface (in local coordinates of the sample). 
For sphers and cylinders, it is the distance from the curved surface 
(the outer one in the case of hollow shapes). 
For complex samples like ETubes, the 'depth' is a position inside the sample 
projected on the scan direction.
"""


"""
Define initial intensity distribution model.

The depth distributions are modelled as a set of points interpolated by 
splines of selected order (1 to 3). Define below a minimum number of depth 
and intensity values which gives a satisfactory estimate of the intensity 
variation. Obviously, the intensity values should be kept constant for 
homogeneous materials.


Define the x,y distribution values and associated fit-flags. For example,
fitx=1 means a free x-variable, 0 means fixed.

x = depth values [mm]
y = intrinsic scattering intensity [rel. units]

"""
distribution = {
    'x':    [0., 1., 2.5, 3.5, 4],
    'fitx': [0, 0, 0, 0, 0],
    'y':    [1., 1., 1., 1., 1.],
    'fity': [0, 0, 0, 0, 0]
    }

"""
Set the method for interpolation between the nodes.
Use one of 'natural','clamped', 'PCHIP', 'Akima'
See documentation in the lmfit package. 
"""
interpolation = 'natural'

"""
Set scaling parameters and corresponding fit-flags.

The scaling formula is y_scaled = A*y(z-zc) + B
"""
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

"""Commands to define the fitting model."""
ifit = comm.define_ifit(scan, distribution, 3000)
comm.define_scaling(ifit, scaling, f_scaling)
ifit.setInterpModel(interpolation)

"""
Define fit options.
"""
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

"""
Guess fit
If necessary, adjust the above parameters to improve the initial model.
The guess fit is the fast method neglecting de-smearing.

If necessary, adjust the above model parameters to improve the initial model.
"""

# Run guess fit with given parameters (see docs for run_fit_guess)
if runIFit:
    comm.run_fit_guess(ifit, maxiter=maxiter, ar=areg, outname='')

#%% Fit intensities - run and report results
    
if runIFit:
    comm.run_fit(ifit, maxiter=maxiter, ar=areg, bootstrap=bootstrap, 
                 loops=loops, outname=scan['intfile'])

#%% Fit strain - setup model

"""
Fitting of strain depth distribution is similar to the above procedure for 
fitting intensities. The scattering probability distribution determined 
above will be automatically taken into account in modelling of 
pseudo-strains below.

The depth distributions are modelled as a set of points 
[depth, strain] interpolated by splines of selected order (1 to 3).
 Define below a minimum number of depth and strain values which is required
 to gives a satisfactory estimate of the strain distribution. 
"""


"""
Define the x, y distribution values and associated fit-flags. For example,
fitx=1 means a free x-variable, 0 means fixed.

x = depth values [mm]
y = strain in 1e-6 units

the 1st and last x-value should always be fixed ...

"""
distribution = {
    'x':   [0., 1., 2.0 , 2.5, 3.0 , 3.5 , 3.7 ,  4.],
    'fitx': [0] + 6*[1] + [0],
    'y':  8*[0],
    'fity': 8*[1]
    }

"""
Set the method for interpolation between the nodes.
Use one of 'natural','clamped', 'PCHIP', 'Akima'
See documentation in the lmfit package. 
"""
interpolation = 'natural'

# use surface position from intensity fit 
zc = ifit.params['xc'].value
print('Using surface position: {:g}\n'.format(zc))

"""Commands to define the fitting model."""
sfit = comm.define_sfit(scan, distribution, 3000, z0=zc)
sfit.setInterpModel(interpolation)

"""
Define fit options.

A bootstrap cycle, if defined, is applied to estimate confidence limit.
A regularization cycle for given number of coefficients is executed if requested. 
"""
# Use bootstrap method for estimation of confidence limits?
bootstrap = True
# Set loops for the number of bootstrap cycles.
loops = 3
# Define a list of regularization factors:
aregs = [1, 2, 3, 4, 5, 6]
# maximum iterations for guess fit
maxguess = 100
# maximum iterations for fit
maxiter = 200
# Run regularization loop?
runReg = False
# Run strain fit?
runSFit = True

"""
Guess fit
If necessary, adjust the above parameters to improve the initial model.
The guess fit is the fast method neglecting de-smearing.

If necessary, adjust the above model parameters to improve the initial model.
"""

# Run guess fit with given parameters (see docs for run_fit_guess)
if runSFit:
    comm.run_fit_guess(sfit, maxiter=maxguess, ar=areg, outname='')


#%% Fit strain - run regularization loop

if runSFit and runReg:
    comm.run_fit_reg(sfit, maxiter=maxiter, ar=aregs, outname='')


#%% Fit strain - run final fit

"""
Choose the optimum regularization value (areg) and run the final fit.
"""
areg = 3
if runSFit:
    comm.run_fit(sfit, maxiter=maxiter, ar=areg, outname=scan['epsfile'], 
                 bootstrap=bootstrap, loops=loops)





