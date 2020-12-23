# STRESSFIT - fitting of residual stress distributions
*Written by*: Jan Saroun, Nuclear Physics Institute CAS, Rez, saroun@ujf.cas.cz  
*Repository*: https://github.com/NPLtools/stressfit

-------------------

STRESSFIT is a Python package for fitting of neutron residual strain scanning data. The primary objective is to restore intrinsic lattice strain distributions measured by neutron diffraction. In an experiment, this distribution is smeared by finite spatial resolution of the instrument. Measured strains are also affected by so called pseudo-strains, which arise in general from inhomogeneous sampling of the material. The best known is the *surface effect*, when the nominal gauge volume of the instrument is only partly immersed under the sample surface. this gives rise to peak shifts (pseudo strains) of similar magnitude as the intrinsic strain. More generally, similar effects occur when the local scattering probability varies significantly over the scale comparable with the gauge volume. Possible causes are gradients in texture, phase composition or abrupt change in beam attenuation near a Bragg edge.

In STRESSFIT, this problem is solved by indirect deconvolution of measured strain profiles: a parametric model of the intrinsic strain distribution is convoluted with the sampling distribution. The free parameters of the model are then fitted to the experimental data - measured strains. The sampling distribution is obtained by neutron ray-tracing simulation of the instrument in given configuration, which export a list of scattering points associated with relevant scattering event information: position, incident and scattered wave vectors, event probability and strain value detected by the instrument for this event. Such a sampling distribution can be simulated rather quickly for any instrument setup, for example by the ray-tracing program SIMRES (http://neutron.ujf.cas.cz/restrax). The convolution is done in 3D by summing up probabilities of the sampling events multiplied by the local property of the material: scattering probability or strain, yielding corresponding folded quantity at given sample position. The sampling distribution is then reused to carry out the convolution at any sample position, with any sample shape, orientation or strain distribution. This decoupling of MC ray-tracing simulation from the convolution procedure permits to write a code which is rather fast (typically about 1s per one convolution of a scan  at 1 CPU). It is therefore possible to use the MC integration as a part of cost function for least squares fitting.
</p><p>
Currently, STRESSFIT enables to model pseudo strains for several basic sample shapes (curved plates, full and hollow cylinders and spheres). It also enables least squares fitting of measured scans for a single strain component. Simultaneous fitting of multiple components with a single stress distribution model is envisaged in future versions.

</p>

## Documentation

The package contains *Jupyter notebook files* with comments that guide user through the fitting process. An example with output of STRESSFIT is also available via Jupyter viewer server:
<p>
<a href='http://nbviewer.jupyter.org/url/neutron.ujf.cas.cz/restrax/download/stressfit/stressfit_example1.ipynb'>
Example 1</a>: ECNS2019, Fitting of strain gradient under the inner surface of a tube, test on synthetic data for STRESS-SPEC.
</p>
<p>
For more information and use examples, see: <br/>
<a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/saroun_MECASENS_2017.pdf'>MECASENS 2017 poster</a><br/>
<a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/ECRS2018_stressfit.pdf'>ECRS10, 2018, slides</a><br/>
<a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/saroun_ECNS2019_poster.pdf'>ECNS 2019, poster</a> <br/>
</p>
<p>
For access to the latest source code, contact the author.
</p>


## Package contents

### Execution scripts:

`templates/stressfit-template.py`:	A commented template  
available also as Jupyter notebook  files, *.ipynb (see http://jupyter.org)  

### Module files:
The STRESSFIT package files are in the subdirectory `stressfit`: 
`sample.py`: handling of sample properties and MC convolution  
`mccfit.py`: classes and functions for data fitting using MC convolution  
`graphs.py`: some plot functions and definitions  
`shapes/*.py`: classes with sample shape definitions and methods  

### Example data:
`input/events_B_1mm.dat`: sampling events for BEER@ESS, Fe(211), gauge volume ~ 1x1x3 mm  
`input/events_S_1mm.dat`: sampling events for STRESS-SPEC, Fe(211), gauge volume  ~ 1x1x3 mm, lambda=1.68 A  
`tables/Fe_mu.dat`: Total removal cross-section table for alpha-Fe
`input/*.dat`: example input data (synthetic data simulated by SIMRES)

### Matlab tools:
the `matlab` folder contains older matlab scripts implementing the analytical method for calculation of pseudo-strains based on instrument setup and diffraction geometry parameters. It employs the matrix description of neutron transport through the instrument and derives analytical formulas for the pseudo strain in Gaussian approximation. See https://doi.org/10.1107/S0021889813008194 for more details. 

### Other:

`tools/fit_centre.py`:	example script for fitting of sample centre from two perpendicular intensity scans.  
`tools/eps2sig.py`: script for conversion of strain distributions to stresses  
`input/strain_table.dat`:   strain and intensity distribution table used to simulate synthetic data in SIMRES
