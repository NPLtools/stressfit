# STRESSFIT - fitting of residual stress distributions
*Written by*: Jan Saroun, Nuclear Physics Institute CAS, Rez, saroun@ujf.cas.cz  
*Repository*: https://github.com/NPLtools/stressfit

-------------------
STRESSFIT is a Python package for evaluation of lattice strain distributions measured by neutron diffraction. It addresses the problem of large neutron sampling volume, which causes a number of undesired effects commonly called pseudo-strains:
- Smearing of measured strain distributions
- False strains observed due to non-uniform sampling, caused by
    - scanning through sample surface or other material boundaries
    - steep gradients in composition and texture
    - variation of beam attenuation during a scan (sample shape effect, Bragg edge). 

On the input, STRESFIT uses a list of neutron scattering events with associated data (position, wave vectors, weight factor and "as-measured" lattice spacing - $d_{hkl}$}) accessible at given instrument setup. This list describes the instrumental sampling distribution independently on the sample. It can be obtained by ray-tracing simulation of the instrument using appropriate software, such as McStas (http://mcstas.org) or SIMRES (https://github.com/saroun/simres). 

STRESSFIT provides tools for 3D convolution of such a sampling list with the sample model and permits to calculate: 

- “centre of gravity” and size of the neutron sampling volume as a function of sample position (in 3D),
- variation of intensity and position of diffraction peaks due to the perturbation of sampling distribution (material boundaries, absorption, composition and texture gradients),
- “as measured” (smeared) intensity and strain distributions including the pseudo-strain effects,
- least-squares fit of intrinsic strain and intensity distributions.

</p><p>
STRESSFIT enables to model pseudo-strains for several sample shapes such as curved plates, cylinders, spheres, tubes (both single- and multi-channel) and polygonal rods. Neutron attenuation tables for several common materials generated with the help of the NCrystal library (https://github.com/mctools/ncrystal) are provided as a part of the package resources. Currently, STRESSFIT enables least squares fitting of measured scans for a single strain component. Simultaneous fitting of multiple strain or stress tensor components is envisaged in future versions.
</p>

-----------------------------

## Documentation

The package contains commented *Jupyter notebook* template which guides users through usual workflow of the fitting process. An example with output of STRESSFIT is also available via Jupyter viewer server:
<p>
<a href='http://nbviewer.jupyter.org/url/neutron.ujf.cas.cz/restrax/download/stressfit/stressfit_example1.ipynb'>
Example 1</a>: ECNS2019, Fitting of strain gradient under the inner surface of a tube, test on synthetic data for STRESS-SPEC.
</p>
<p>
For more information and use examples, see: <br/>
<a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/ECRS2018_stressfit.pdf'>ECRS10, 2018, slides</a><br/>
<a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/saroun_ECNS2019_poster.pdf'>ECNS 2019, poster</a> <br/>
<a href='http://neutron.ujf.cas.cz/restrax/download/stressfit/stressfit_ECNS2023_poster.pdf'>ECNS 2023, poster</a> <br/>
</p>

-----------------------------

## Installation

Get the package source either by cloning from git repository

<code>git clone https://github.com/NPLtools/stressfit</code>

or download it directly from github at the link above.

Unpack it and, from the root package directory, run

<code>pip install -e .</code>

You can of course use any other alternative to install stressfit into a suitable python environment. It is recommended to use Anaconda distribution, which should provide most of the required packages by default. Yet you may need to add some other packages manually like lmfit.

-----------------------------

## Graphical user interface

An experimental version of GUI is available. It runs in Jupyter notebook using [ipywidgets](https://ipywidgets.readthedocs.io).

To start it, execite the following code in the Jupyter notebook in your browser:
<code>  
import stressfit.ui.notebook as nb  
ui = nb.UI()  
ui.display()  
</code>  

It is assumed that the Jupyter server runs on localhost (running on remote server is not supported) and stressfit has to be installed in the same environment. The GUI implements all main features of the STRESSFIT package. On startup, it creates default workspace in the user's profile and loads test data from resources. It is thus ready to run all commands, such as plotting the 2D view o experiment geometry with the sampling distribution, evaluate the spatial resolutiona and pseudo-strain characteristics, compare the test data with pseudo-strain distributions and also fit intensity and strain distributions. The GUI also allows to save/load complete setup data for later use.

-----------------------------

## Package contents

### Execution scripts:

`stressfit/resources/data`: Test input data and simulated neutron sampling lists.
`stressfit/resources/tables`: Atttenuation tables or some common materials
`stressfit/resources/conf`: Default configuration data used by the GUI
`templates/stressfit-template.py`:	A commented template  
available also as Jupyter notebook  files, *.ipynb.  

### Example data:
In `stressfit/resources/data`:

- `ievents_B_1mm.dat`: sampling events for BEER@ESS, Fe(211), gauge volume ~ 1x1x3 mm  
- `events_S_1mm.dat`: sampling events for STRESS-SPEC, Fe(211), gauge volume  ~ 1x1x3 mm, lambda=1.68 Ang
- `*.dat`: example input data (synthetic data simulated by SIMRES)

In `stressfit/resources/tables`:
`tables/Fe_mu.dat`: Neutron atttenuation tables or some common materials 

In `stressfit/resources/conf`:
`*.json`: Default configuration data used by the GUI (required when starting GUI, but also usable as templates for user projects)

### Matlab tools:
the `matlab` folder contains older matlab scripts implementing the analytical method for calculation of pseudo-strains based on instrument setup and diffraction geometry parameters. It employs the matrix description of neutron transport through the instrument and derives analytical formulas for the pseudo strain in Gaussian approximation. See https://doi.org/10.1107/S0021889813008194 for more details. 


