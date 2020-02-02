%% Analytical model for neutron diffraction peak shifts due to the surface effect
% written by: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
% version: Time-of-Flight
% last update: Sept 2016

clear all;
close all;
format shortG;
inppath='./input/';
outpath='./output/';
addpath('libmat-tof','lib','libtools');
% dhkl values
Fe_110=2.0269;
Fe_200=1.4332;
Fe_211=1.1702;
Fe_220=1.0135;
Fe_310=0.9075;
Fe_321=0.7654;

%% Distribution shape factors and other constants

% uniform distribution                 
g_uni=sqrt(1/6);
% circular distribution                 
g_circ=sqrt(1/8);
% triangular distribution
g_tri=sqrt(1/3);
% gaussian distribution
g_norm=1/sqrt(4*log(2));
% degree
deg=pi/180;
% Set excluded parameters to infinity
INF=1e10;

%% USER INPUT

% prefix for output files
prefix='enx_S14_Fe211_';

% reflection
d_sam=Fe_211;

% Instrument setup

% SOURCE
pulse_width=88.4;  % pulse width [us] 
pulse_shape=g_norm;   % pulse shape factor
src_dist=49730;  % distance [mm] from the sample
% get tables with wavelength dependence of flux and pulse width
% 2-column format with wavelength [A], flux or fwhm [us] (gaussian width)
pulse_table=readXY([inppath,'pulse.dat']);
flux_table=readXY([inppath,'spectrum.dat']);


% NEUTRON GUIDE
% included only if s0_on=0 (divergence slit is not used)
gd_width=25;     % width [mm]
gd_dist=1500;    % distance [mm] of the guide exit from the sample
gd_m=3;          % coating m-value

% DIVERGENCE SLIT
s0_width=14;    % width [mm]
s0_dist=4000;   % distance [mm] from the sample
%s0_width=12;    % width [mm]
%s0_dist=1500;   % distance [mm] from the sample
s0_on=1;        % = 1 to include this slit instead of the neutron guide

% PRIMARY SLIT
s1_width=2;    % width [mm]
s1_dist=60;    % distance [mm] from the sample

% SAMPLE
extObj=Extinction([0 0.034, 0.121, 7.82e-2]); % extinction function (see Extinction class for help)
% mu=[0.034, 0.121, 7.82e-2]; % attenuation coeff. mu = mu[1] + mu[2]*lambda +  mu[2]*lambda^2 in [cm^-1]
thick=10;       % sample thickness [mm]

% RADIAL COLLIMATOR or SECONDARY SLIT
s2_width=2.0;  % gauge width [mm]
s2_dist=0;      % set to zero for a radial collimator

% DETECTOR
d_binwidth=3;     % resolution fwhm [mm]
d_binshape=g_uni; % uni for uniform, norm for gaussian bin
d_dist=1530;      % distance from the sample
d_angle_min=75*deg;  % minimum angle (integration range)
d_angle_max=105*deg; % maximum angle (integration range)
d_resol=1*g_norm;        % time resolution [us]

% WAVELENGTH

% depth scan range [mm] (negative = outside the sample)
zmin=-2;
zmax=thick+2;
nz=71;

%% Notes
%
% To consider radial collimators instead of slits, set corresponding
% distances to zero (e.g. s2_dist=0 for secondary radial collimator. The width
% is then FWHM of the spatial distribution of neutron flux transmitted by
% the collimator at the focal plane (assumed to be triangular).
% collimator.
% 
% <html>
% <b>NOTE</b><br/>
% transmission widths of individual components should be weighted by 
% factors relating <i>fwhm</i> to the <i>variance</i> of the corresponding
% distribution. This is done in the objects implementing the caclulation
% method, using constants g_norm, g_tri, g_step.
% </html>
%   
% 

%% Instrument setup
% Note the indexing and signs of distances: 
% Before the sample, tracing is assumed in up-stream
% direction, therefore the distances are negative and refer to the next
% component along the neutron stream.

% mean Bragg angle
theta=0.25*(d_angle_min+d_angle_max);
% mean wavelength
lam0=2*d_sam*sin(theta);

i=0;
% SOURCE
% parameters: distance [mm], wavelength [A], tau [us], shape (see NOTE)
if (s0_on) 
srcL=-src_dist+s0_dist;
else
srcL=-src_dist+gd_dist;  
end
i=i+1;src(i)=Pulse(srcL,pulse_width, pulse_shape);
%
i=i+1;
% DIVERGENCE SLIT
% parameters: distance [mm], wavelength [A], width[mm]
if (s0_on) 
i=i+1;src(i)=Slit(-s0_dist+s1_dist,s0_width);
else
i=i+1;src(i)=Guide(-gd_dist+s1_dist,gd_width,gd_m);    
end
% INPUT SLIT or RADIAL COLLIMATOR
% parameters:  distance [mm], wavelength [A], width[mm]
i=i+1;  
src(i)=Slit(-s1_dist,s1_width);

% SAMPLE
gamma_R=MX_model_tof.g_ref(theta);     % surface angle [rad]
gamma_T=MX_model_tof.g_tra(theta);     % surface angle [rad]
% parameters: Bragg angle [rad], surface angle [rad], attenuation coeff. [mm^-1]
i=i+1;  src(i)=Sample(d_sam,gamma_R,thick);
i_sam=i;
fprintf('Sample reflection, 2 thetaS=%f, gammaR=%f, gammaT=%f\n',2*theta/deg, gamma_R/deg, gamma_T/deg);

% OUTPUT SLIT or RADIAL COLLIMATOR
% parameters: distance [mm], wavelength [A], width[mm]
i=i+1;  
src(i)=Slit(s2_dist, s2_width);

% DETECTOR
i=i+1;  src(i)=TofDetector(d_dist-s2_dist, d_binwidth, d_binshape,d_resol,d_angle_min, d_angle_max);
i_det=i;

% Define component sequence for the special case 
% of a thin bent monochromator and open in-pile collimation
% requires: monochromator, slit 0, slit 1, slit 2 and detector
% ATTENTION: depends on the order and type of components in the src array !
%spec=[src(i_mon),src(i_mon+1),src(i_mon+2),src(i_sam+1),src(i_sam+2)];

fprintf('wavelength = %f\n',lam0);
fprintf('wavelength range = %f, %f\n',2*d_sam*sin(d_angle_min/2),2*d_sam*sin(d_angle_max/2));

% Create the object which implements the matrix method
G=MX_model_tof(flux_table,pulse_table,extObj);
% set new instrument sequence and calculate T and C matrices
G.setParam(src,theta,i_sam,i_det);  

%%
% <html>
% <b>depth-scale and peak shift</b><br\>
% </html>
% 
clear z;
clear zscale;
dz=(zmax-zmin)/(nz-1);
z=zmin:dz:zmax;
zscale=zeros(nz,2);
na=11;
da=(d_angle_max-d_angle_min)/(na-1);
trange=0.5*(d_angle_min:da:d_angle_max);

gwidth=zeros(nz,2);
pshift=zeros(nz,2);

PS = PSModel();
PS.setParam(thick,0.0,0.0);



% transmission
src(i_sam).gamma=gamma_T;
G.setParam(src,theta,i_sam,i_det);
[DSE, beta, zeta, p] = G.scan_theta(src,trange); 
mDSE=PSModel.meanp(DSE,p);
mbeta=PSModel.meanp(beta,p);
mzeta=PSModel.meanp(zeta,p);
fprintf('Transmission:\n');
fprintf('\tDSE = %f\n',mDSE);
fprintf('\tbeta = %f [1/mm]\n',mbeta);
fprintf('\tzeta = %f\n',mzeta);
fprintf('\tgauge fwhm [mm] = %f\n',mbeta/g_norm);
zscale(1:nz,1)=PS.get_ctr_ave(z',beta,zeta, p);
pshift(1:nz,1)=PS.get_shift_ave(z',DSE,beta,zeta, p);
gwidth(1:nz,1)=PS.get_gwidth_ave(z',beta,zeta, p);
table=[mDSE, mbeta, mzeta];
% reflection
src(i_sam).gamma=gamma_R;
G.setParam(src,theta,i_sam,i_det);
[DSE, beta, zeta, p] = G.scan_theta(src,trange);
mDSE=PSModel.meanp(DSE,p);
mbeta=PSModel.meanp(beta,p);
mzeta=PSModel.meanp(zeta,p);
fprintf('Reflection:\n');
fprintf('\tDSE = %f\n',mDSE);
fprintf('\tbeta = %f [1/mm]\n',mbeta);
fprintf('\tzeta = %f\n',mzeta);
fprintf('\tgauge fwhm [mm] = %f\n',mbeta/g_norm);
zscale(1:nz,2)=PS.get_ctr_ave(z',beta,zeta, p);
pshift(1:nz,2)=PS.get_shift_ave(z',DSE,beta,zeta, p);
gwidth(1:nz,2)=PS.get_gwidth_ave(z',beta,zeta, p);
table=cat(1,table,[mDSE, mbeta, mzeta]);
dlmwrite([outpath,prefix,'table.dat'],table,' ');

% plot options
opt=Utils.plotOptions();
opt.xauto=false;
opt.xlimits=[zmin zmax];
opt.yauto=false;
opt.LineWidth=2;
opt.groups=[1 2];
opt.types=[0 0];
leg={'transmission','reflection'};

opt.xlabel='surface position, mm';
opt.ylabel='gauge shift, mm';
opt.LegendPosition='SouthEast';
str1=sprintf('Depth scale conversion');
ymin=min(0,min(min(zscale)));
ymax=max(max(zscale));
ymax=ymax +0.1*(ymax-ymin);
opt.ylimits=[ymin ymax];
myfigure(z,zscale,str1,opt,leg);
dlmwrite([outpath,prefix,'zscale.dat'],cat(2,z',zscale),' ');

opt.LegendPosition='South';
opt.ylabel='fwhm, mm';
str1=sprintf('Gauge volume width');
ymin=min(0,min(min(gwidth)));
ymax=max(max(gwidth));
ymax=ymax +0.1*(ymax-ymin);
opt.ylimits=[ymin ymax];
myfigure(z,gwidth,str1,opt,leg);
dlmwrite([outpath,prefix,'gwidth.dat'],cat(2,z',gwidth),' ');

opt.LegendPosition='South';
opt.ylabel=sprintf('Spurious strain, \\mue ');
str1=sprintf('Spurious strain');
ymin=min(0,min(min(pshift)));
ymax=max(max(pshift));
ymax=ymax +0.1*(ymax-ymin);
opt.ylimits=[ymin ymax];
myfigure(z,pshift,str1,opt,leg);
dlmwrite([outpath,prefix,'strain.dat'],cat(2,z',pshift),' ');

%%
% <html>
% <b>Resolution</b><br\>
% </html>
% 
drange=0.5:0.1:2.4;
fwhm = G.scan_dhkl_fwhm(src,drange); 
opt.showlegend=0;
opt.xlabel=sprintf('dhkl, A');
opt.ylabel=sprintf('\\Deltad/d, %%');
str1=sprintf('Resolution');
ymin=min(0,min(min(fwhm)));
ymax=max(max(fwhm));
ymax=ymax +0.1*(ymax-ymin);
opt.ylimits=[ymin ymax];
opt.xlimits=[min(drange) max(drange)];
myfigure(drange',fwhm,str1,opt,{'resolution'});
dlmwrite([outpath,prefix,'fwhm.dat'],cat(2,drange',fwhm),' ');


%%
% <html>
% <b>Orientation dependence</b><br\>
% </html>
% 
opt.xauto=0;
opt.xlimits=[-90, 90];
opt.yauto=1;
opt.showlegend=1;
opt.groups=[1 2 2];
opt.types=[0 0 0];
colors=[[0 0 0]' [0 0 0]']';
leg={'2\theta=90 deg','2\theta=90 +- 15 deg'};
opt.xlabel=sprintf('\\chi, deg');
opt.ylabel=sprintf('\\Delta_{SE}, \\mu\\epsilon/mm');
opt.LegendPosition='NorthEast';
arange=-90*deg:2*deg:+90*deg;

str1=sprintf('Orientation dependence');
grange=gamma_R+arange;
G.setParam(src,theta,i_sam,i_det);
[DSE, beta, zeta, p] = G.scan_gamma(src,grange);
G.setParam(src,d_angle_min/2,i_sam,i_det);
[DSE1, beta1, zeta1, p1] = G.scan_gamma(src,grange);
G.setParam(src,d_angle_max/2,i_sam,i_det);
[DSE2, beta2, zeta2, p2] = G.scan_gamma(src,grange);
myfigure(arange'/deg,cat(2,DSE,DSE1,DSE2),str1,opt,leg);
dlmwrite([outpath,prefix,'gamma.dat'],cat(2,grange'/deg,DSE, DSE1, DSE2, beta, zeta, p),' ');

%%
% <html>
% <b>dhkl dependence</b><br\>
% </html>
% 
drange=0.5:0.1:2.4;
src(i_sam).gamma=gamma_R;
[DSER, betaR, zetaR] = G.scan_dhkl(src,drange);
src(i_sam).gamma=gamma_T;
[DSET, betaT, zetaT] = G.scan_dhkl(src,drange);

str1=sprintf('dhkl dependence');
opt.groups=[ 1 2 ];
opt.types=[ 2 2 ];
opt.showlegend=1;
opt.xlimits=[min(drange),max(drange)];
leg={'transmission','reflection'};
opt.xlabel=sprintf('\\chi, deg');
opt.ylabel=sprintf('\\Delta_{SE}, \\mu\\epsilon/mm');
opt.LegendPosition='East';
opt.xlabel=sprintf('dhkl, A');
opt.ylabel=sprintf('\\Delta_{SE}, \\mu\\epsilon/mm');
opt.yauto=false;
ymin=min(0,min(min(DSET,DSER)));
ymax=max(max(DSET,DSER));
ymax=ymax +0.1*(ymax-ymin);
ymin=ymin -0.1*(ymax-ymin);
opt.ylimits=[ymin ymax];

myfigure(drange',cat(2,DSET,DSER),str1,opt,leg);
dlmwrite([outpath,prefix,'dhkl.dat'],cat(2,drange',DSET, DSER, betaT, betaR, zetaT,zetaR),' ');


return;


