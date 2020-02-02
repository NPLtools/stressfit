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
dhkl=Fe_110;

%% Instrument setup
% angles
angle_min=75;
angle_max=105;
theta=45;
% sample thickness
d_sam=10;
% collimation
s0w=14;
s0d=4000;
s1w=2;
s1d=60;
s2w=2;
% extinction
extObj=Extinction([0, 0.034, 0.121, 7.82e-2]);

% depth scan range [mm] (negative = outside the sample)
zmin=-2;
zmax=d_sam+2;
nz=71;

%% Initialize
INS=ENGINX([inppath,'spectrum.dat'],[inppath,'pulse.dat'],'ENGINX.par',extObj);
na=11;
da=(angle_max-angle_min)/(na-1);
arange=angle_min:da:angle_max;

fprintf('wavelength = %f\n',2*dhkl*sin(theta/2*deg));
fprintf('wavelength range = %f, %f\n',2*dhkl*sin(angle_min/2*deg),2*dhkl*sin(angle_max/2*deg));

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
gwidth=zeros(nz,2);
pshift=zeros(nz,2);

% transmission
omega=90+theta;
INS.init(arange, omega, dhkl, d_sam, s0w,s0d,s1w,s1d,s2w);
[DSE, beta, zeta, res] = INS.getPeakParam(z);
fprintf('Transmission:\n');
fprintf('\tDSE = %f\n',DSE);
fprintf('\tbeta = %f [1/mm]\n',beta);
fprintf('\tzeta = %f\n',zeta);
fprintf('\tgauge fwhm [mm] = %f\n',beta/g_norm);
zscale(1:nz,1)=res(:,2);
pshift(1:nz,1)=res(:,3);
gwidth(1:nz,1)=res(:,4);
table=[DSE, beta, zeta];
% reflection
omega=theta;
INS.init(arange, omega, dhkl, d_sam, s0w,s0d,s1w,s1d,s2w);
[DSE, beta, zeta, res] = INS.getPeakParam(z);
fprintf('Reflection:\n');
fprintf('\tDSE = %f\n',DSE);
fprintf('\tbeta = %f [1/mm]\n',beta);
fprintf('\tzeta = %f\n',zeta);
fprintf('\tgauge fwhm [mm] = %f\n',beta/g_norm);
zscale(1:nz,2)=res(:,2);
pshift(1:nz,2)=res(:,3);
gwidth(1:nz,2)=res(:,4);
table=cat(1,table,[DSE, beta, zeta]);
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
fwhm = INS.getResolution(drange); 
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
arange=-90:2:+90;

str1=sprintf('Orientation dependence');
grange=theta+90+arange;
INS.setTheta(theta);
[DSE, beta, zeta, p] = INS.scan_gamma(grange);
INS.setTheta(theta-15.0/2);
[DSE1, beta1, zeta1, p1] = INS.scan_gamma(grange);
INS.setTheta(theta+15.0/2);
[DSE2, beta2, zeta2, p2] = INS.scan_gamma(grange);
myfigure(arange',cat(2,DSE,DSE1,DSE2),str1,opt,leg);
dlmwrite([outpath,prefix,'gamma.dat'],cat(2,grange',DSE, DSE1, DSE2, beta, zeta, p),' ');
INS.setTheta(theta);

%%
% <html>
% <b>dhkl dependence</b><br\>
% </html>
% 
drange=0.5:0.1:2.4;
INS.setGamma(theta+90);
[DSER, betaR, zetaR] = INS.scan_dhkl(drange);
INS.setGamma(theta+180);
[DSET, betaT, zetaT] = INS.scan_dhkl(drange);

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


