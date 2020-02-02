%% Surface effect fit
% Fit the peak shift due to the surface effect measured as a function of 
% depth under the surface.
%
% see also J. Saroun et al., J. Appl. Cryst. (2013) 46, 628-638.
%
% version: 1.21
% date: 22/6/2016
%
close all;
clear all;

%% Definitions
addpath('lib');
inpath='./input/';
outpath='./output/';
deg=pi/180;

%% Options and settings
% Bragg angle
theta_B=90.0/2*deg;
% estimate of the attenuation coefficient for symmetric reflection geometry
mu=1.23; % [1/cm]

%% External input data
% Provide 3-column data with: encoder position [mm], shift [microstrain=10^-6], error
inpfile='shift_sim_S14R.dat';
% transmission geometry?
trans_geom=0;
% deviation from the symmetric geometry (transmission or reflection)
gamma=0.0*deg;
%
logfile=strcat(inpfile,'.log');

%% Initial values. 
% a and z0 should usually be determined independently and fixed.
% beta can be determined from intensity scan or calculated analytically.
%
% offset in microstrain (correction for d0 or aberation effects etc.)
y0=0;           
% amplitude of the peak shift function [microstrain/mm]
%DSE=198;  % tra 
DSE=433.4; % ref
% gauge volume width parameter [1/mm]
beta=0.982; % tra
%beta=1.021; % ref
% zero position [mm]
z0=-5; % tra
%z0=-0.024; % ref

% sample thickness [mm]
d=10.0;

%% Fit data
% read data
shift=readXYZ(strcat(inpath,inpfile));

% attenuation parameter [1/mm]
g=trans_geom*pi/2;
a0=mu*0.1*abs(1/sin(theta_B+gamma+g)+1/sin(theta_B-gamma-g));    
if (trans_geom); a=0; else a=a0; end;  

% assign initial parameter values and set fixed attributes
% Usually, beta,a,z0 and d are determined independently and can be fixed.
param=[y0 DSE beta a z0 d];
fixed=[0 0 0 1 0 1];
% fit
[ output, err, shiftFit ] = fit_peakshift( shift, param,fixed,1 );
resultFit=cat(2,shift(:,2),shiftFit(:,1));

% assign output parameters back to local variables
pCell = num2cell(output(1:6));
[y0 DSE beta a z0 d] = pCell{:};
pCell = num2cell(err(1:6));
[erry0 errDSE errbeta erra errz0 errd] = pCell{:};

% Start log file
flog=strcat(outpath,logfile);
ReportInit(flog,mfilename('fullpath'));
% Log results
names={'offset' 'D_SE' 'beta' 'a' 'z0' 'thickness'};
header='Fit of the surface effect';
ReportFitResults(flog,1,header,names,output,err,fixed);

%% Plot result 
% conversion between surface position and true depth of the gauge volume
% center
opt.xauto=true;
opt.yauto=true;
opt=Utils.plotOptions();
opt.xlimits=[-20 20];
opt.ylimits=[0 12000];
opt.xlabel='surface position, mm';
opt.ylabel='peak shift [microstrain]';
opt.LegendPosition='South';
opt.LineWidth=2;
opt.XMinorTick='on';
opt.YMinorTick='on';
opt.groups=[1 2];
opt.types=[1 0];
leg={'data', 'fit'};


fprintf('Attenuation coefficient for symmetric reflection geometry: \n');
fprintf('\ta=%f \n',a0);

myfigure(shift(:,1),resultFit,'Fit of the surface effect',opt,leg);
dlmwrite([outpath, 'fit_',inpfile],cat(2,shift(:,1),resultFit),' ');

