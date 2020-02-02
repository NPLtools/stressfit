%% Intensity scan fit
% Fit the intensity measured as a function of depth under the surface.
%
% See also J. Saroun et al., J. Appl. Cryst. (2013) 46, 628-638.
%
% updated for layers (two parallel surfaces)
%
% version: 1.31
% date: 22/6/2016
%
close all;
clear all;

%% Definitions
addpath('lib', 'libtools');
inpath='./input/';
outpath='./output/';
deg=pi/180;


%% User Input
%---------------------------------------------------------
% prefix for input/output files
prefix='int_211_S14';

logfile=cat(2,prefix,'.log');
% Bragg angle [deg]
theta_S=90.0/2*deg;

% The script can fit multiple data sets in a single run
% Provide 3-column data: encoder position [mm], intensity, error
inpfiles={cat(2,prefix,'R.dat'), cat(2,prefix,'T.dat')};

% fix parameters b,c
fxbc=0;

% For each input file, provide additional parameters: 
% Bragg angle
theta=[theta_S, theta_S];
% geometry: reflection (0) or transmission (1)
trans_geom=[0, 1];
% deviation from symmetric condition
gamma=[0, 0]*deg; 

% provide initial parameters for each data file (one column per input file).
% The order of rows is:
% param=[y0, I0, z0, beta, mu, d, b, c]'
% y0   ... background
% I0   ... intensity
% z0   ... surface position [mm]
% beta ... gauge volume size parameter [mm-1]
% mu   ... Attenuation coefficient [1/cm], estimate
% d    ... sample thickness [mm] 
% b    ... linear inhomogeneity coefficient
% c    ... quadratic inhomogeneity coefficient -0.0425, 0.00568

% mu values for 110, 200, 211, 220, 310, 311
m=[1.30, 1.25, 1.147, 1.156, 1.093, 1.019 ];
% parameters
param=[1, 3000, -5.0, 1.0, 1.147, 10, 0, 0;...
       1, 1000, -5.0, 1.0, 1.147, 10, 0, 0]';

% 110, S25
%param=[1, 3000, -5.0, 1.0, 1.30, 10, -0.0425, 0.00568;...
%       1, 1000, -5.0, 1.0, 1.30, 10, 0.24481, -0.0238]';

% 200, S25
%param=[1, 3000, -5.0, 1.0, 1.25, 10, 0.121595, -0.011620;...
%       1, 1000, -5.0, 1.0, 1.25, 10, -0.065362, 0.007320]';

% 211, S25
%param=[1, 3000, -5.0, 1.0, 1.147, 10, 0.026856, -0.003454;...
%       1, 1000, -5.0, 1.0, 1.147, 10, -0.004307, 0.002216]';

% provide flags indicating which parameters are fixed (0|1)
fixed=[1, 0, 0, 0, 1, 1, fxbc, fxbc;...
       1, 0, 0, 0, 1, 1, fxbc, fxbc]';
%---------------------------------------------------------

%% Calculate dependent parameters
% ag is the geometry factor for depth-dependent transmission term
% cg is a geometry factor for the constant transmission term (set 0 for
% reflection geometry)
nfiles=numel(inpfiles);
ag=zeros(nfiles);
cg=zeros(nfiles);
for i=1:nfiles
    if (trans_geom(i)>0)
        ag(i)=0.1*(1/cos(theta(i)+gamma(i))-1/cos(theta(i)-gamma(i)));
        cg(i)=0.1/cos(theta(i)-gamma(i));
    else
        ag(i)=0.1*(1/sin(theta(i)+gamma(i))+1/sin(theta(i)-gamma(i))); 
        cg(i)=0.0;
    end;    
    if (fixed(5,i)==0)
       param(5,i)=max(param(5,i),0.01);
    end;        
end;


%% Fit the intensity scan
% The dependence of scattered intensity on the surface position with
% respect to the gauge volume center can be expressed analytically,
% assuming the Gaussian model from the JAC paper:
%
% $$ I(z) = I_{bcg} + A F_0(\beta(z-z_0), \zeta, d)  $$
%
% where $d=\beta thickness$ and
%
% $$F0(h, \zeta, d) = \frac{\sqrt{\pi}}{2}  e^{-\zeta(2h-\zeta)} [erf(h-\zeta) - erf(h -\zeta - d)]$$
%
% We can fit this function to the measured depth dependence of scattered
% intensity to obtain:
% $z_0$ (surface position), $\beta^{-1}$ (the effective size of the gauge volume) 
% and $\zeta$ (absorption parameter). It is also possible to fit the layer
% thickness _t_, but it should better be fixed if the value is already
% known.
%
% $\zeta$ is proportional to the absorption coefficient. For the symmetric 
% reflection geometry, the absorption coefficient [1/mm] is:
%
% $$ \mu = \beta \zeta \sin(\theta_B) $$
% 
% It is also related to the absorption coefficient _a_ used by Hutchings:
%
% $$\zeta = \frac{a}{2 \beta}$$
%

% initialize fitting model
M=IntensityModel();

% Start log file
flog=strcat(outpath,logfile);
ReportInit(flog,mfilename('fullpath'));

% Plot options - fit
opt=Utils.plotOptions();
opt.xauto=true;
opt.yauto=true;
opt.xlimits=[-20 20];
opt.ylimits=[0 12000];
opt.xlabel='nominal depth, d [mm]';
opt.ylabel='Intensity';
opt.LegendPosition='South';
opt.LineWidth=2;
opt.XMinorTick='on';
opt.YMinorTick='on';
opt.groups=[1 1 2];
opt.types=[1 0 0];
opt.styles={'-' '--'};

% Plot options - z-scale
opt1=Utils.plotOptions();
opt1.xauto=true;
opt1.yauto=true;
opt1.xlimits=[-3 13];
opt1.ylimits=[-3 13];
opt1.xlabel='nominal depth, d [mm]';
opt1.ylabel='true depth [mm]';
opt1.LegendPosition='NorthWest';
opt1.LineWidth=2;
opt1.XMinorTick='on';
opt1.YMinorTick='on';
opt1.groups=[1 2];
opt1.types=[0 0];
opt1.styles={'-' '--'};



%% Read and process input data
% files with 3-columns: encoder position [mm], intensity, error
nfiles=numel(inpfiles);

for i=1:nfiles
fname= inpfiles{i}; 
S=strrep(fname,'_','\_');
hdr=sprintf('Fitting file: %s',fname);
% load data
data=readXYZ(strcat(inpath,fname));
nz=size(data,1);
% set Bragg angle and geometry factors for absorption
M.theta_B=theta(i);
M.ag=ag(i);
M.cg=cg(i);

% run fit
M.setParam1(param(:,i),fixed(:,i));
M.setX(nz,data(:,1));
[ output,errors,fit,chi2,Cov ] = fit_model(data,M,1);
[y0, I0, z0, beta, mu, t, b, c] = M.getResult(1);

% print report
M.reportResults(flog,1, hdr,output,errors,chi2);

% calculate fit without inhomogeneity correction
if (b~=0 || c~=0)
  param1=output';
  param1(7:8)=0.0;
  M.setParam1(param1,fixed(:,i));
  [ output,errors,fith,chi2,Cov ] = fit_model(data,M,0);
  % fith = M.fitFnc();  
end

% plot fit
leg={'data','fit','fit without inhomogeneity'};
if (b~=0 || c~=0)
  fdata=cat(1,data(:,1:2),cat(2,data(:,1),fit(:,1)),cat(2,data(:,1),fith(:,1)));
  result=cat(2,data(:,1:2),fit(:,1),fith(:,1));
  ndat=[nz nz nz];
else
  fdata=cat(1,data(:,1:2),cat(2,data(:,1),fit(:,1)));
  result=cat(2,data(:,1:2),fit(:,1));
  ndat=[nz nz];
end
plotData(ndat,fdata,S,leg,'model',opt);
% plotFit(nz,data(:,1),data(:,2),nz,data(:,1),fit,S,leg,'measurement','fit',opt);
dlmwrite([outpath,'fit_',fname],result,' ');

% save fit on equidistant scale
zmin=data(1,1);
zmax=data(nz,1);
z=zmin:(zmax-zmin)/100:zmax;
neq=101;
M.setX(neq,z');
fite=M.fitFnc();
dlmwrite([outpath,'fit_eq_',fname],cat(2,z',fite),' ');

% plot corrected z-scale
zeta=param(5)*ag(i)/2/beta;
zscale=zprobeabs(nz,data(:,1)-z0,t,beta,mu/2/beta,b, c);
if (b~=0 || c~=0)
  zscale1=zprobeabs(nz,data(:,1)-z0,t,beta,mu/2/beta,0, 0);
  zsc=cat(1,zscale,zscale1);
  dlmwrite([outpath,'zscale_',fname],cat(2,zscale,zscale1),' ');
  leg={'homogeneous', 'general'};
  plotData([nz,nz],zsc,S,leg,' ',opt1);
else
  dlmwrite([outpath,'zscale_',fname],zscale,' ');
  leg={''};
  plotData(nz,zscale,S,leg,' ',opt1);
end

% plot inhomogeneity function (if defined)
if (b~=0 || c~=0)
nh=101;
fh=zeros(nh,2);
for j=1:nh
  z=(j-1)*t/(nh-1);
  fh(j,1)=z;
  fh(j,2)=1.0 + b*beta*z + c*(beta*z)^2;
end;
dlmwrite([outpath,'homog_',fname],fh,' ');
end


end;


%zmin=-2;
%zmax=12;
%z=zmin:(zmax-zmin)/70:zmax;
%h=z';
%tmp=SurfEffect(h,0,423.41941,0.93297,0.197,10,-0.0425,0.00568);
%dlmwrite([outpath,'shift_','Fe110_S14R.dat'],cat(2,h,tmp-tmp(35)),' ');
%tmp=SurfEffect(h,0,180.34703,0.97863,0,10,0.24481,-0.0238);
%dlmwrite([outpath,'shift_','Fe110_S14T.dat'],cat(2,h,tmp-tmp(35)),' ');
%tmp=SurfEffect(h,0,723.43405,0.88891,0.197,10,-0.0425,0.00568);
%dlmwrite([outpath,'shift_','Fe110_S25R.dat'],cat(2,h,tmp-tmp(35)),' ');
%tmp=SurfEffect(h,0,-147.50473,0.94381,0,10,0.24481,-0.0238);
%dlmwrite([outpath,'shift_','Fe110_S25T.dat'],cat(2,h,tmp-tmp(35)),' ');
%tmp=SurfEffect(h,0,401.18174,0.93111,0.1747,10,0.02813,-0.00347);
%dlmwrite([outpath,'shift_','Fe211_S14R.dat'],cat(2,h,tmp-tmp(35)),' ');
%tmp=SurfEffect(h,0,182.49185,0.96871,0,10,9.21E-4,0.00238);
%dlmwrite([outpath,'shift_','Fe211_S14T.dat'],cat(2,h,tmp-tmp(35)),' ');




