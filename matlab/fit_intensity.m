function [ fit par epar z_scale chi2 ] = fit_intensity(inpath,outpath, prefix, fname, cfg, par0, fix)
%fit_intensity encapsulates the script for fitting intensity scans 
% suitable for use in an automation script for processing many data files
%
% parameters:
% inpath = input path
% outpath = subdirectory where the output should go
% prefix = string prefix for input filename, log file name and output directory
% fname = input file name WITHOUT EXTENSION (.dat is assumed).
% cfg(1) = ag (linear coefficent for attenuation depth dependence)
% cfg(2) = cg (constant coefficent for attenuation depth dependence)
% cfg(3) = sign of front surface orientation: +1 if the scan direction is from outside into the
% par0 = initial parameters value
% fix = flag for fixed parameters
%
% Returns:
% fit = fitted intensity curves
% par = fitted parameters.
% epar = parameters error estimates.
% zeta = mu*a/2/beta derived from the 
% z_scale = depth scale conversion 
%           array(encoder, corrected for z0, gauge centre, dtto without inhomogeneity)
% version: 1.1
% date: 20/6/2017
%


% create output directory if necessary
[status, msg, msgID] = mkdir(outpath); %#ok<NASGU,ASGLU>

% prefix for input/output files
logfile=cat(2,prefix,'_',fname,'.log');

% The script can fit multiple data sets in a single runx
% Provide 3-column data: encoder position [mm], intensity, error
inpfile=cat(2,prefix,'_',fname,'.dat');

ag=cfg(1);
cg=cfg(2);
xsgn=cfg(3);
% front surface orientation: xsng=+1 if the scan direction from outside into the
% bulk is positive


%% provide initial parameters for each data file (one column per input file).
% The order of rows is:
% param=[y0, I0, z0, beta, mu, d, b, c]'
% y0   ... background
% I0   ... intensity
% z0   ... surface position [mm]
% beta ... gauge volume size parameter [mm-1]
% mu   ... Attenuation coefficient [1/cm], estimate
% d    ... sample thickness [mm] 
% b    ... linear inhomogeneity coefficient
% c    ... quadratic inhomogeneity coefficient
param=par0;
% provide flags indicating which parameters are fixed (0|1)
fixed=fix;
 
%% Fit the intensity scan

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

S=strrep(inpfile,'_','\_');
hdr=sprintf('Fitting file: %s',inpfile);
fprintf('Fitting file: %s',inpfile);
% load data
data1=readXYZ(strcat(inpath,inpfile));
% select the range here ...
data = data1(1:size(data1,1),:);
nz=size(data,1);
% set geometry factors for absorption
M.ag=ag*0.1;
M.cg=cg*0.1;

% run fit
M.setParam1(param,fixed);
M.setX(nz,data(:,1),xsgn);
[ output,errors,fit,chi2,Cov ] = fit_model(data,M,1); %#ok<NASGU>
par=M.a;
epar=errors;
[y0, I0, z0, beta, mu, t, b, c] = M.getResult(1); %#ok<ASGLU>

% print report
M.reportResults(flog,1, hdr,output,errors,chi2);

% calculate fit without inhomogeneity correction
if (b~=0 || c~=0)
  param1=output';
  param1(7:8)=0.0;
  M.setParam1(param1,fixed);
  [ output,errors,fith,chi20,Cov ] = fit_model(data,M,0); %#ok<NASGU,ASGLU>
  % fith = M.fitFnc();  
end

% plot fit
leg={'data','fit','fit without inhomogeneity'};
if (b~=0 || c~=0)
  fdata=cat(1,data(:,1:2),cat(2,data(:,1),fit(:,1)),cat(2,data(:,1),fith(:,1)));
  result=cat(2,data(:,1)-z0,data(:,1:2),fit(:,1),fith(:,1));
  ndat=[nz nz nz];
else
  fdata=cat(1,data(:,1:2),cat(2,data(:,1),fit(:,1)));
  result=cat(2,data(:,1)-z0,data(:,1:2),fit(:,1));
  ndat=[nz nz];
end
plotData(ndat,fdata,S,leg,' ',opt);
% plotFit(nz,data(:,1),data(:,2),nz,data(:,1),fit,S,leg,'measurement','fit',opt);
out=cat(2,outpath,prefix,'_',fname,'_fit.dat');
dlmwrite(out,result,' ');

% save fit on equidistant scale
zmin=data(1,1);
zmax=data(nz,1);
z=zmin:(zmax-zmin)/100:zmax;
neq=101;
M.setX(neq,z',xsgn);
fite=M.fitFnc();
out=cat(2,outpath,prefix,'_',fname,'_fit_eq.dat');
dlmwrite(out,cat(2,z',fite),' ');

% plot corrected z-scale

zeta=M.ag*mu/2/beta;
zscale=zprobeabs(nz,data(:,1)-z0,t,beta,zeta,b, c,xsgn);
out=cat(2,outpath,prefix,'_',fname,'_zscale.dat');
if (b~=0 || c~=0)
  zscale1=zprobeabs(nz,data(:,1)-z0,t,beta,zeta,0, 0,xsgn);
  zsc=cat(1,zscale,zscale1);
  z_scale=cat(2,data(:,1),zscale,zscale1(:,2));
  dlmwrite(out,z_scale,' ');
  leg={'homogeneous', 'general'};
  plotData([nz,nz],zsc,S,leg,' ',opt1);
else
  z_scale=cat(2,data(:,1),zscale);
  dlmwrite(out,z_scale,' ');
  leg={''};
  plotData(nz,zscale,S,leg,' ',opt1);
end

% plot inhomogeneity function (if defined)
if (b~=0 || c~=0)
  nh=101;
  fh=zeros(nh,2);
  for j=1:nh
    z=(j-1)*t/(nh-1);
    fh(j,1)=z*xsgn;
    fh(j,2)=1.0 + b*beta*z + c*(beta*z)^2;
  end;
  out=cat(2,outpath,prefix,'_',fname,'_homog.dat');
  dlmwrite(out,fh,' ');
  opt1.xlabel='true depth, d [mm]';
  opt1.ylabel='inhomogeneity weight factor';
  plotData(nh,fh,S,leg,' ',opt1);
end;

end

