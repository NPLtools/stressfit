%% Script for fitting through surface scans - intensity
%
%
close all;
clear all;
deg=pi/180;
addpath('lib','libtools');
% Define input / outut paths
inpath='./input/int/';
outpath='./output/';
% Input files:
inpTable=[inpath,'int_scan_input.txt']; 
outTable=[outpath,'int_scan_input.out'];

%% Instrument setup
% range of scattering angles to average over
arange=75:5:105;
% sample thickness
d_sam=10;
% collimation
s0w=14;
s0d=4000;
s1w=2;
s1d=60;
s2w=2;
% extinction
extab=readXY('./input/extinction.dat');
extObj=Extinction([0 0.034, 0.121, 7.82e-2],extab);
%% Initialize
% Create object implementing analytical calculation of the gauge parameters
% It is specific for each instrument, but encapsulates the general implementation
% of the matrix model in MX_model_tof.m
INS=ENGINX('./input/spectrum.dat','./input/pulse.dat','./input/ENGINX.par',extObj);

% Initial parameters for each data file (one column per input file).
% The order is:
% param=[y0, I0, z0, beta, mu, d, b, c]'
% y0   ... background
% I0   ... intensity
% z0   ... surface position [mm]
% beta ... gauge volume size parameter [mm-1]
% mu   ... Attenuation coefficient [1/cm], estimate
% d    ... sample thickness [mm] 
% b    ... linear inhomogeneity coefficient
% c    ... quadratic inhomogeneity coefficient
par0=[1, 3000, -5.0, 1.0, 1.147, 10, 0, 0]';
% flags indicating which parameters are fixed (0|1)
fix0=[1, 0, 0, 0, 1, 1, 0, 0]';

%% Read input table
% read list of input data in the format (name, dhkl, omega, s0w, s0d, z0, xsign)
% xsign = scan sign: positive if increasing position goes under surface
% z0 = estimate position of the front surface (facing primary beam)
% Output: names = ID, data = (dhkl, omega, s0w, s0d, z0, xsign);
[ names data ] = readTable(inpTable,6);
ninp=size(names,1);

% prepare output table
outtab=cell(ninp,1);

%% start loop with input data
for i=1:ninp
  close all;
  % decompose name to get directory, prefix and filename
  ss=names{i};
%  s2=strrep(ss,'_',' ');
%  s3=textscan(s2, '%s');
%  s31=s3{1,1};
%  dir=s31{1,1};
%  prefix=['int_' s31{1,1} '_'  s31{2,1}];
%  fname=[sprintf('%s_',s31{3:end-1}),s31{end}];
  dir='int';
  prefix='int';
  fname=ss;
% get paraneters from the input table 
  dhkl=data(i,1);
  omega=data(i,2);
  s0w=data(i,3);
  s0d=data(i,4);
  z0=data(i,5);
  xsgn=data(i,6);
  
  % FUNCTION for calculating gauge parameters analytically
  INS.init(arange, omega, dhkl, d_sam, s0w, s0d, s1w, s1d, s2w);
  res = INS.getAveragedParam();
  DSE=res(1);
  beta=res(2);
%  zeta=res(3);
  ag=res(4);
  cg=res(5);
  mu=res(6);
  % prepare input for fitting
  % 
  cfg=[ag,cg,xsgn];
  fix=fix0;
  par0(3)=z0;
  par0(4)=beta;
  par0(5)=mu;
  if (omega>90)
      fix(5)=1;
  end;
  % FUNCTION for intensity fit:
  [ fit par epar z_scale chi2] = fit_intensity(inpath,[outpath,dir,'/'],prefix, fname, cfg, par0, fix);
  zeta=0.1*mu*ag/2/beta;
  % add output table row
  tmp=sprintf('%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g',ss,chi2,dhkl,omega,DSE,beta,zeta,par);
  outtab{i}=sprintf('%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g',tmp,epar);
end

%% save result table
fout=fopen(outTable,'wt'); 
fprintf(fout,'name\tchi2\tdhkl\tomega\tDSE\tbeta\tzeta\ty0\tI0\tz0\tbeta\tmu\td\tb\tc\t');
fprintf(fout,'err_y0\terr_I0\terr_z0\terr_beta\terr_mu\terr_d\terr_b\terr_c\n');
for i=1:size(outtab,1)
    fprintf(fout,'%s\n',outtab{i,1});
end;
fclose(fout);


  
  