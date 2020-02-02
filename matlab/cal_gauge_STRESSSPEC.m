function [ beta DSE ag cg G src] = cal_gauge_STRESSSPEC(theta_M, theta_S, omega, thick, mu, s1w,s1d,s2w,s2d )
% Calculates gauge parameters analytically 
% theta_M = monochromator Bragg angle [deg]
% theta_S = sample Bragg angle [deg]
% omega = sample rotation [deg] (=theta_S for symmetric reflection)
% thick = sample thickness [mm]
% mu = absotption coefficient [cm^-1]
% s1w, s1d = width and distance of the input slit [mm]
% s2w, s2d = width and distance of the output slit [mm]
%
% RETURNS:
% beta = gauge width parameter [1/mm]
% DSE = peak shift coefficient [microstrain/mm]
% ag, cg = geometrical factors describing beam path as a function of scan depth, z: 
%          path/thick = ag*z[mm] + cg; 
% G = MAT_general object
% src = list of components

format shortG;
addpath('lib','libmat_mon','libtools');


%% Distribution shape factors and other constants
% gaussian distribution
g_norm=1/sqrt(4*log(2));
% degree
deg=pi/180;
% Set excluded parameters to infinity
INF=1e10;
% dhkl values
Si_400=1.3576;

%% User input

% Instrument parameters
% Components before the monochromator are defined directly 
% in the instrument setup section

% MONOCHROMATOR
takeoff_M=-2*theta_M;   % take-off angle [deg]
dM=Si_400;              % d-spacing [A]
radM=5.5;               % radius [m-1]
m_dist=2100;            % monochromator-sample distance
m_thick=11.6;           % monochroator thickness
m_chi=0;                % monochroator cutting angle [deg], 0 for symmetric reflection

% DIVERGENCE SLIT (after monochromator)
% included only if s0_on=1 (divergence slit is used)
s0_width=INF;    % width [mm]
s0_dist=950;    % distance [mm] from the sample
s0_on=0;        % = 1 to include this slit

% RADIAL COLLIMATOR or PRIMARY SLIT
s1_width=s1w;      % width [mm]
s1_dist=s1d;       % distance [mm] from the sample

% SAMPLE
gamma=(90+omega)*deg; % sample surface angle w.r.t. ki

% RADIAL COLLIMATOR or SECONDARY SLIT
s2_width=s2w;     % gauge width [mm]
s2_dist=s2d;      % set to zero for a radial collimator

% DETECTOR
d_binwidth=2;       % resolution fwhm [mm]
d_binshape=g_norm;  % uni for uniform, norm for gaussian bin
d_dist=935;        % distance from the sample

%% Instrument setup
% Note the convention for distances input (1st parameter) below:
% Signs of distances: Before the sample, tracing is assumed in up-stream
% direction, therefore the distances are negative. 
% Distances refer to the previous/next component along the neutron beam.

% wavelength
lambda=2*dM*sin(abs(takeoff_M*deg/2));
% dhkl
dhkl=lambda/2/sin(abs(theta_S*deg));

i=0;
% SOURCE - comment out if not included
% parameters: distance [mm], width[mm]
%i=i+1;src(i)=Slit(-3097.0,80);

% SOURCE SLIT - comment out if not included
% parameters:  distance [mm], width[mm], divergence[rad]
i=i+1;  src(i)=Slit(-1150.0,50);

% SOURCE COLLIMATOR - comment out if not included
% parameters:  distance [mm], width[mm], divergence[rad]
%i=i+1;  src(i)=Soller(-1150.0,INF,INF);

% NEUTRON GUIDE - comment out if not included
% parameters:  distance [mm], width[mm], m-value, wavelength [A]
% i=i+1;  src(i)=Guide(-500,50,3,lambda);  

% MONOCHROMATOR
% parameters:  distance [mm], thickness [mm], tak-off angle [rad], cutting angle [rad], curvature [1/mm], Poisson constant
if (s0_on==1) 
   LM=m_dist-s0_dist;
else
   LM=m_dist-s1_dist; 
end
i=i+1;  src(i)=CrystalBent(-LM,m_thick,takeoff_M*deg,m_chi*deg,0.001/radM,0.3);
i_mon=i;

% DIVERGENCE SLIT
% parameters: distance [mm], width[mm]
if (s0_on==1)
i=i+1;  src(i)=Slit(-(s0_dist-s1_dist), s0_width);
end;

% INPUT SLIT or COLLIMATOR
% parameters: distance [mm], width[mm]
i=i+1;  src(i)=Slit(-s1_dist, s1_width);

% SAMPLE
% parameters: Bragg angle [rad], surface angle [rad], attenuation coeff. [mm^-1]
i=i+1;  src(i)=Sample(theta_S*deg,gamma,mu*0.1,thick);
i_sam=i;
fprintf('Sample reflection, dhkl=%f, 2 thetaS=%f, gamma=%f\n',dhkl, 2*theta_S, gamma/deg);

% OUTPUT SLIT
% parameters: distance [mm], width[mm]
i=i+1;  src(i)=Slit(s2_dist, s2_width);

% DETECTOR
% parameters: distance [mm], resolution [mm]
LD=d_dist-s2_dist; % detection distance
i=i+1;  src(i)=Detector(LD, d_binwidth*d_binshape/g_norm);

%----------------------------------------
% End of instrument definition
%----------------------------------------
fprintf('wavelength = %f\n',lambda);
% Create the object which implements the matrix method
G=MAT_general();
% set new instrument sequence and calculate T and C matrices
G.setParam(src,i_sam,i_mon,3,4);

%% Results
% geometrical factors describing beam path as a function of scan depth,z
% path/thick = ag*z[mm] + cg;
cos1=cos(gamma);
cos2=cos(2*theta_S*deg-gamma);
ag=0.1*(1/cos2-1/cos1);
cg=0.1*0.5*(sign(cos2)/cos2 + sign(cos1)/cos1 - ag*10);
% conversion detecor coordinate -> microstrain
shift_fact=-1e6/2/LD/tan(theta_S*deg);
% peak shift factor (NOTE: oposite sign convention is used in MAT_general) 
DSE=-G.DeltaX*shift_fact;
% gauge size parameter
beta=G.beta;

end

