%% ENGINX
% Analytical model of ENGINX for spurious peak shift corrections

classdef ENGINX < handle 
    properties(Constant)        
    end
    properties
        MX; % MS_model_tof class instance
        EXT; % Extinction class instance
        par; % default parameters
        arange; % range of scattering angles (array of values) in [rad]
        src; % list of instrument components
        i_sam; % references to sample and detector objects
        i_det;
        theta; % mean Bragg angle [rad]
        gamma; % sample orientation angle
        lam0; % mean wavelength
        dhkl; % [A]
        thick; % sample thickness [mm]
    end
    
methods
       % constructor 
function X = ENGINX(flux,pulse,param,extinction)  
           pulse_table=readXY(pulse);
           flux_table=readXY(flux);
           X.par=readHash(param);
           X.MX = MX_model_tof(flux_table, pulse_table, extinction);
end
       
%% instrument setup initialization
% arange = range of scattering angles (array of values) in [deg]
% omega = sample rotation [deg] (=theta_S for symmetric reflection)
% dhkl = sample dhkl [A]
% thick = sample thickness [mm]
% s0w, s0d = divergence slit width and distance [mm]
% s1w, s1d = width and distance of the input slit [mm]
% s2w, s2d = width and distance of the output slit [mm]
function init(X, arange, omega, dhkl, thick, s0w,s0d,s1w,s1d,s2w)  
    % degree
    deg=pi/180;         
%% User input of instrument parameters        
    % DIVERGENCE SLIT
    X.par('s0_width')=s0w;    % width [mm]
    X.par('s0_dist')=s0d;   % distance [mm] from the sample

    % PRIMARY SLIT
    X.par('s1_width')=s1w;    % width [mm]
    X.par('s1_dist')=s1d;     % distance [mm] from the sample

    % SAMPLE
    X.gamma=(90+omega)*deg; % sample surface angle w.r.t. ki
    X.thick = thick;
    X.dhkl=dhkl;

    % RADIAL COLLIMATOR or SECONDARY SLIT
    X.par('s2_width')=s2w;   % gauge width [mm]

    % DETECTOR
    d_angle_min=arange(1)*deg;  % minimum angle (integration range)
    d_angle_max=arange(end)*deg; % maximum angle (integration range)   

    X.arange=arange*deg;
    % mean Bragg angle
    X.theta=0.25*(d_angle_min+d_angle_max);
    % mean wavelength
    X.lam0=2*X.dhkl*sin(X.theta);
    X.setParam();

end
       
function setParam(X)
 %% Creation of the instrument array
% Note the convention for distances input (1st parameter) below:
% Signs of distances: Before the sample, tracing is assumed in up-stream
% direction, therefore the distances are negative. 
% Distances refer to the previous/next component along the neutron beam.
    g_norm=1/sqrt(4*log(2));  % gaussian distribution
    g_uni=sqrt(1/6);  % uniform distribution
    i=0;
    % SOURCE
    % parameters: distance [mm], wavelength [A], tau [us], shape (see NOTE)
    if (X.par('s0_on')) 
        srcL=-X.par('src_dist')+X.par('s0_dist');
    else
        srcL=-X.par('src_dist')+X.par('gd_dist');  
    end
    p1=X.par('pulse_width');
    p2=X.par('pulse_shape')*g_norm;
    i=i+1;sr(i)=Pulse(srcL,p1, p2);
%
    i=i+1;
    % DIVERGENCE SLIT
    % parameters: distance [mm], wavelength [A], width[mm]
    if (X.par('s0_on')) 
        i=i+1;sr(i)=Slit(-X.par('s0_dist')+X.par('s1_dist'),X.par('s0_width'));
    else
        i=i+1;sr(i)=Guide(-X.par('gd_dist')+X.par('s1_dist'),X.par('gd_width'),X.par('gd_m'));    
    end
    % INPUT SLIT or RADIAL COLLIMATOR
    % parameters:  distance [mm], wavelength [A], width[mm]
    i=i+1;  
    sr(i)=Slit(-X.par('s1_dist'),X.par('s1_width'));


    % parameters: dhkl [A], surface angle [rad], thickness [mm]
    i=i+1;  sr(i)=Sample(X.dhkl,X.gamma,X.thick);
    X.i_sam=i;

    % OUTPUT SLIT or RADIAL COLLIMATOR
    % parameters: distance [mm], wavelength [A], width[mm]
    i=i+1;  
    sr(i)=Slit(X.par('s2_dist'), X.par('s2_width'));

    % DETECTOR
    i=i+1;  sr(i)=TofDetector(X.par('d_dist')-X.par('s2_dist'), ...
    X.par('d_binwidth'), X.par('d_binshape')*g_uni,X.par('d_resol'),...
    X.arange(1), X.arange(end));
    X.i_det=i;  
    X.src=sr;
           
    % set new instrument sequence and calculate T and C matrices
    X.MX.setParam(X.src,X.theta,X.i_sam,X.i_det); 
end

%% set theta
function setTheta(X,theta)
    X.theta=theta*pi/180;
    X.MX.setParam(X.src,X.theta,X.i_sam,X.i_det);
end
    
%% set gamma
function setGamma(X,gamma)
    X.gamma=gamma*pi/180;
    X.src(X.i_sam).gamma=X.gamma;
    X.MX.setParam(X.src,X.theta,X.i_sam,X.i_det);
end

%% gauge parameters averaged over the detector angle range
function [DSE, beta, zeta] = getGaugeParam(X)
        [xDSE,xbeta, xzeta, p] = X.MX.scan_theta(X.src,0.5*X.arange);
        DSE=xDSE'*p;   
        beta=xbeta'*p;
        zeta=xzeta'*p;
end

%% extinction parameters averaged over the detector angle range
function [res] = getAveragedParam(X)
        [res] = X.MX.getAverageParam(X.src,0.5*X.arange);
end

%% gauge parameters averaged over the detector angle range
% and peak parameters as a table(z,true_z,pshift,gwidth)
function [DSE, beta, zeta, res] = getPeakParam(X,z)
        d=X.src(X.i_sam).thickness;
        PS = PSModel();
        PS.setParam(d,0,0);
        [xDSE,xbeta, xzeta, p] = X.MX.scan_theta(X.src,0.5*X.arange);
        DSE=xDSE'*p;   
        beta=xbeta'*p;
        zeta=xzeta'*p;
        zscale=PS.get_ctr_ave(z',xbeta,xzeta, p);
        pshift=PS.get_shift_ave(z',xDSE,xbeta,xzeta, p);
        gwidth=PS.get_gwidth_ave(z',xbeta,xzeta, p);
        res=cat(2,z',zscale,pshift,gwidth);        
end
    
%% Resolution curve
function [ res] = getResolution(X,drange)
        res = X.MX.scan_dhkl_fwhm(X.src,drange); 
end
    
%% Scan sample orientation
% return gauge parameters for each gamma value in rows
function [ DSE, beta, zeta, p ] = scan_gamma(X,gamma)
        [DSE,beta, zeta, p] = X.MX.scan_gamma(X.src,gamma*pi/180);              
end
    
    
%% Scan dhkl
% return gauge parameters for each dhkl value in rows
function [ DSE, beta, zeta ] = scan_dhkl(X,dhkl)
        [DSE,beta, zeta] = X.MX.scan_dhkl(X.src,dhkl);             
end
       
    
end
    
end

