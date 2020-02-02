%%MX_model_tof
% Matrix model for a neutron ToF diffractometer.
% Provides gauge volume parameters needed for calculation of spurious strains.
% Written by: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
% Last update: Sept 2016
%------------------------------------------------------------------------
% Usage:
% 1) create an instance, X=MX_model_tof()
% 2) initialize for given setup parameters, X.setParam(inst,theta,i_sam,i_det)
%   inst = array of instrument components
%   theta = nominal Bragg angle [rad]
%   i_sam,i_det = indexes of the Sample and TofDetector components in
% 3) use functions scan_theta, scan_dhkl and get_average to calculate the
% gauge volume parameters:
% DX = spurious strain sensitivity 
% beta = gauge volume size parameter
% zeta = attenuation parameter (depth dependence)
% fwhm = instrument resolution (delta_d/d) 
% 

classdef MX_model_tof < handle
    properties(Constant)
%% constants
% uniform distribution                 
g_uni=sqrt(1/6);
% circular distribution                 
g_circ=sqrt(1/8);
% triangular distribution
g_tri=sqrt(1/3);
% gaussian distribution
g_norm=1/sqrt(4*log(2));
r8ln2=sqrt(8*log(2));
% degree
deg=pi/180;
% h/m [Ang*m/ms]
hm=3.95603;
%   
    end
    properties  

%% local variables
    DeltaX; % spurious strain sensitivity, factor at the shift function, eq.(14)
    beta;   % gauge volume size parameter
    zeta;   % attenuation term in eq. (12), zeta=a(1)/2beta
	thick;  % sample thickness
    a;      % attenuation, depth dependence coefficients
    mu;     % attenuation coefficient [1/cm]
    m;      % index of the time variable
    k;      % index of the depth variable
    T;      % total transmission matrix
    C;      % total transport matrix
    isam;   % index to the sample
    idet;   % index to the detector
    lam0;   % nominal wavelength
    Ltof;   % path length
    tof0;   % nominal flight tim [us]
    dtau;   % detector time resolution [us]
    pulsetab; % source pulse width table (2-columns with wavelength [A] and pulse width [us])
    fluxtab;  % source flux table (2-columns with wavelength [A] and flux)
    ext;     % Extinction object to calculate attenuation cross-sections
    
    end
    
    methods(Static)
        function [res] = g_ref(theta)
        % sample rotation angle - symmetric reflection
           res=theta+pi/2; 
        end        
        function [res] = g_tra(theta)
        % sample rotation angle - symmetric transmission
           res=theta+pi; 
        end
                
        % Normalized gauge volume distribution
        % z = coordinate along surface nomal
        % u = depth
        % t = layer thickness
        % all in normalized length units
        function [res] = gauge_dist(z,u,t)
            res=2/sqrt(pi)*exp(-z^2)/(erf(u)-erf(u-t));        
        end
        
        % SGV center as a function of the nominal depth, u
        % t = thickness, all in normalized length scale
        function [res] = zscale(u,t)
            n=length(u);
            res=zeros(n,1);
            for i=1:n            
              res(i)=(u(i)-f(u(i),t));
            end            
        end
        
        % width (fwhm) of the sampled gauge volume 
        % u = nominal depth
        % t = thickness, all in normalized length scale
        function [res] = gwidth(h,t)  
            n=length(h);
            res=zeros(n,1);
            for i=1:n 
              u=h(i);
              fu=f(u,t);
              gd=MX_model_tof.gauge_dist(u-t,u,t);
              sig2=0.5 + u*fu-fu^2-0.5*t*gd;
              res(i)=sqrt(sig2);
            end;
        end

    end
    
    methods 
        function G = MX_model_tof(flux, pulse, ext) 
            G.m=5;
            G.k=2;    
            G.pulsetab=pulse;
            G.fluxtab=flux;
            if (nargin>2) 
                G.ext=ext;
            else 
                % default attenuation cross section coefficients
                G.ext=Extinction([0.954 0.033, 0.127, 7.82e-2]);
            end;
            %G.pulsetab=[0.1 80.0;20 80.0];
            %G.fluxtab=[0.1 1;20 1];
        end
        
        % get tables provided by user
       % function setFlux(flux)
       %     G.fluxtab=flux;
       % end
       % function setPulse(pulse)
       %     G.pulsetab=zeros(size(pulse,1),size(pulse,2));
       %     G.pulsetab(:,:)=pulse(:,:);
       % end
        
        % set new instrument sequence and calculate T and C matrices
        function setParam(G,inst,theta,i_sam,i_det)            
             n=length(inst); 
             G.Ltof=0.0;  
             G.T=zeros(5,5);             
             G.lam0=2*inst(i_sam).dhkl*abs(sin(theta));
             G.isam=i_sam;
             G.idet=i_det;
			 G.thick=inst(i_sam).thickness;
             % get pulse width for given wavelength
             inst(1).tau=interp1q(G.pulsetab(:,1),G.pulsetab(:,2),G.lam0)*MX_model_tof.g_norm;
             % start initialization from sample:
             inst(i_sam).initialize(eye(5,5),G.lam0);
             % inst(i_sam).display();
             C1=inst(i_sam).getCin();
             C2=inst(i_sam).getCout();
             % sample -> detector
             for i=i_sam+1:n
                 inst(i).initialize(C2,G.lam0);  
                 TT=inst(i).getT();
                 G.T = G.T + TT;                 
                 C2 =  inst(i).getCout();
                 G.Ltof = G.Ltof+inst(i).L;
             end
             % sample -> source
             for i=i_sam-1:-1:1
                 inst(i).initialize(C1,G.lam0);
                 TT=inst(i).getT();
                 G.T = G.T + TT;
                 C1 =  inst(i).getCout();
                 G.Ltof = G.Ltof-inst(i).L;
             end
             G.tof0=G.Ltof*G.lam0/MX_model_tof.hm;
             G.C=C2;
             G.a=inst(G.isam).a;
             G.dtau=inst(G.idet).dtau;
             G.init();
        end    
         
        % Check transmission and transport matrices. Return true if they are
        % valid (det|T|>0 and det|C|=1). 
        % Calculate dependent arrays and width parameters
        % it = index of time-coordinate
        % iz = index of z-coordinate
        % Calculate beta and the factor at <x_D>, eqs. (9) and (14)
        function res = init(G)
           DT=det(G.T);
           DC=det(G.C);
           if (DT<1e-8) || (abs(DC-1)>1e-6)
               res=false;
               return
           end
           res=true;
           % eq. (7), define W and V as submatrices of T
           W=G.T;         
           %fprintf('init m,k=%d,%d\n',ix,iz);
           %fprintf('W rank=%d\n',rank(W));           
           W(G.k,:)=[];
           W(:,G.k)=[];        
           V=G.T(:,G.k);
           V(G.k,:)=[];
           WV=W\V; % = inv(W)*V
           beta2=G.T(G.k,G.k)-V'*WV;
           %fprintf('beta2=%f\n',beta2);
          % display(W);
           %display(V);
           Cm=G.C(G.m,:);
           Cm(:,G.k)=[];
           G.beta=sqrt(beta2);           
           G.DeltaX=-1.0e6*(G.C(G.m,G.k)-Cm*WV)/G.tof0; 
           G.mu=G.ext.getMu(G.lam0);
           G.zeta=0.1*G.a(1)*G.mu/(2*G.beta);
        end
        
        % width (fwhm) of the sampled gauge volume [mm]
        function [res] = gauge_width(G)
            res=MX_model_tof.r8ln2/G.beta;         
        end
        
        % peak width in delta_d/d 
        function [res] = peak_width(G)
            TC=G.C*(G.T\G.C');
            % add detector resolution
            res=1e2*sqrt(TC(G.m,G.m)+G.dtau^2)/G.tof0/MX_model_tof.g_norm;
        end

        % scan gauge parameters for a range of Bragg angles
        % assumes unchanged setup (src) sequence and a previous initialization by
        % setParam. p contains weight factors derived from the flux table
        function [DX, beta, zeta, p] = scan_theta(G,src,theta)
            n=length(theta);
            beta=zeros(n,1);
            DX=zeros(n,1);
            zeta=zeros(n,1);
            th0=src(G.isam).theta; 
            p=zeros(n,1);
            for i=1:n                
              th=theta(i);
              G.setParam(src,th,G.isam,G.idet); 
              p(i)=interp1q(G.fluxtab(:,1),G.fluxtab(:,2),G.lam0);
              beta(i)=G.beta;
              DX(i)=G.DeltaX;
              zeta(i)=G.zeta;
            end 
            p0=sum(p);
            p=p/p0;
            % reset to the previous state
            G.setParam(src,th0,G.isam,G.idet); 
        end
        
        % Return array with parameters averaged over angular range
        % res = [DX,beta,zeta, ag,cg, mu]
        function [res] = getAverageParam(G,src,theta)
            n=length(theta);
            p_DX=zeros(n,1);
            p_beta=zeros(n,1);
            p_zeta=zeros(n,1);
            p_ag=zeros(n,1);
            p_cg=zeros(n,1);
            p_mu=zeros(n,1);
            th0=src(G.isam).theta; 
            p=zeros(n,1);
            for i=1:n                
              th=theta(i);
              G.setParam(src,th,G.isam,G.idet); 
              p(i)=interp1q(G.fluxtab(:,1),G.fluxtab(:,2),G.lam0);
              p_DX(i)=G.DeltaX;
              p_beta(i)=G.beta;
              p_zeta(i)=G.zeta;
              p_ag(i)=G.a(1);
              p_cg(i)=G.a(2);
              p_mu(i)=G.mu;
            end 
            p0=sum(p);
            p=p/p0;
            res(1)=p_DX'*p;
            res(2)=p_beta'*p;
            res(3)=p_zeta'*p;
            res(4)=p_ag'*p;
            res(5)=p_cg'*p;
            res(6)=p_mu'*p;
            % reset to the previous state
            G.setParam(src,th0,G.isam,G.idet); 
        end
        
        % scan resolution for a range of Bragg angles
        % assumes unchanged setup (src) sequence and a previous initialization by
        % setParam
        function [fwhm, p] = scan_theta_fwhm(G,src,theta)
            n=length(theta);
            fwhm=zeros(n,1);
            th0=src(G.isam).theta; 
            p=zeros(n,1);
            for i=1:n                
              th=theta(i);
              G.setParam(src,th,G.isam,G.idet);
              p(i)=interp1q(G.fluxtab(:,1),G.fluxtab(:,2),G.lam0);
              fwhm(i)=G.peak_width();
            end
            p0=sum(p);
            p=p/p0;
            % reset to the previous state
            G.setParam(src,th0,G.isam,G.idet); 
        end

 
        % scan gauge parameters for a range of sample orientations
        % assumes unchanged setup (src) sequence and a previous initialization by
        % setParam
        function [DX, beta, zeta, p] = scan_gamma(G,src,gamma)
            n=length(gamma);
            beta=zeros(n,1);
            DX=zeros(n,1);
            zeta=zeros(n,1);
            g0=src(G.isam).gamma; 
            th0=src(G.isam).theta; 
            p=zeros(n,1);
            for i=1:n                
              src(G.isam).gamma=gamma(i);              
              G.setParam(src,th0,G.isam,G.idet); 
              p(i)=interp1q(G.fluxtab(:,1),G.fluxtab(:,2),G.lam0);
              beta(i)=G.beta;
              DX(i)=G.DeltaX;
              zeta(i)=G.zeta;
            end 
            p0=sum(p);
            p=p/p0;
            % reset to the previous state
            src(G.isam).gamma=g0;
            G.setParam(src,th0,G.isam,G.idet); 
        end
         

        % scan gauge parameters for a range of dhkl
        % assumes unchanged setup (src) sequence and a previous initialization by
        % setParam
        function [DX, beta, zeta] = scan_dhkl(G,src,dhkl)
            n=length(dhkl);
            beta=zeros(n,1);
            DX=zeros(n,1);
            zeta=zeros(n,1);
            na=11;
            th1=src(G.idet).amin/2;
            th2=src(G.idet).amax/2;
            da=(th2-th1)/(na-1);
            theta=th1:da:th2;
            th0=src(G.isam).theta; 
            d0=src(G.isam).dhkl;
            for i=1:n                
              src(G.isam).dhkl=dhkl(i);      
              [xDX, xbeta, xzeta, p]=scan_theta(G,src,theta);
              DX(i)=xDX'*p;
              beta(i)=xbeta'*p;
              zeta(i)=xzeta'*p;
            end 
            % reset to the previous state
            src(G.isam).dhkl=d0;
            G.setParam(src,th0,G.isam,G.idet); 
        end
        
        % scan resolution for a range of dhkl
        % assumes unchanged setup (src) sequence and a previous initialization by
        % setParam
        function [ fwhm ] = scan_dhkl_fwhm(G,src,dhkl)
            n=length(dhkl);
            fwhm=zeros(n,1);
            na=11;
            th1=src(G.idet).amin/2;
            th2=src(G.idet).amax/2;
            da=(th2-th1)/(na-1);
            theta=th1:da:th2;
            th0=src(G.isam).theta; 
            d0=src(G.isam).dhkl;
           % fprintf('scan_dhkl_fwhm:\ndhkl, lmin, lmax, fwhm\n');
            for i=1:n                
              src(G.isam).dhkl=dhkl(i);
              [w, p]=scan_theta_fwhm(G,src,theta);
              fwhm(i)=w'*p;
              %lam1=2*dhkl(i)*sin(theta(1));
              %lam2=2*dhkl(i)*sin(theta(na));
             % fprintf('%f %f %f  %f\n',dhkl(i),lam1,lam2,fwhm(i));
            end 
            % reset to the previous state
            src(G.isam).dhkl=d0;
            G.setParam(src,th0,G.isam,G.idet); 
        end
        
        % get gauge parameters averaged over the detector range
        % assumes unchanged setup (src) sequence and a previous initialization by
        % setParam
        function [DX, beta, zeta] = get_average(G,src)
            n=11;
            theta=zeros(n,1);
            th1=src(G.idet).amin/2;
            th2=src(G.idet).amax/2;
            dth=(th2-th1)/(n-1);
            for i=1:n
                theta(i)=th1+(i-1)*dth;
            end;
            [xDX, xbeta, xzeta] = G.scan_theta(src,theta);           
            DX=mean(xDX);
            beta=mean(xbeta);
            zeta=mean(xzeta);
        end

    end
end

