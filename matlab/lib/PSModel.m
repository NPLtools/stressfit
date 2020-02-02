%% PSModel
% Gaussian model for calculation of spurious strains due to perturbation of 
% sampling gauge volume (e.g. surface effect).
% Calculates spurious strains and gauge senter as a function of scan depth 
% for given gauge volume parameters:
%    DSE = spurious strain sensitivity [mm-1]
%    beta = gauge volume size parameter [mm-1]
%    zeta = attenuation parameter (depth dependence)
% Use the matrix model (like MX_model_tof) to get these parameters
% analytically.
%
% Written by: J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
% Last update: Sept 2016
%------------------------------------------------------------------------
% Usage:
% 1) create an instance, X=PSModel()
% 2) obtain the gauge parameters DSE, beta, zeta for given setup
%    a) analytically: DSE, beta, zeta are calculated by MX_model_tof
%    b) experimentally: beta and zeta can be fitted on intensity vs. depth scan
%    c) simulation: DSE, beta can be obtained by fitting simulated depth scans
% 3) use functions get_shift and get_ctr to calculate spurious stran and 
% sampling centre of mass for given DSE, beta, zeta and an array of depth values.
% NOTE: Give DSE, beta, zeta as arrays to calculate multiple shift functions
% simultaneously.
% 

classdef PSModel < handle
    properties(Constant)  
    end
    properties  
%% local variables
   % sample thickness [mm]   
      dsam; 
   % inhomogeneity parameters, assumig quadratic weight function:
   %  f(u) = 1 + b*u + c*u^2       
      b;
      c;
    end
    
methods(Static)

    function [res] = meanp(A,p)
       res=A'*p; 
    end
    
% Fn(u,d) functions are n-th moments of the gauge volume distribution:
% Integral{z^n exp(-(z-u)^2)}dz from 0 to d
% in relative length units:
% h=depth*beta, u = h - zeta
% zeta=a/2/beta
% d=thickness*beta

    function [ val ] = F0( u, d )
        val=sqrt(pi)/2*(erf(u)-erf(u-d));
    end    

    function [ val ] = F1( u, d )
        val=0.5*(exp(-u^2)-exp(-(u-d)^2)) + u*PSModel.F0(u,d);
    end

    function [ val ] = F2( u, d )
        val= 0.5*PSModel.F0(u,d) + u*PSModel.F1(u,d) - d*0.5*exp(-(u-d)^2) ;
    end

    function [ val ] = F3( u, d )
        val= 0.5*u*PSModel.F0(u,d) + (u^2+1)*PSModel.F1(u,d) - 0.5*d*(d+u)*exp(-(u-d)^2) ;
    end

% width (fwhm) of the sampled gauge volume 
% u = nominal depth
% d = thickness, all in normalized length scale
	function [res] = gwidth(u,d)
        f2=PSModel.F2(u,d);
        f1=PSModel.F1(u,d);
        f0=PSModel.F0(u,d);
        sig2=f2/f0-(f1/f0)^2;
        res=sqrt(sig2);
    end
    
    
% volume of the sampled gauge volume = intensity
% u = nominal depth
% d = thickness, all in normalized length scale
	function [res] = gvol(u, d, zeta, b, c)
        f0=PSModel.F0(u,d);
        val=f0;
        if (b~=0)
          f1=PSModel.F1(u,d);
          val = val + b*f1;
        end;
        if (c~=0)
          f2=PSModel.F2(u,d); 
          val = val + c*f2;
        end;        
        res=exp(-zeta*(2*u+zeta))*val;        
    end
         

% Centre of mass of the sampling volume.
% Includes inhomogeneity correction.
% Assuming relative length units, 1/beta    
% u = h - zeta, where 
%    h is distance from the front surface in (< 0 outside the sample), 
%    zeta is the attenuation coefficient (incl. geometry factor)
% d = sample thickness
% b, c = inhomogeneity coefficients: p = 1 + b*u + c*u^2
    function [ val ] = gcentre( u, d, b, c )
        % approximation for large u
        if (u-d>5)
            val=d;
        elseif (u<-5)
            val=0;
        else        
            f0=PSModel.F0(u,d);
            f1=PSModel.F1(u,d);            
            if (b==0 && c==0)
                val = f1/f0;
            else
                f2=PSModel.F2(u,d);
                f3=PSModel.F3(u,d);
                val = (f1+b*f2+c*f3)/(f0+b*f1+c*f2);
            end
        end;
    end

end
    
methods 
	function G = PSModel() 
	end
        
% set sample parameters:
% d = thickness [mm]
% b,c = inhomogeneity parameters, assumig quadratic weight function:
%       f(u) = 1 + b*u + c*u^2 
	function setParam(G,dsam,b,c)
        G.dsam=dsam;
        G.b=b;
        G.c=c;
	end;

% Spurious strain as a function of scan depth (h>0 is under front surface)
% Parameters
% h ... depth values [mm] (absolute position)
% DSE ... surface effect amplitude [micrstrain/mm]
% beta ... gauge volume size parameter
% zeta ... attenuation parameter
	function [ eshift ] = get_shift(G,h,DSE,beta,zeta)
        n=size(h,1);
        nb=size(DSE,1);
        eshift=zeros(n,nb);
        for i=1:n
        for j=1:nb
            hh=h(i)*beta(j)-zeta(j);
            d=G.dsam*beta(j);
            eshift(i,j)=DSE(j)*(PSModel.gcentre(hh,d,G.b,G.c)/beta(j) - h(i));
        end
        end
    end;

    function [ eshift ] = get_shift_ave(G,h,DSE,beta,zeta,p)
        n=size(h,1);
        nb=size(DSE,1);
        es=zeros(n,nb);
        for i=1:n
        for j=1:nb
            hh=h(i)*beta(j)-zeta(j);
            d=G.dsam*beta(j);
            es(i,j)=DSE(j)*(PSModel.gcentre(hh,d,G.b,G.c)/beta(j) - h(i));
        end
        end
        eshift=zeros(n,1);
        for i=1:n
            %eshift(i)=sum(es(i,:))/nb;
            eshift(i)=es(i,:)*p;
        end
    end;
    
% Centre of mass of the sampling volume as a function of scan depth (h>0 is under front surface)
% Parameters
% h ... depth values [mm] (absolute position)
% beta ... gauge volume size parameter
% zeta ... attenuation parameter
	function [ ctr ] = get_ctr(G,h,beta,zeta)
        n=size(h,1);
        nb=size(beta,1);
        ctr=zeros(n,nb);
        for i=1:n
        for j=1:nb
            hh=h(i)*beta(j)-zeta(j);
            d=G.dsam*beta(j);
            ctr(i,j)=PSModel.gcentre(hh,d,G.b,G.c)/beta(j);
        end
        end
    end;
    function [ ctr ] = get_ctr_ave(G,h,beta,zeta,p)
        n=size(h,1);
        nb=size(beta,1);
        xctr=zeros(n,nb);
        for i=1:n
        for j=1:nb
            hh=h(i)*beta(j)-zeta(j);
            d=G.dsam*beta(j);
            xctr(i,j)=PSModel.gcentre(hh,d,G.b,G.c)/beta(j);
        end
        end
        ctr=zeros(n,1);
        for i=1:n
            % ctr(i)=sum(xctr(i,:))/nb;
            ctr(i)=xctr(i,:)*p;
        end
    end; 
        
% Gauge width size (fwhm) as a function of scan depth (h>0 is under front surface)
% Parameters
% h ... depth values [mm] (absolute position)
% beta ... gauge volume size parameter
% zeta ... attenuation parameter
	function [ gwidth ] = get_gwidth(G,h,beta,zeta)
        n=size(h,1);
        nb=size(beta,1);
        gwidth=zeros(n,nb);
        for i=1:n
        for j=1:nb
            hh=h(i)*beta(j)-zeta(j);
            d=G.dsam*beta(j);
            gwidth(i,j)=PSModel.gwidth(hh,d)/beta(j);
        end
        end
    end;   
    function [ gwidth ] = get_gwidth_ave(G,h,beta,zeta, p)
        n=size(h,1);
        nb=size(beta,1);
        gw=zeros(n,nb);
        for i=1:n
        for j=1:nb
            hh=h(i)*beta(j)-zeta(j);
            d=G.dsam*beta(j);
            gw(i,j)=PSModel.gwidth(hh,d)/beta(j);
        end
        end
        gwidth=zeros(n,1);
        for i=1:n
           %  gwidth(i)=sum(gw(i,:))/nb;
            gwidth(i)=gw(i,:)*p;
        end
    end;   
    end
end

