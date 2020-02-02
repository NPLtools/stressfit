%% Gauge centre function
% Centre of mass of the sampling volume.
% Includes inhomogeneity correction
% 
% u = h - zeta, where 
%    h is distance from the front surface in (< 0 outside the sample), 
%    zeta is the attenuation coefficient (incl. geometry factor)
% d = sample thickness
% b, c = inhomogeneity coefficients: p = 1 + b*z + c*z^2
%
% Assuming relative length units, 1/beta
%
% version: 1.3
% date: 2/8/2016
%

function [ val ] = gcentre( u, d, b, c )
% approximation for large u
 if (u-d>5)
     val=d;
 elseif (u<-5)
     val=0;
 else        
    f0=F0(u,d);
    f1=F1(u,d);
    if (b==0 && c==0)
        val = f1/f0;
    else
        f2=F2(u,d);
        f3=F3(u,d);
        val = (f1+b*f2+c*f3)/(f0+b*f1+c*f2);
    end
    % zjm= erf(x) - erf(x-d);
    % val =  -1/sqrt(pi)*(exp(-x^2)-exp(-(x-d)^2))/zjm ;
 end;
end

