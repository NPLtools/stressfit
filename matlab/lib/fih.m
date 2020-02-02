%% Shift function
% Analytical function describing the difference between true and nominal gauge volume center of mass
% Version with inhomogeneity correction
% 
% $x$ = distance from surface (< 0 outside the sample), 
% $d$ = sample thickness
% ihb, ihc = inhomogeneity coefficients (f=1 + ihb*z + ihc*z^2)
%
% version: 1.3
% date: 2/8/2016
%

function [ val ] = fih( x, d, ihb, ihc )
% approximation for large x
 if (x-d>5)
     val=x-d;
 elseif (x<-5)
     val=x;
 else        
    f0=F0(x,d);
    f1=F1(x,d);
    if (ihb==0 && ihc==0)
        val = x - f1/f0;
    else
        f2=F2(x,d);
        f3=F3(x,d);
        val = x - (f1+ihb*f2+ihc*f3)/(f0+ihb*f1+ihc*f2);
    end
    % zjm= erf(x) - erf(x-d);
    % val =  -1/sqrt(pi)*(exp(-x^2)-exp(-(x-d)^2))/zjm ;
 end;
end

