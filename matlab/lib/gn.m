function [ val ] = gn(x, u,  d )
%G0 Integral of the gauge volume distribution function
%%
%
% $$ gn(x) \equiv \frac {2} {\sqrt(\pi}}  \frac {e^{-x^2}} {erf(u) - erf(u-d} $$
% Assumes reduced length units:
%
% $$ u = depth \beta - \zeta $$
% $$ x = z \beta - \zeta $$
% $$ d = \beta thickness $$
% 
%
% version: 1.3
% date: 28/2/2015
%
 if ((u>-5) && (u<d+5))
    v0=2/sqrt(pi)/(erf(u)-erf(u-d));
 else
    if (u<=-5)
      c=u;
    else
      c=u-d;
    end;
    v0=2*c/(exp(-(u-d)^2)-exp(-u^2));
 end
 val=v0*exp(-x^2);

