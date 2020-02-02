function [ val ] = G3(x,u,a,d)
%G1 Integral of the gauge volume distribution function
%%
%
% $$ G3(z, h, d) = [3a(a-x) + x^2 + 1] \frac {1}{2} gn(x,u,d) + (a^2 + \frac {3}{2}) aG0(x,u,d)$$
%
% $$ G3(x,u,a,d) \equiv \int { (a-x)^3 gn(x) } d x $$
%
% $$ gn(x) \equiv \frac {2} {\sqrt(\pi}}  \frac {e^{-x^2}} {erf(u) - erf(u-d} $$
%
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
val =(3*a*(a-x)+x^2+1)*0.5*gn(x,u,d) + (a^2+1.5)*a*G0(x,u,d);


