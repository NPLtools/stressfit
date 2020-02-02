function [ val ] = G2(x,u,a,d)
%G1 Integral of the gauge volume distribution function
%%
%
% $$ G2(z, h, d) = (2a-x) \frac {1}{2} gn(x,u,d) + (a^2 + \frac {1}{2}) G0(x,u,d)$$
%
% $$ G2(x,u,a,d) \equiv \int { (a-x)^2 gn(x) } d x $$
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
val = (2*a-x)*0.5*gn(x,u,d) + (a^2+0.5)*G0(x,u,d);


