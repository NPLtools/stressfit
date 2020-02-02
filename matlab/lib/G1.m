function [ val ] = G1(x,u,a,d)
%G1 Integral of the gauge volume distribution function
%%
%
% $$ G1(z, h, d) = a G0(x,u,d) + \frac {1}{2} gn(x,u,d) $$
%
% $$ G1(x,u,a,d) \equiv \int { (a-x) gn(x) } d x $$
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
val = a*G0(x,u,d)+0.5*gn(x,u,d);


