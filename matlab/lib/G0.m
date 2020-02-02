function [ val ] = G0(x, u,  d )
%G0 Integral of the gauge volume distribution function
%%
%
% $$ G0(x, h, d) = \frac{erf(x)}{erf(u) - erf(u-d)} $$
%
% $$ G0(z, h, d) \equiv \int {gn(x) } d x $$
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
 zjm=erf(u)-erf(u-d);
 if (abs(zjm)>0)
    val=erf(x)/zjm;
 else
    val=0;
 end

