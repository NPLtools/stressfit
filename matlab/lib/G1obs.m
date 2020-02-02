function [ val ] = G1obs(z, h, zeta,  d )
%G1 Integral of the gauge volume distribution function
%%
%
% $$ G1(z, h, d) = \frac{ \zeta  erf(z)  - \frac{1}{\sqrt{\pi} e^-{z^2} } {erf(h) - erf(h-d)} $$
%
% $$ G1(z, h, d) \equiv \int_{-\inf}^{z} { u g(u) \PHI(u) } d u $$
%
% $$ g(u) \equiv \frac {e^{-u^2}} {erf(h) - erf(h-d} $$
%
% where PHI(u)=1 inside the interval (h-d,h)
% Assumes reduced length units:
%
% $$ h = depth \beta - \zeta $$
% $$ z = z' \beta - \zeta $$
% $$ d = \beta thickness $$
% 
%
% version: 1.2
% date: 8/4/2014
%
 zjm=erf(h)-erf(h-d);
 if (abs(zjm)>0)
    val=(zeta*erf(z)-exp(-z^2)/sqrt(pi))/zjm;
 else
    val=0;
 end

