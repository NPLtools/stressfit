function [ val ] = M0(h, omega, d)
%
% $$ M0(h,\omega, d)  =  \frac{ erf(h -\omega -d) - erf(h -\omega )} {erf(h -d) - erf(h)} \exp(-\omega (2(h-\zeta) - \omega)) $$ 
%
% h ... reduced depth coordinate: $$ h = \beta depth - \zeta $$
% d ... sample thickness , $$ d = \beta thickness $$
% omega ... exponential coefficient of the stress depth profile function (in units of beta)
% version: 1.2
% date: 4/4 2014
%
 val=(erf(h-omega-d) - erf(h - omega))/(erf(h-d) - erf(h))*exp(-omega*(2*h-omega));
end

