function [ output ] = sigconv(h,a,omega,d)
%sigconv Principal components of the stress tensor as a function of depth
% smeared by the gauge volume
% For the surface effect and meaning of the zeta,beta parameters, see
% J. Saroun et al. J. Appl. Cryst. (2013) 46, 628-638.
% The true stress components in the principal directions, j, are modeled as
% $$ \sigma(h) = (a_{1} + a_{2}h a_{3}h^2) exp(-2 omega h) $$
% h ... reduced depth coordinate: $$ h = \beta depth - \zeta $$
% a .. vector with the 3 stress profile coeficients 
% omega .. dim=3 vector with exponential decay coefficients of the stress profile
% d  .. sample thickness
% returns: value of the smeared strain at given depth
% NOTE: if the dimension of beta, zeta and S(:,1)=1, then only one strain component is calculated.
% This function assumes dimensionless length unit (scaled in $\beta^{-1}$)
% version: 1.2
% date: 4/4 2014
%
	 u=h-omega;
     m0=M0(h,omega,d);
     m1=u-f(u,d);
     m2=0.5 + u*(u-f(u,d));
	 ex=2*u/d-1;
	 if (abs(ex)<0.01)
		ad=exp(-(d/2)^2)/2/sqrt(pi)/erf(d/2)*(1+u*d-d^2/2);
	 else
		ad=f(u,d)/(1-exp(d*(d-2*u)));
	 end;
	 m2=m2+ad;
     output=m0*(a(1)+a(2)*m1+a(3)*m2);
end

