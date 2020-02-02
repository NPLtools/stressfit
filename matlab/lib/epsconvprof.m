function [ output ] = epsconvprof(d,S,a,omega,beta,zeta,dsam)
%epsconvprof Smeared depth profile of one strain component
% Assumes:
% a) isotropic material
% b) stress depth profile modeled as
% $$ \sigma(z) = (a_{1} + a_{2}z a_{3}z^2) exp(-2 omega z) $$
% h ... depth coordinate [mm]
% S .. one row of the compliance matrix (in units of E .. Young modulus)
% a .. dim=3 vector with stress profile coeficients 
% omega .. exponential decay coefficient of the stress profile
% dsam  .. sample thickness
% returns: depth profile of the smeared strain in microstrain units [10^-6]
%
% version: 1.2
% date: 4/4 2014
%
n=size(d,1);
output=zeros(n,1);
mm=zeros(3,1);
bb=[1  beta beta^2];
for k=1:n 
   h=d(k)*beta-zeta; 
   for j=1:3
     om=omega(j)/beta;     
     aa=a(j,:).*bb;
     mm(j)=sigconv(h,aa,om,dsam*beta);
   end;
   output(k)=S*mm*1e6;
end
