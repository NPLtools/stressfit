function [ output ] = epsconv(ncom,d,S,a,omega,beta,zeta,dsam)
%epsconv Principal components of the strain tensor as a function of depth
% smeared by the gauge volume, without (!) the surface effect
% For the surface effect and meaning of the zeta,beta parameters, see
% J. Saroun et al. J. Appl. Cryst. (2013) 46, 628-638.
% Calculates depth profile of strain assuming 
% a) isotropic material
% b) stress components in the principal directions, i, modeled as
% $$ \sigma_{j}(z) = (a_{j,1} + a_{j,2}z a_{j,3}z^2) exp(-2 omega_j z) $$
% d ... depth coordinate [mm]
% S .. compliance matrix (in units of E .. Young modulus)
% a .. 3x3 matrix with stress profile coeficients 
% omega .. dim=3 vector with exponential decay coefficients of the stress profile
% dsam  .. sample thickness
% returns: depth profile of the smeared strain in microstrain units [10^-6]
%
% ncom ... vector of integers: partitioning of d into strain components
% size(ncom) gives the number of components and also dimension of beta,
% zeta, ...
% ncom(i) = number of consecutive d-values corresponding to the i-th strain component 
%
% version: 1.2
% date: 4/4 2014
%

ndim=numel(ncom);

% if ndim=1, return the components in columns 1..3
% otherwise, return everything in one sequence with given partitioning

if (ndim==1)
    nd=size(d,1);
    ndat=[nd nd nd];
    x=cat(1,d,d,d);
    nc=3;
    nmax=3*nd;
else
    nd=min(size(d,1),sum(ncom));
    ndat=ncom;
    x=d;
    nc=ndim;
    nmax=nd;
end;

econv=zeros(nmax,1);
mm=zeros(3,1);

n1=1;
% loop through components
for i=1:nc    
   bb=[1  beta(i) beta(i)^2];
   n2=min(n1-1+ndat(i),nmax);
% loop through data points
for k=n1:n2 
   h=x(k)*beta(i)-zeta(i); 
   for j=1:3
     om=omega(j)/beta(i);     
     aa=a(j,:).*bb;
     mm(j)=sigconv(h,aa,om,dsam*beta(i));
   end;
   econv(k)=S(i,:)*mm*1e6;
end;
   n1=n1+ndat(i);
end;

if (ndim==1)
    output=reshape(econv,nd,3);
else
    output=econv;
end;


