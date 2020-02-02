function [ output ] = strain(npart,z,S,a,omega)
%strain Principal components of the strain tensor as a function of depth
% Calculates depth profile of strain assuming 
% stress components in the principal directions, i, modeled as
% $$ \sigma_{j}(z) = (a_{j,1} + a_{j,2}z a_{j,3}z^2) exp(-2 omega_j z) $$
% z .. depth under surface
% S .. compliance matrix (in units of E .. Young modulus)
% a(j,k) .. 3x3 matrix with stress profile coeficients (in units of E)
% omega .. dim=3 vector with exponential decay coefficients 
%
% NOTE: The model uses absolute length units [mm]
%
% version: 1.2
% date: 4/4 2014
%


ndim=numel(npart);

% if ndim=1, return the components in columns 1..3
% otherwise, return everything in one sequence with given partitioning

if (ndim==1)
    nd=size(z,1);
    ndat=[nd nd nd];
    x=cat(1,z,z,z);
    nc=3;
    nmax=3*nd;
else
    nd=min(size(z,1),sum(npart));
    ndat=npart;
    x=z;
    nc=ndim;
    nmax=nd;
end;


e=zeros(nmax,1);
mm=zeros(3,1);
% loop through components
n1=1;
for i=1:nc
  n2=min(n1-1+ndat(i),nmax);
for k=n1:n2
    zz=x(k);   
    for j=1:3
       d=exp(-2*omega(j)*zz)*1e6;
       mm(j)=d*(a(j,1)+a(j,2)*zz+a(j,3)*zz^2);
    end
    e(k)=S(i,:)*mm;     
end
  n1=n1+ndat(i);
end
if (ndim==1)
    output=reshape(e,nd,3);
else
    output=e;
end;

