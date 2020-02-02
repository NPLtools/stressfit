function [output] = zprobeabs(npart,h,d,beta,zeta,ihb, ihc, xsgn)
%zprobeabs Conversion from nominal position to the true probe position
% Calculates the true depth of the center of mass of the gauge volume 
% for given nominal depth
% Works in absolute length coordinates (unit=mm)
% For the surface effect and meaning of the zeta,beta parameters, see
% J. Saroun et al. J. Appl. Cryst. (2013) 46, 628-638.
% h ... nominal (encoder) depth coordinate [mm]
% d ... sample thickness [mm]
% beta  .. gauge volume size parameter [1/mm]
% zeta  .. attenuation coefficient , zeta = a/2/beta
% ihb, ihc .. inhomgeneity coefficients, f = 1 + ihb*z + ihc*z^2
% xsgn ... optional array with x-scale signs (one for each beta value)
% 
% Updated: added inhomogeneity terms 
% version: 1.3
% date: 2/8/2016
%


ndim=numel(npart);

% if ndim=1, return the components in columns 1..3
% otherwise, return everything in one sequence with given partitioning


if (ndim==1)
    nd=size(h,1);
    nc=numel(beta);
    ndat=zeros(1,nc)+nd;
    x=h;
    for i=2:nc
        x=cat(1,x,h);
    end;
    nmax=nc*nd;
else
    nd=min(size(h,1),sum(npart));
    ndat=npart;
    x=h;
    nc=ndim;
    nmax=nd;
end;
if (nargin<8)
    xsg=ones(1,nc);
else
    xsg=xsgn;
end;

if (nargin<7)
    ihcc=zeros(1,nc);
else
    ihcc=ihc;
end;

if (nargin<6)
    ihbb=zeros(1,nc);
else
    ihbb=ihb;
end;

zdep=zeros(nmax,2);
% loop through components
n1=1;
for j=1:nc
  n2=min(n1-1+ndat(j),nmax);
  for k=n1:n2
    bb=beta(j);
    xx=xsg(j)*bb*x(k)-zeta(j);
    dd=bb*d;
    zdep(k,1)=x(k);    
    zdep(k,2)=xsg(j)*PSModel.gcentre(xx,dd,ihbb(j),ihcc(j))/bb;
    % zdep(k,2)=(xx-f(xx,dd))/bb;
  end
  n1=n1+ndat(j);
end
if (ndim==1)
    output=reshape(zdep,nd,nc*2);
else
    output=zdep;
end;
