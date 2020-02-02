function [ eshift ] = SurfEffect(h,y0,DSE,beta,zeta,dsam,b,c)
%SurfEffect Calculate surface effect corrections
% Parameters
% h ... depth values [mm] (absolute position)
% y0 ... strain offset [microstrain], e.g. due to misalignment or aberation 
% DSE ... surface effect amplitude [micrstrain/mm]
% beta ... gauge volume size parameter
% zeta ... attenuation parameter
% dsam ... sample thickness
% b,c ... linear and quadratic inhomogeneity coefficients
%

n=size(h,1);
eshift=zeros(n,1);
bb=0;
cc=0;
if (nargin>6)
    bb=b;
end
if (nargin>7)
    cc=c;
end


for i=1:n
	hh=h(i)*beta-zeta;
	d=dsam*beta;
	eshift(i)=y0 + DSE*(PSModel.gcentre(hh,d,bb,cc)/beta - h(i));
end

