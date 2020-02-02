function [out] = zprobe(h,d)
%zprobe Conversion from nominal position to the true probe position
% Calculates the true depth of the center of mass of the gauge volume 
% for given nominal depth
% Works in reduced length coordinates (units of beta^-1).
% h ... reduced depth coordinate: $$ h = \beta depth - \zeta $$
% d ... sample thickness , $$ d = \beta thickness $$
%
% version: 1.2
% date: 4/4 2014
%
n=size(h,1);
n2=size(h,2);
out=zeros(n,n2);
for j=1:n2
 for k=1:n
    xx=h(k,j);
	dd=d(j);
    out(k,j)=xx-f(xx,dd);
 end
end
end
