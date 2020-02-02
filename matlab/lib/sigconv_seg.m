function [ Y ] = sigconv_seg(h,zs,sig,d)
%sigconv_seg Value of smeared stress profile for given depth
% Depth (z) sress profile is defined by straight lines on given intervals.
% zs(i) ... node positions along z [1/beta units]
% sig(i) ... stress values at the nodes [0.001*E]
% Dimension of zs determines the number of segments, ns=numel(zs)-1
% Other arguments:
% h ... gauge depth (in 1/beta units)
% d ... sample thickness (in 1/beta units)
% Works for a single branch only. Provide L and sig as (:,1) arrays.
% version: 1.2
% date: 7/4 2014
%

% number of nodes
nz=numel(zs);
nsig=size(sig,1);
% number of branches
nbr=size(sig,2);
% number of segments
ns=nz-1;
% initialize result array
Y=zeros(nbr,1);
if (nsig ~= nz)
    fprintf('sigconv_seg: Dimensions of zs and sig arrays do not agree: nz=%d, nsig=%d\n',nz,nsig);
    return;
end;

% calculate gradients
alpha=zeros(ns,nbr);
for k=1:nbr
for i=1:ns
    alpha(i,k)=(sig(i+1,k)-sig(i,k))/(zs(i+1)-zs(i));
end;
end;
GG0=zeros(nsig,1);
GG1=zeros(nsig,1);
for i=1:nsig
    GG0(i)=G0(h-zs(i),h,d);
    GG1(i)=G1(h-zs(i),h,d);
end;

for k=1:nbr
for i=1:ns
    a=GG0(i)-GG0(i+1);
    b=GG1(i)-GG1(i+1);
    Y(k)=Y(k)+a*(sig(i,k)+alpha(i,k)*(h-zs(i))) - b*alpha(i,k);
end;
end;


