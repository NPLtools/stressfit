function [ sigma ] = stress_seg(z,L,sig,dsam)
%stress_seg Segmented stress profile model
% Depth (z) profile is defined by straight lines on given intervals.
% It is therefore defined by:
% L(i) ... inner node positions along z [mm]
% sig(i) ... stress values at the nodes [0.001*E]
% dsam ... sample thickness [mm]
% Number of sig rows determines the number of segments, ns=size(sig,1)-1.
% Number of sig columns determines the number of components.
% L(i) are only the inner nodes, the first and last ones have positions 0
% and dsam, respectively.
% version: 1.2
% date: 7/4 2014
%

% number of data points
nd=size(z,1);
% number of branches
nbr=size(sig,2);
% number of nodes
nsig=size(sig,1);
% number of segments
ns=nsig-1;
% number of inner node positions
nL=numel(L);
% initialize result array
sigma=zeros(nd,nbr);
if (ns ~= nL+1)
    fprintf('stress_seg: Dimensions of L and sig arrays do not agree: nL=%d, nsig-2=%d\n',nL,nsig-2);
    fprintf('Expected: nL = nsig-2\n');
    return;
end;

% partitionling into segments of the interval (0,dsam)
zs=zeros(nL+2,1);
zs(2:nL+1)=L(:);
zs(nL+2)=dsam;
% calculate gradients
alpha=zeros(ns,nbr);
for k=1:nbr
for i=1:ns
    alpha(i,k)=(sig(i+1,k)-sig(i,k))/(zs(i+1)-zs(i));
end;
end;
for i=1:nd
  zz=z(i); 
  for k=1:nbr
  for is=1:ns
      if (zz>=zs(is) && zz<=zs(is+1))
          sigma(i,k)=sig(is,k)+alpha(is,k)*(zz-zs(is));
      end;
  end;  
  end;
end;



    
