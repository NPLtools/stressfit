function [ output ] = epsconv_seg(ncom,d,S,L,sigma,beta,zeta,dsam)
%epsconv_seg Principal components of the strain tensor as a function of depth
% smeared by the gauge volume, without (!) the surface effect
% Depth (z) sress profile is defined by straight lines on given intervals.
% It is defined by:
% L(i) ... node positions along z [mm]
% sig(i) ... stress values at the nodes [0.001*E]
% Free variables: ns-1 L values and ns+1 sig values
% Other arguments:
% S ... compliance matrix [units of E]
% dsam  .. sample thickness [mm]
% returns: depth profile of the smeared strain in microstrain units [10^-6]
%
% ncom ... vector of integers: partitioning of d into strain components
% ncom(i) = number of consecutive d-values corresponding to the i-th strain component 
% size(ncom) gives the number of strain components = nc
% number of rows of beta and zeta must equal nc
% Dimension of S must be compatible with the number of required stress and strain
% components
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
% initialize result array
econv=zeros(nmax,1);


% number of branches
nbr=size(sigma,2);
% number of nodes
nsig=size(sigma,1);
% number of segments
ns=nsig-1;
% number of inner node positions
nL=numel(L);
if (ns ~= nL+1)
    if (ndim==1)
        output=reshape(econv,nd,3);
    else
    output=econv;
    end;
    fprintf('epsconv_seg: Dimensions of L and sig arrays do not agree: nL=%d, nsig-2=%d\n',nL,nsig-2);
    fprintf('Expected: nL = nsig-2\n');
    return;
end;
if (size(S,2) ~= nbr)
    fprintf('epsconv_seg: Dimensions of S and sig do not agree: cols(S)=%d, cols(sig)=%d\n',size(S,2),nbr);
    return;
end;
if (size(S,1) ~= nc)
    fprintf('epsconv_seg: Dimensions of S and sig do not agree: cols(S)=%d, cols(sig)=%d\n',size(S,1),nc);
    return;
end;

% partitionling into segments of the interval (0,dsam)
nL=numel(L);
zs=zeros(nL+2,1);
zs(2:nL+1)=L(:);
zs(nL+2)=dsam;

n1=1;
% loop through components
for i=1:nc    
   n2=min(n1-1+ndat(i),nmax);
   % convert zs and dsam to reduced length scale
   dd=beta(i)*dsam;
   zsn=beta(i)*zs;
% loop through data points
for k=n1:n2 
   % convert depth positions to reduced length scale 
   h=x(k)*beta(i)-zeta(i); 
   mm=sigconv_seg(h,zsn,sigma,dd);
   econv(k)=S(i,:)*mm*1e6;
end;
   n1=n1+ndat(i);
end;

if (ndim==1)
    output=reshape(econv,nd,3);
else
    output=econv;
end;


