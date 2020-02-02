function [ eps ] = strain_seg(z,S,L,sig)
%strain_seg Segmented stress profile model
% Depth (z) profile is defined by straight lines on given intervals.
% It is therefore defined by:
% L(i) ... node positions along z [mm]
% sig(i) ... stress values at the nodes [0.001/E]
% ns ... number of segments
% Free variables: ns-1 L values and ns+1 sig values
% Other arguments:
% S ... Compliance tensor (or its submatrix)
% dsam ... sample thickness [mm]
%
% version: 1.2
% date: 7/4 2014
%

% number of data points
nz=size(z,1);
% number of branches
nbr=size(sig,2);
% number of S rows (strain components)
ne=size(S,1);
% number of S columns (stress components)
ns=size(S,2);
eps=zeros(nz,ne);
if (nbr ~= ns)
    fprintf('strain_seg: compliance matrix (%d x %d) is not compatible with the number of stress components (%d) !\n', ne,ns,nbr);    
    return;
end;

sigma=zeros(nz,ns);

for i=1:nbr
  sigma(:,i)=stress_seg(z,L,sig(:,i));
end;
for i=1:nz
    eps(i,:)=S*sigma(i,:)';
end;


