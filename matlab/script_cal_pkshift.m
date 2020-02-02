%% Calculate shifts for the gauge parameters given as a table 
close all;
clear all;

%% Definitions
addpath('lib','libtools');
inpath='./input/';
outpath='./output/';
deg=pi/180;


%% User Input
%---------------------------------------------------------
% parameters
% assumed format is: name,DSE,beta,zeta,thickness,b,c
ftable='param_sim.dat';
% prefix and suffix for output files
% output file name
output='cal_shift_sim.dat';
% scan range 
zmin=-2;
zmax=11.8;
nd=70;
z=zmin:(zmax-zmin)/(nd-1):zmax;
% subtract d0? set 0 or 1
sub=0;

%% read input table
% data = [ DSE	beta zeta	d	b	c	z0 ]
[ hrows data ] =readTable(strcat(inpath,ftable),7);
nr=size(data,1);

M=PSModel();
out=z';
for j=1:nr
    M.setParam(data(j,4),data(j,5),data(j,6));
    tmp=M.get_shift(z',data(j,1),data(j,2),data(j,3));
    out=cat(2,out,tmp-sub*tmp(fix(nd/2)));
end;

nd=size(out,1);
fid = fopen([outpath,output], 'w');
fprintf(fid, '# calculated spurious microstrains\n');
fprintf(fid, '# parameter file: %s\n', ftable);
fprintf(fid, '%s\t%s\t', 'z', hrows{:});
fprintf(fid,'\n');
for row=1:nd
    fprintf(fid,'%f\t',out(row,:));
    fprintf(fid,'\n');
end
fclose(fid);



