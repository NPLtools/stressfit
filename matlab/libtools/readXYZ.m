function [ data ] = readXYZ( filename )
%readXYZ Read 3-column data from an ASCII file
fid = fopen(filename);
data = fscanf(fid, '%g %g %g',[3 inf])';
fclose(fid);
end

