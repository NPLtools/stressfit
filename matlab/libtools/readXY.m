function [ data ] = readXY( filename )
%readXYZ Read 2-column data from an ASCII file
fid = fopen(filename);
data = fscanf(fid, '%g %g',[2 inf])';
fclose(fid);
end

