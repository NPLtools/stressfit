function [ hrows data ] = readTable(filename, ncols)
%readTable Read numeric table with row headers
fid = fopen(filename);
cells = textscan(fid, '%s %f %f %f %f %f %f %f','CollectOutput', 1, 'CommentStyle', '#')';
hrows=cells{1};
d=cells{2};
data=d(:,1:ncols);
fclose(fid);
end

