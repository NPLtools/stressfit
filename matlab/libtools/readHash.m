function [  hash ] = readHash( filename )
%% readHash reads a simple hash table in the format <name> <value>
fid = fopen(filename);
cells = textscan(fid, '%s %f','CollectOutput', 1, 'CommentStyle', '#')';
keys=cells{1};
d=cells{2};
values=d(:,1);
fclose(fid);
hash = containers.Map(keys, values);

end

