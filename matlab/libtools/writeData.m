function  writeData( filename, header, ndat,data )
%writeData Write data into an ASCII file
% given in the 'data' argument as 2-column array X,Y. 
% A 3rd column with errors is optional.
% The data can contain multiple sets, partition into these sets is then defined
% by the ndat argument.
%
% Use empty filename for output in the screen window
% version: 1.2
% date: 17/4/2014
%


% analyze data structure
nrow=max(ndat);
nset=numel(ndat);
ncol=min(3,size(data,2));
if (ncol>2)
    fmt='%f\t%f\t%f';
    labels='X%d\tY%d\tERR%d';
else
    fmt='%f\t%f';
    labels='X%d\tY%d';
end;

% initialize text lines
lines=cell(nrow,1);
for i=1:nrow
  lines{i}='';
end;

i1=1;
for k=1:nset
    i2=i1-1+ndat(k);
    for i=i1:i2
       s=sprintf(fmt,data(i,1:ncol)); 
       if (k<nset)
           s=sprintf('%s\t',s);
       end;
       iL=i-i1+1;
       lines{iL}=cat(2,lines{iL},s);
    end;
    i1=ndat(k)+1;
end;


% open file
if (~isempty(filename))
    fid=fopen(filename,'wt'); 
else
    fid=1;
end;


% write header
if (~isempty(header))
  fprintf(fid,'# %s\n',header);
end;

% write labels
s='';
for k=1:nset
    if (ncol>2)
         s1=sprintf(labels,k,k,k);
    else
         s1=sprintf(labels,k,k);
    end;
    s=cat(2,s,s1);
    if (k<nset)
        s=sprintf('%s\t',s);
    end;        
end;
fprintf(fid,'%s\n',s);
%write table
for i=1:nrow
  fprintf(fid,'%s\n',lines{i});
end;

% close file
if (~ (fid==1))
   fclose(fid);
end;


