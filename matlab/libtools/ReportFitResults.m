function  ReportFitResults( filename, log2scr, header,names, values, errors, fixed )
%ReportFitResults Print fit results
% Use empty filename for output in the screen window
% set log2scr=1 to write also on the screen
% version: 1.2
% date: 4/4/2014
%


what=zeros(1,2);
if (log2scr)
    what(1)=1;
end;

if (~isempty(filename))
    what(2)=1;
end;

% ------ start loop with log output ------
for iw=1:2
if (what(iw)>0)


if (iw==1)
   fid=1;
else
   fid=fopen(filename,'at'); 
end;

if (~isempty(header))
  fprintf(fid,'%s\n',header);
end;

n= size(fixed,2);
nf=0;
nv=0;
for i=1:n
    if (fixed(i))
        nf=nf+1;
    else
        nv=nv+1;
    end;
end;

if (nv>0) 
 fprintf(fid,'Fitted variables:\n');
for i=1:n
    if (~ fixed(i))
       fprintf(fid,'\t%s = %f +- %f\n',names{i},values(i),errors(i));
    end;
end;
end;

if (nf>0)
 fprintf(fid,'Fixed variables:\n');
for i=1:n
    if (fixed(i))
        fprintf(fid,'\t%s = %f\n',names{i},values(i));
    end;
end;
end;

if (~ (fid==1))
   fclose(fid);
end;

end;
end;
% ------ end loop with log output ------


