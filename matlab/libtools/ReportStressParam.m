function  ReportStressParam( filename, log2scr, header, values, errors )
%ReportStressParam Print fit results for the stress parameters
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

fprintf(fid,'Fitted values:\na [x1000] =\n');
fprintf(fid,'  %f    %f    %f\n',values(1,:)*1e3);
fprintf(fid,'  %f    %f    %f\n',values(2,:)*1e3);
fprintf(fid,'  %f    %f    %f\n',values(3,:)*1e3);
fprintf(fid,'omega = \n');
fprintf(fid,'  %f    %f    %f\n',values(4,:));

fprintf(fid,'Errors:\na [x1000] =\n');
fprintf(fid,'  %f    %f    %f\n',errors(1,:)*1e3);
fprintf(fid,'  %f    %f    %f\n',errors(2,:)*1e3);
fprintf(fid,'  %f    %f    %f\n',errors(3,:)*1e3);
fprintf(fid,'omega = \n');
fprintf(fid,'  %f    %f    %f\n',errors(4,:));

if (~ (fid==1))
   fclose(fid);
end;

end;
end;
% ------ end loop with log output ------


