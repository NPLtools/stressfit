function ReportInit( filename, header )
%REPORTINIT Writes logfile header in a new file

if (~isempty(filename))
    fid=fopen(filename,'wt');
    fprintf(fid,'Date: %s\n',datestr(now));
    fprintf(fid,'%s\n',header);
    fprintf(fid,'%s\n','------------------------------------------------');
    fclose(fid);
end;

end

