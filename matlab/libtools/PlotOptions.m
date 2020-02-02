classdef PlotOptions
%% Data structure for common figure options
% only line plots
    
    properties
       title='';
       xlabel=' ';
       ylabel=' ';
       xlimits=[0,10];
       ylimits=[0,10];
       colors=[[0 0 0]' [1 0 0]'  [0 0 1]' [0 1 0]' [0.5 0 0.5]' [0.5 0.5 0]' [0 0.5 0.5]' ]';
       styles={'-' '--' '-.' ':'};
       autox=true;
       autoy=true;
       gridx=true;
       gridy=true;
       ptsize=3;
       ptstyles = {'+' 'o' '*' '.' 'x'};  
       showlegend=1;
    end
    
    methods 
             
    end
    
end

