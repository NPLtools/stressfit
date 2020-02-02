%% Static functions and global constants
classdef Utils
    %% Distribution shape factors and other constants
    properties(Constant)
        % uniform distribution                 
        g_uni=sqrt(1/6);
        % circular distribution                 
        g_circ=sqrt(1/8);
        % triangular distribution
        g_tri=sqrt(1/3);
        % gaussian distribution
        g_norm=1/sqrt(4*log(2));
        % degree
        deg=pi/180;
    end
    
    methods (Static)
        
        %% Analytical function describing dependence of peak shift vs. surface position
        % $$f(x)  =  \zeta - \pi ^{-{1/2}}\frac{  \exp(-{(x-\zeta )}^{2})}
        % {1 + erf(x-\zeta )}  $$ 
        % see equation (12)
        % 
        % $x$ = distance from surface (< 0 outside the sample), 
        % $\zeta$ = attenuation factor
        %
        function [res] = fsh(x,zeta)
            res = zeta - sqrt(pi)*exp(-(x-zeta)^2)/(1+erf(x-zeta)) ;   
        end    
        
        %% Generate default plot optios to be used with |myfigure.m|.
        function [opt] = plotOptions()
            opt.groups=[1 2 3 4 5];
            opt.colors=[[0 0 0]' [1 0 0]'  [0 0 1]' [0 1 0]' [0.5 0 0.5]' [0.5 0.5 0]' [0 0.5 0.5]' ]';
            opt.styles={'-' '--' '-.' ':' };
            opt.XGrid='on';
            opt.YGrid='on';
            opt.XMinorTick='on';
            opt.YMinorTick='on';
            opt.xauto=true;
            opt.yauto=true;
            opt.xlimits=[0 10];
            opt.ylimits=[0 10];
            opt.LineWidth=1;
            opt.xlabel='X-axis';
            opt.ylabel='Y-axis';
            opt.LegendPosition='Best';
            opt.ptsize=3;
            opt.ptstyles = {'o' 'x' '*' '.' '+'};  
            opt.showlegend=1;
        end;
        
        %% Generates legend as a set of strings from numerical data.
        %
        function [leg] = getLegend(data, format)
            n=length(data);
            leg=cell(1,n);    
            for i=1:n
                s=sprintf(format,data(i));
                leg(i)={s};
            end;
        end
        
    end
    
end

