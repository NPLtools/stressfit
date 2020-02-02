%%Guide
% Matrix model of a neutron guide
% Only the exit window width and critical angle are taken as constraints,
% assuming 100% reflectivity.
% 

classdef Guide < Comp
    properties(Constant) 
        GammaNi = 1.731e-3; %  Ni critical angle [rad/A]
    end
    properties  
        w;  % width [mm] 
        m;  % m-value
    end

    methods              
        function X=Guide(distance,width, m)
           X = X@Comp(distance);
           X.w=width*Comp.g_uni;     
           X.m=m;
           X.id='Guide';
        end        
    
        % Transmission matrix
        function [T] =getT(X)
           V=X.Cin(1,:);
           T=getT@Comp(X);
           if (X.w>0)
             T=T+(X.w^-2)*kron(V',V);
           end
           if (X.m>0)
             V=X.Cin(3,:);
             thc=2*X.m*Guide.GammaNi*X.lam0*X.g_uni;
             T=T+(thc^-2)*kron(V',V);
           end
           
        end
        
    end
end

