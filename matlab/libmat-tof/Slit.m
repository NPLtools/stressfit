%%Slit
%Matrix model of a slit or radial collimator (distance=0)
% 

classdef Slit < Comp
    properties(Constant)     
    end
    properties  
        w;  % width [mm] 
    end

    methods              
        function X=Slit(distance,width)
           X = X@Comp(distance);
           if (distance == 0) 
             X.w=width*Comp.g_tri;
           else
             X.w=width*Comp.g_uni;     
           end  
           X.id='Slit';
        end        
    
        % Transmission matrix
        function [T] =getT(X)
           V=X.Cin(1,:);
           T=getT@Comp(X);
           if (X.w>0)
             T=T+(X.w^-2)*kron(V',V);
           end
        end
        
    end
end

