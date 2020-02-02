%%Pulse
%Matrix model of a pulsed source 

classdef Pulse < Comp
    properties(Constant)     
    end
    properties  
        tau;  % time width [us]       
    end

    methods              
        function X=Pulse(distance,tau,shape)
           X = X@Comp(distance);
           if (nargin > 2) 
             X.tau=tau*shape;
           else
             X.tau=tau*Comp.g_norm;    
           end  
           X.id='Pulse';
        end
    
        % Transmission matrix
        function [T] =getT(X)
           V=X.Cin(5,:);
           T=getT@Comp(X);
           if (X.tau>0)
             T=T+(X.tau^-2)*kron(V',V);
           end
        end
        
    end
end

