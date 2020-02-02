%% Comp
% Matrix model of an abstract component, ToF version
% Phase space variables are indexed in following order:
% x, z, angle, d_lambda/lambda, time
% This component represent just a transparent plane
% Actual components like slits or collimators should be derived from Comp
% 

classdef Comp < matlab.mixin.Heterogeneous & handle
    properties(Constant)
% shape factors (e.g.sqrt(1/6) for a rectangular slit, sqrt(1/3) for radial collimator, sqrt(1/8) for circular slit)              
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
    properties  
        L;              % distance from previous [mm]  
        lam0;           % nominal wavelength
        Cin;            %  transport matrix to local coordinates
        Cpre;           %  transport matrix to local coordinates        
        id='Comp';      % name
    end

    methods      
        function X=Comp(distance)
            X.Cin=eye(5,5);
            X.Cpre=eye(5,5);
            if (nargin>0) 
              X.L=distance; 
            else  
              X.L=0.0;
            end
        end        
       
        % recalculate dependent data after change of input parameters
        function initialize(X,C,lam0)  
            X.Cpre=C;
            X.lam0=lam0;
            X.Cin=X.getCin();
        end 
        
        % transform to the component local coordinates from the output
        % coordinates of the preceding one
        function [C] =getCin(X)            
           C=CT( X.L, 0.0, X.lam0 )*X.Cpre;
        end
        
        % Transmission matrix
        function [T] =getT(X) %#ok<MANU>
           T=zeros(5,5);
        end
        
         % Get transformation to output coordinates. It may include
         % scattering for some components.
        function [Cout] =getCout(X) 
           Cout=X.Cin;
        end
        
        function display(X) 
            fprintf('%s.Cin=\n',X.id);
            disp(X.Cin);           
        end;
        
    end
end

