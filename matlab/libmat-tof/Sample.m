%% Sample
% Matrix model of a sample
% 

classdef Sample < Comp
    properties(Constant)
    end
    properties  
        % input parameters
        dhkl;          % d-spacing [A]
        gamma;         % surface angle [rad]
      %  mu;            % absorption coefficients, linear wavelength dependence
		thickness;     % layer thickness in mm
        % other calculated parameters
        theta;         % Bragg angle [rad]
        a;             % attenuation coefficients (linear depth dependence)
      %  mu_tot;        % absorption coefficient for given wavelength
    end
    methods      
        function X =Sample(dhkl,gamma,thickness)  
            X = X@Comp();
            X.dhkl=dhkl;
            X.gamma=gamma;    
           % X.mu=mu;
		    X.thickness=thickness;            
            X.id='Sample';
        end
        
        % recalculate dependent data after change of input parameters
        function initialize(X,C,lam0)  
            X.initialize@Comp(C,lam0);
            X.theta=asin(lam0/2/X.dhkl);
          %  X.mu_tot=X.mu(1)+X.mu(2)*X.lam0+X.mu(3)*X.lam0^2;           
            % attenuation - depth dependence: a_tot=a(1)*d + a(2)
            cos1=cos(X.gamma);
            cos2=cos(2*X.theta-X.gamma);
            X.a(1)=1/cos2-1/cos1;
            X.a(2)=0.5*(sign(cos2)/cos2 + sign(cos1)/cos1 - X.a(1));            
        end 
        
        % transform to the component local coordinates from the output
        % coordinates of the preceding one
        function [C] =getCin(X)            
           C=CT(0.0,-X.gamma, 0.0);           
        end
        
         % Get transformation to output coordinates.
         % Includes scattering and rotation by 2*theta-gamma
        function [Cout] =getCout(X) 
            S=[eye(2,2), zeros(2,3); zeros(2,2), Sigma(X.theta), zeros(2,1); zeros(1,4), 1.0 ] ;
            Cout=CT(0.0, 2*X.theta-X.gamma, 0.0)*S;
        end
        
    end
end
