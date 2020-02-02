%% Extinction
% Describes beam extinction in a polycrystal
% Constructor parameters:
% sigma = attenuation coefficients [cm^-1]: 
%      sigma(1) = free atom scattering cross section in 1/cm
%      sigma(2) = incoherent scattering cross section in 1/cm
%      sigma(3) = absorption + single phonon scattering cross section in 1/cm/A
%      sigma(4) = quadratic term in 1/cm/A^2 (only used without table)
% Optional:
% table = lookup table with lambda,sigma_tot columns [A, cm^-1]
% 
% Outside the table range, the extrapolation is used:
% sigma(1) for short wavelengths
% sigma(2) + sigma(3)*lambda for long wavelengths
% Without the table, the long-wavelength approx. is used
%

classdef Extinction < handle 
    properties(Constant)
    end
    properties  
       % input parameters
        sig;   % attenuation coefficients
        tab;   % table lambda,sigma_tot in {A, cm^-1]
        nl;    % table number of rows and range
        lmin;
        lmax;
    end
    
    methods
        function X =Extinction(sigma,table)  
            X.sig=sigma;              
            if (nargin>1) 
              X.tab=table; 
              X.nl=size(table,1);
              X.lmin=table(1,1);
              X.lmax=table(X.nl,1);
            else 
              X.tab=0;
              X.nl=0;
              X.lmin=0;
              X.lmax=0;
            end;
        end 
        % return attenuation cross section [1/cm] for given wavelength
        function [mu] =getMu(X,lambda)  
            if (X.nl<1)
               mu=X.sig(2)+X.sig(3)*lambda+X.sig(4)*lambda^2;
            elseif (lambda<X.lmin) 
               mu=X.sig(1);
            elseif (lambda>X.lmax) 
               mu=X.sig(2)+X.sig(3)*lambda;
            else
               mu=interp1q(X.tab(:,1),X.tab(:,2),lambda);
            end                
        end
        % mu averaged over given spectrum spec=table(lambda,flux)
        function [mu] =getMuAve(X,spec)  
            nb=size(spec,1);
            mu=zeros(nb);
            sum=0;
            for j=1:nb
                sum=sum+spec(j,2);
                mu=X.getMu(spec(j,1))*spec(j,2);
            end
            mu=mu/sum;               
        end
    end
    
end

