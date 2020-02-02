%%PKS_bent_thin
% Analytical model for neutron diffraction peak shifts due to the surface
% effect
% J. Saroun, Nuclear Physics Institute, Rez, saroun@ujf.cas.cz
% submitted to Journal of Applied Crystallography,  25 Jan 2013
%------------------------------------------------------------------------
% Calculates peak shift and optimum curvature for a thin bent
% monochromator. 
%
% see Section 4.1
%


classdef PKS_bent_thin    < handle
    
    properties      
    % dimensions
    dM;     % monochromator thickness [mm]
    w0;     % slit 0: divergence slit width [mm]
    w1;     % slit 1: input slit width [mm]
    w2;     % slit 2: output slit width [mm]
    % distances
    LM;     % distance slit0 .. monochromator [mm]
    L0;    % distance slit0 .. slit 1 [mm]
    L1;	% distance slit1 .. sample axis [mm]
    L2;	% distance sample .. slit 2 [mm]
    LD;     % distance slit 2 .. detector [mm]
    % scattering angles
    thetaM;     % take-off angle, monochromator [deg]
    gammaM;     % surface orientation angle
    % other
    chi;        % monochromator cutting angle (0 for symmetric reflection)
    nu;         % Poisson elastic constant
    mu;         % absorption coefficient [1/mm]
    end
    
    methods(Static)
        function [res] = g_ref(theta)
        % sample rotation angle - symmetric reflection
           res=theta+pi/2; 
        end        
        function [res] = g_tra(theta)
        % sample rotation angle - symmetric transmission
           res=theta+pi; 
        end
    end
    
    methods        
        % Set input parameters: instrument components for this special case:
        % monochromator, slit 0, slit 1, slit 2 and detector, represented
        % by respectivce classes (CrystalBent, Slit)        
        function setParam(X,spec, mu)
        % Instrument setup data
            % distances
            mon=spec(1);
            s0=spec(2);
            s1=spec(3);
            s2=spec(4);
            det=spec(5);
            X.LM=-mon.L;        % distance slit0 .. monochromator [mm]
            X.L0=-s0.L;         % distance slit0 .. slit 1 [mm]
            X.L1=-s1.L;     	% distance slit1 .. sample axis [mm]
            X.L2=s2.L;      	% distance sample .. slit 2 [mm]
            X.LD=det.L;     	% distance slit 2 .. detector [mm]
            % dimensions
            X.dM=mon.w;         % monochromator thickness [mm]
            X.w0=s0.w;          % slit 0: divergence slit width [mm]
            X.w1=s1.w;          % slit 1: input slit width [mm]
            X.w2=s2.w;          % slit 2: output slit width [mm]                       
            % other
            X.thetaM=mon.thetaB;    % take-off angle, monochromator [deg] 
            X.chi=mon.chi;          % monochromator cutting angle (0 for symmetric reflection)
            X.nu=mon.nu;            % Poisson elastic constant
            X.gammaM=X.gamma_M();   % surface orientation angle
            X.mu=mu;
        end
        
        function [res] = a(X,theta)
            res=-tan(theta)/tan(X.thetaM);
        end
        
        function [res] = gamma_M(X)
        % angle between monochromator surface normal and axis of the
        % reflected beam)
            res=sign(X.thetaM)*pi/2 -(X.thetaM+X.chi);
        end                        

        function [res] = x2dd(X,theta)
        % Convert peak shift at the detector to delta_d/d
            res=-1/(2*X.LD*tan(theta));
        end  
        
        function [a] = getAbs(X,theta,gamma)
            a=X.mu*(1/cos(2*theta-gamma) - 1/cos(gamma));          
        end
        
        %% calculate coefficients P,D of the quadratic forms  
        % in the formula for peak shift coefficients, eqs. A28-A30.
        function [P,D] = shift_coeff(X,theta,gamma)
            b=-X.a(theta);
            CG=cos(gamma);
            CT=cos(2*theta-gamma); 
            cgm=cos(X.gammaM);
            LM0=X.LM+X.L0;
            L10=X.L1+X.L0;
            Lam0=(1+2*b)*X.L2*CG + X.L1*CT;
            Lam1=Lam0 + X.L0*CT;            
            Z=2*b/cgm;            
            P=zeros(1,3);
            D=zeros(1,3);
           % display(tmp);
            %tmp=(Lam0*X.w0^2 + Lam1*X.w1^2)*(1+2*b)*X.LD - (X.L0*X.w2)^2*CG;
            
            P(1)=(Lam0*X.w0^2 + Lam1*X.w1^2)*(1+2*b)*X.LD - (X.L0*X.w2)^2*CG;
            P(2)=-(LM0*(2*Lam0-X.L1*CT)*X.w0^2 + X.LM*(2*Lam1-L10*CT)*X.w1^2)*X.LD*Z;
            P(3)=Z^2*((X.LM*X.w1)^2 + (LM0*X.w0)^2)*X.LD*X.L2*CG;            
            D(1)=(Lam0*X.w0)^2 + (Lam1*X.w1)^2 + (X.L0*CG*X.w2)^2;
            D(2)=-( LM0*Lam0*X.w0^2 + X.LM*Lam1*X.w1^2 )*2*Z*X.L2*CG;
            D(3)=((X.LM*X.w1)^2 + (LM0*X.w0)^2)*(Z*X.L2*CG)^2;                                              
        end
        
        function [B,beta] = shift_param(X,theta,gamma,ro)
        % Calculate peak shift parameters used in A26. 
        % B = Delta_x/beta, where Delta_x and beta arr calculated according to
        % A27-A29. 
            [P,D]=X.shift_coeff(theta,gamma); % get polynomial coefficients from A30.
            SD=D(1)+D(2)*ro+D(3)*ro^2;
            SP=P(1)+P(2)*ro+P(3)*ro^2;
            Z=sqrt(abs(SD)/6); % = sqrt(Z(ro)) in A29
            S2T=sin(2*theta);
            B=SP/(6*X.L0*Z);
            beta=X.L0*S2T/Z;
        end
        
        function [res] = shift(X,delta,theta,gamma,ro)
        % calculate peak shift for given delta, theta and gamma, eq. A26
        % delta ... surface position
        % theta ... Bragg angle for sample
        % gamma ... surface orientation (theta+pi/2 for symmetric reflection)
        % ro    ... monochromator curvature [m-1]
            n=length(delta);
            [B,beta]=X.shift_param(theta,gamma,ro/1000);            
            res=zeros(n,1);
            zeta=X.getAbs(theta,gamma)/beta;
            for i=1:n
              res(i)=B*f(beta*delta(i),zeta);
            end
        end
        
        function [res] = rho_opt(X,theta,gamma)
        % calculate optimum monochromator curvature [m^-1] for peak shift
        % general formula
        % returns two solutions, ordered by magnitude
        % S are coefficients returned by |shift_coeff(theta,gamma)|
            eps=1e-10;
            res=zeros(1,2);
            [S,dum]=X.shift_coeff(theta,gamma); %#ok<NASGU>
            if (abs(S(3))<eps)
                res(:)=-S(1)/S(2)*1000;
                return;
            end;
            D=S(2)^2 - 4*S(3)*S(1);
            if (D>eps)
                DD=sqrt(D);
                a1=(-S(2)+DD)/(2*S(3))*1000;
                a2=(-S(2)-DD)/(2*S(3))*1000;
                if (abs(a1)>abs(a2))
                    res=[a2,a1];
                else
                    res=[a1,a2];
                end;
            end;
        end
        
        function [res] = rho_res(X,theta)
        % calculate optimum monochromator curvature [m^-1] for resolution
            cgm=cos(X.gammaM);
            a=X.a(theta);
            res=cgm*(2*a-1)/((X.LM+X.L0)*2*a)*1000;            
        end
        
        function [res] = rho_anal(X,theta,gamma)
        % calculate optimum monochromator curvature [m^-1] for peak shift
        % analytical formula for L2=0 (secondary soller collimator)
        % returns two solutions
        LM0=X.LM+X.L0;
        L10=X.L1+X.L0;
        K=-cos(gamma)/cos(2*theta-gamma);        
        a=X.a(theta);        
        C1=X.L0*L10*X.w1^2 - K*LM0*(X.L0*X.w2)^2/(X.LD*(2*a-1));
        C2=X.L1*LM0*X.w0^2 + L10*X.LM*X.w1^2;
        rho1=rho_res(theta);
        res=rho1*(1+C1/C2);        
        end
        
    end
    
end

