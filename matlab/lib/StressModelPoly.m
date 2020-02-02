%%StressModelSgm -  Segmented stress profile model
%
% Implements piece-wise polynomial model
% Depth (z) profile is defined by the nodes z_i and corresponding stress
% values sigma_i. 
% The model stress profile is defined by various interpolation models,
% ityp=
% 0) linear
% 1) cubic spline with linear ends (spline)
% 2) Piecewise cubic Hermite polynomial (pchip)
%% Fitted parameters:
% L(i) ... inner node positions along z [mm]
% sig(i) ... stress values at the nodes [0.001*E]
%% Other required parameters: 
% dsam ... sample thickness [mm]
% 
% Number of sig rows determines the number of segments, ns=size(sig,1)-1.
% Number of sig columns determines the number of stress components.
% L(i) are only the inner nodes, the first and last ones have positions 0
% and dsam, respectively.
% version: 1.1
% date: 26/2 2015
%

classdef StressModelPoly < FittingModel
    properties(Constant)
%% constants
%   
    end
    properties  


%% local variables
    % number of branches
    nbr;
    % number of segments
    nseg;
    % all nodes positions
    zs;
    % polynomial coefficients
    % rows: segments + stress components; columns: coefficients
    coef;
    % type of interpolation: linear(0), spline(1), pchip(2)
    ityp;
    % compliance matrix
    S;
    % instrumental parameters
    beta,zeta,dsam;
    % regularization
    regVal;
    end
    
    methods(Static)
    end
    
    methods 
        function G = StressModelPoly() 
            G=G@FittingModel();
        end                        
        
        % Regularization term (only for linear interpolation)
        % minimizes 2nd derivatives
        function R=RegValue(G)
            if (G.ityp==1); 
                sum=0;
                for k=1:G.nbr;
                  ib=(k-1)*G.nseg;
                  for i=1:G.nseg-1;                      
                      sum = sum+(G.coef(ib+i+1,3)-G.coef(ib+i,3))^2;
                  end;
                end;
                R=sum;
            else
                R=0;
            end; 
        end;
             
        
        % calculate polynom parameters according to the ityp value
        function calInterp(G)
            G.coef=zeros(G.nbr*G.nseg,4);
            %fprintf('calInterp, nbr=%d\n',G.nbr);
            for k=1:G.nbr;
              if (G.ityp==1);                
                %cs=interp1(G.zs(:,k),G.b(:,k),'spline','pp');
                cs=spline(G.zs(:,k),[0 G.b(:,k)' 0 ]);
              elseif(G.ityp==2);   
                cs=interp1(G.zs(:,k),G.b(:,k),'pchip','pp');
              else
                cs=interp1(G.zs(:,k),G.b(:,k),'linear','pp');
              end;
             % disp(cs.coefs);
              ib=(k-1)*G.nseg;
              np=min(size(cs.coefs,2),4);
              ic=5-np;
              %fprintf('calInterp, nseg=%d, ib=%d, ic=%d, np=%d\n',G.nseg,ib,ic,np);
              G.coef(ib+1:ib+G.nseg,ic:ic+np-1)=cs.coefs(1:G.nseg,1:np);
            end;
            G.regVal=G.RegValue();
        end  
        
        % set values for physical parameters       
        function setStressParam(G,L, sigma,fixL,fixSigma,dsam,ityp)
            G.ityp=ityp;
            G.setParam(L,sigma,fixL,fixSigma);
            % number of branches (stress components)
            G.nbr=size(sigma,2);
            % number of nodes
            nsig=size(sigma,1);
            % number of segments
            G.nseg=nsig-1;
            G.dsam=dsam;
            G.setVariables(G.vars);        
        end
              
        % set values for physical parameters       
        function setInstrParam(G,beta,zeta)
           G.beta=beta;
           G.zeta=zeta;
        end
        
         % set all variables as a single array      
        function setVariables(G,p)
            % parent method: set a,b arrays from p
           G.setVariables@FittingModel(p);
           L=G.a;
           % partitionling into segments of the interval (0,dsam)
           dimL=size(L,2);
           nL=size(L,1);
           G.zs=zeros(nL+2,G.nbr);
           % single set of nodes for all branches
           if (dimL==1);
             for k=1:G.nbr;
               G.zs(2:nL+1,k)=L(:);
               G.zs(nL+2,k)=G.dsam;
             end;
           % different nodes for each branch
           else
             for k=1:G.nbr;
               G.zs(2:nL+1,k)=L(:,k);
               G.zs(nL+2,k)=G.dsam;
             end;  
           end;
           % calculate interpolation coefficients.
           %fprintf('setVariables, nbr=%d, nL=%d\n',G.nbr,nL);
           %fprintf('setVariables, zs(%d,%d), b(%d,%d)\n',size(G.zs,1),size(G.zs,2),size(G.b,1),size(G.b,2));
           G.calInterp();
        end  
         
        % calculate stress values at given z depths
        function sig=getSigma(G,z)
            nd=size(z,1);
            sig=zeros(nd,G.nbr);
            for i=1:nd;
            for k=1:G.nbr;
            ib=(k-1)*G.nseg;
            for is=1:G.nseg;
                zz=z(i)-G.zs(is,k);                      
                if ((zz>=0) && (zz<=(G.zs(is+1,k)-G.zs(is,k))))
                    p=G.coef(ib+is,1:4);
                    sig(i,k)=p(1)*zz^3+p(2)*zz^2+p(3)*zz+p(4);
                end;
            end;  
            end;      
            end;
        end
        
        % calculate strain values at given z depths
        function e=getEps(G,z)
            nd=size(z,1);
            e=zeros(nd,G.nbr);
            sig=G.getSigma(z);
            for k=1:G.nbr
            for i=1:nd               
               e(i,k)=G.S(k,:)*sig(i,:)'*1e6;
            end;
            end;
        end     
        
        % calculate all components of smeared stress at a given depth
        function sigc=getSigmaConv(G,depth,beta, zeta)
            ns=G.nseg;
            A0=zeros(ns,G.nbr);
            A1=zeros(ns,G.nbr);
            A2=zeros(ns,G.nbr);
            A3=zeros(ns,G.nbr);
			% conversion to reduced length scale 
            x=beta*G.zs;
            u=beta*depth-zeta;
            d=beta*G.dsam;
            beta2=beta^2;
            beta3=beta^3;
            for k=1:G.nbr
              for i=1:ns
                ux=u-x(i,k);
                ux1=u-x(i+1,k);
                A0(i,k)=G0(ux,u,d)-G0(ux1,u,d);
                A1(i,k)=(G1(ux,u,ux,d)-G1(ux1,u,ux,d))/beta;
                if (G.ityp>0)
                  A2(i,k)=(G2(ux,u,ux,d)-G2(ux1,u,ux,d))/beta2;
                  A3(i,k)=(G3(ux,u,ux,d)-G3(ux1,u,ux,d))/beta3;
                end;
              end;
            end;
            sigc=zeros(G.nbr,1);
            for k=1:G.nbr
              ib=(k-1)*G.nseg;
              for i=1:G.nseg  
                p=G.coef(ib+i,1:4);
                sigc(k)=sigc(k)+p(1)*A3(i,k)+p(2)*A2(i,k)+p(3)*A1(i,k)+p(4)*A0(i,k);
              end;
            end;
        end
               
        % return smeared strain ic-th component and ix-th data point
        function y = modelFnc(G,ic,ix)
           mm=G.getSigmaConv(G.xval(ix),G.beta(ic),G.zeta(ic));
           y=G.S(ic,:)*mm*1e6;
           % fprintf('modelFnc x=%f, mm=%f, y=%f\n',G.xval(ix),mm,y);
        end
        
        % generate stress and strain curves with error bars for given
        % covariance matrix by MC simulation for nt trials
        function [stressErr strainErr] = genErrors(G, z, Cov, nt)
            [U,W,~] = svd(Cov);
            dev=(diag(W)).^0.5;
            p0=G.getFree();
            nz=numel(z);
            strainErr1=zeros(nz,G.nbr);
            strainErr2=zeros(nz,G.nbr);
            stressErr1=zeros(nz,G.nbr);
            stressErr2=zeros(nz,G.nbr);
            for it=1:nt
                rx=normrnd(0.0,dev);
                ry=U*rx;
                p=p0+ry';                
                G.setFree(p);
                tmp=G.getEps(z);
                tmp1=G.getSigma(z);
                strainErr1=strainErr1+tmp;
                strainErr2=strainErr2+tmp.^2; 
                stressErr1=stressErr1+tmp1;
                stressErr2=stressErr2+tmp1.^2;
            end;
            e2=strainErr2./nt-(strainErr1./nt).^2;
            strainErr=(e2.^0.5);
            e2=stressErr2./nt-(stressErr1./nt).^2;
            stressErr=(e2.^0.5);   
        end
        
        
        
         function res = validated(G)
            res=G.validated@FittingModel();
            a=(size(G.S,1) == G.nset);
            if (~ a)
                fprintf('Number of data sets, %d <> %d , S rows\n',G.nset,size(G.S,1) );
                res = res && a;
            end;
            a=(size(G.S,2) == G.nbr);
            if (~ a)
                fprintf('Number of stress components, %d <> %d , S columns\n',G.nbr,size(G.S,2) );
                res = res && a;
            end;
            a = (size(G.zs,1) == G.nseg+1);
            if (~ a)
                fprintf('Numbers of segments for L and do not agree: rows(L)+2=%d, rows(sigma)=%d\n',size(G.zs,1),G.nseg+1);
                res = res && a;
            end;
            a = (numel(G.beta) == G.nset && numel(G.zeta) == G.nset );
            if (~ a)
                fprintf('Dimensions of beta, zeta must agree with the number of data sets: nBeta=%d, nZeta=%d, nset=%d\n',...
                    numel(G.beta),numel(G.zeta),G.nset);
                res = res && a;
            end;
        end;
        
        function nam = getNames(G)
          dimL=size(G.a,2);
          nn=0;
          if (dimL>1);
            nc=2*G.nseg*G.nbr;
            nam=cell(1,nc);  
            for k=1:G.nbr
                for i=1:G.nseg-1  
                    nn=nn+1;
                    nam{nn}=sprintf('L(%d,%d)',i,k);             
                end
            end
          else
            nc=G.nseg-1+(G.nseg+1)*G.nbr;
            nam=cell(1,nc);            
            for i=1:G.nseg-1  
                nn=nn+1;
                nam{nn}=sprintf('L(%d)',i);             
            end
          end;
          for k=1:G.nbr
            for i=1:G.nseg+1
                nn=nn+1; 
                nam{nn}=sprintf('sigma(%d,%d)',i,k);             
            end
          end
        end;
    end
    
end

