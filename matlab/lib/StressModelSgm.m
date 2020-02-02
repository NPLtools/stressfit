%%StressModelSgm -  Segmented stress profile model
%
% Implements Segmented stress profile model
% Depth (z) profile is defined by straight lines on given intervals.
% It is therefore defined by:
% L(i) ... inner node positions along z [mm]
% sig(i) ... stress values at the nodes [0.001*E]
% dsam ... sample thickness [mm]
% Number of sig rows determines the number of segments, ns=size(sig,1)-1.
% Number of sig columns determines the number of components.
% L(i) are only the inner nodes, the first and last ones have positions 0
% and dsam, respectively.
% version: 1.2
% date: 7/4 2014
%

classdef StressModelSgm < FittingModel
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
    % gradients
    alpha;
    % nodes values
    sigma;
    % compliance matrix
    S;
    % instrumental parameters
    beta,zeta,dsam;
    

    end
    
    methods(Static)
    end
    
    methods 
        function G = StressModelSgm() 
            G=G@FittingModel();
        end                        
        
        % set values for physical parameters       
        function setStressParam(G,L, sigma,fixL,fixSigma,dsam)
            G.setParam(L, sigma,fixL,fixSigma);
            % number of branches
            G.nbr=size(sigma,2);
            % number of nodes
            nsig=size(sigma,1);
            % number of segments
            G.nseg=nsig-1;
            G.sigma=sigma;
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
           G.setVariables@FittingModel(p);
           G.sigma=G.b;
           L=G.a;
           % partitionling into segments of the interval (0,dsam)
           nL=numel(L);
           G.zs=zeros(nL+2,1);
           G.zs(2:nL+1)=L(:);
           G.zs(nL+2)=G.dsam;
           % calculate gradients
           G.alpha=zeros(G.nseg,G.nbr);
           for k=1:G.nbr
           for i=1:G.nseg
              G.alpha(i,k)=(G.sigma(i+1,k)-G.sigma(i,k))/(G.zs(i+1)-G.zs(i));
           end;
           end;   
        end  
         
        % calculate stress values at given z depths
        function sig=getSigma(G,z)
            nd=size(z,1);
            sig=zeros(nd,G.nbr);
            for i=1:nd
            for k=1:G.nbr
            for is=1:G.nseg
                zz=z(i);
                if ((zz>=G.zs(is)) && (zz<=G.zs(is+1)))
                    sig(i,k)=G.sigma(is,k)+G.alpha(is,k)*(zz-G.zs(is));
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
        
        % calculate values of smeared stress at a given depth, h
        % assumed reduced depth, h = beta*depth-zeta
        % zn=zs*beta, d = dsam*beta
        function sigc=getSigmaConv(G,depth,beta, zeta)
            nsig=G.nseg+1;
            GG0=zeros(nsig,1);
            GG1=zeros(nsig,1);
            zn=beta*G.zs;
            h=beta*depth-zeta;
            d=beta*G.dsam;
            for i=1:nsig
                GG0(i)=G0(h-zn(i),h,d);
                GG1(i)=G1obs(h-zn(i),h,zeta,d);
            end;
            sigc=zeros(G.nbr,1);
            for k=1:G.nbr
            for i=1:G.nseg
                a=GG0(i)-GG0(i+1);
                b=(GG1(i)-GG1(i+1))/beta;
                sigc(k)=sigc(k)+a*(G.sigma(i,k)+G.alpha(i,k)*(depth-G.zs(i))) - b*G.alpha(i,k);
            end;
            end;  
        end
               
        % return model value for is-th data set and ix-th data point
        function y = modelFnc(G,is,ix)
           mm=G.getSigmaConv(G.xval(ix),G.beta(is),G.zeta(is));
           y=G.S(is,:)*mm*1e6;
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
                fprintf('Dimensions of L and sigma arrays do not agree: rows(L)+2=%d, rows(sigma)=%d\n',size(G.zs,1),G.nseg+1);
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
          nc=G.nseg-1+(G.nseg+1)*G.nbr;
          nam=cell(1,nc);            
          for i=1:G.nseg-1              
             nam{i}=sprintf('L(%d)',i);             
          end
          nn=G.nseg-1;
          for k=1:G.nbr
          for i=1:G.nseg+1
             nn=nn+1; 
             nam{nn}=sprintf('sigma(%d,%d)',i,k);             
          end
          end
        end;
    end
    
end

