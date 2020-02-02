%%StressModelWell -  Well stress profile model
%
% Implements stress profile as a well function.
% Free parameters are:
% L(1) ... depth of the stressed layer
% L(2) ... width of the relaxation layer
% sig  ... stress inside the stressed layer 
% Number of sig columns determines the number of components.
% version: 1.2
% date: 14/4 2014
%

classdef StressModelWell < StressModelSgm
    properties(Constant)
%% constants
%   
    end
    properties  


%% local variables
    

    end
    
    methods(Static)
    end
    
    methods 
        function G = StressModelWell() 
            G=G@StressModelSgm();
        end                        
                      
        
        % set values for physical parameters       
        function setStressParam(G,L, sigma,fixL,fixSigma,dsam)
            G.setStressParam@StressModelSgm(L, sigma,fixL,fixSigma,dsam);
            G.nseg=2;
        end
        
         % set all variables as a single array      
        function setVariables(G,p)
           G.setVariables@FittingModel(p);
           G.sigma=G.b;
           L=G.a;
           % partitionling into segments of the interval (0,dsam)
           G.zs=zeros(3,1);
           G.zs(2)=L(1);
           G.zs(3)=L(1)+L(2);
           % calculate gradients
           G.alpha=zeros(G.nbr,1);
           for k=1:G.nbr
                if (L(2)>0)               
                    G.alpha(k)=-G.sigma(k)/L(2);
                else
                    G.alpha(k)=0;
                end;
           end;
        end  
         
        % calculate stress values at given z depths
        function sig=getSigma(G,z)
            nd=size(z,1);
            sig=zeros(nd,G.nbr);
            for i=1:nd
            for k=1:G.nbr
                zz=z(i);
                if ((zz>=0) && (zz<=G.zs(2)))
                    sig(i,k)=G.sigma(k);
                elseif ((zz>G.zs(2)) && (zz<=G.zs(3)))
                    sig(i,k)=G.sigma(k)+G.alpha(k)*(zz-G.zs(2));
                end;  
            end;      
            end;
        end
               
        % calculate values of smeared stress at a given depth, h
        % assumed reduced depth, h = beta*depth-zeta
        % zn=zs*beta, d = dsam*beta
        function sigc=getSigmaConv(G,depth,beta, zeta)
            GG0=zeros(3,1);
            GG1=zeros(3,1);
            zn=beta*G.zs;
            h=beta*depth-zeta;
            d=beta*G.dsam;
            for i=1:3
                GG0(i)=G0(h-zn(i),h,d);
                GG1(i)=G1(h-zn(i),h,zeta,d);
            end;
            sigc=zeros(G.nbr,1);
            for k=1:G.nbr
                a=GG0(2)-GG0(3);
                b=(GG1(2)-GG1(3))/beta;
                sigc(k)=sigc(k)+(GG0(1)-GG0(2))*G.sigma(k);
                sigc(k)=sigc(k)+a*(G.sigma(k)+G.alpha(k)*(depth-G.zs(2))) -  b*G.alpha(k);
            end;  
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
                disp(G.S);
                res = res && a;
            end;            
        end;
        
        function nam = getNames(G)
          nc=2+G.nbr;
          nam=cell(1,nc);            
          for i=1:2            
             nam{i}=sprintf('L%d',i);             
          end
          nn=2;
          for k=1:G.nbr
             nn=nn+1; 
             nam{nn}=sprintf('sigma(%d)',k);             
          end
        end;
    end
    
end

