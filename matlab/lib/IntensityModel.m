%%IntensityModel -  basic class with a fitting model
%
% model: $$ y = y_0 + I_0\,T\,Intensity(\beta (z-z_0)-\zeta,\beta t) $$
% where $$ Intensity(u,\tau) = F_0(u,\tau) + b\,F_1(u,\tau) + c\,F_2(u,\tau) $$ 
% and $$ T = exp(-\mu \, t \, c_g)  $$ 
% data:   (x,y,error) table
%
% properties:
% ag  .. geometry factors for absorption
% cg  .. geometry factors for through-thickness transmission (=0 for reflection geometry)
%
% param:  [y0, I0, z0, beta, mu, t, b, c]
% y0   ... background
% I0   ... intensity
% z0   ... surface position [mm]
% beta ... gauge volume size parameter [mm-1]
% mu   ... attenuation coefficient [cm-1]
% t    ... sample thickness [mm] 
% b    ... linear inhomogeneity coeficient
% c    ... quadratic inhomogeneity coefficient
%
% call setParam1(a,afix) to set variables and fix flags
% a,afix are arrays of variables and flags values in the above order  
%
% provide a(i) as array columns 
% use multiple columns of a() for fitting multiple data sets simultaneously
%
% version: 1.3
% date: 2/8/2016
%


classdef IntensityModel < FittingModel
    properties(Constant)  
    end
    
% local variables
    properties  
      ag;      % geometry factors for absorption
      cg;      % geometry factor for through-thickness transmission (=0 for reflection geometry)
    end
    
    methods(Static)
    end
    
    methods 
        function G = IntensityModel() 
            G=G@FittingModel();
        end        
        
        function [y0, I0, z0, beta, mu, t, b, c] = getResult(G,is)          
           y0=G.a(1,is);
           I0=G.a(2,is);
           z0=G.a(3,is);
           beta=G.a(4,is);
           mu=G.a(5,is);
           t=G.a(6,is);
           b=G.a(7,is);
           c=G.a(8,is);
        end
        
        % return model value for is-th data set and ix-th data point
        function y = modelFnc(G,is,ix)
          [y0, I0, z0, beta, mu, t, b, c] = G.getResult(is);
           zeta=G.ag(is)*mu/2/beta;
           tau=t*beta;
           T=exp(-mu*t*G.cg(is)); % in transmission geometry, add constant transmission factor
           u=G.xsgn(is)*beta*(G.xval(ix)-z0)-zeta;                
           y = y0+I0*T*PSModel.gvol(u,tau, zeta, b, c); 
           % y = y0+I0*T*P0(u,zeta)*(F0(u,tau)+b*F1(u,tau)+c*F2(u,tau)); 
        end

        function nam = getNames(G)
          nc=6*G.nset;
          nam=cell(1,nc);   
          nn=0;
          vnam={'y0','I0','z0','beta','mu','thickness','b','c'};
          for i=1:G.nset  
             for j=1:8
                nn=nn+1;  
                if (G.nset>1)
                    nam{nn}=sprintf('%s(%d)',vnam{j},i); 
                else
                    nam{nn}=sprintf('%s',vnam{j});
                end;
             end;
          end;
        end;
        
        function printReport(G, fid, header,names,values,errors, chi2)
            G.printReport@FittingModel(fid, header,names,values,errors, chi2);
            [va, ~] = G.vars2param(values);
            [ea, ~] = G.vars2param(errors);
            betaZ=va(4,:);
            aZ=va(5,:);
            betaZerr=ea(4,:);
            aZerr=ea(5,:);
            ns=size(aZerr,2);
            for i=1:ns               
                zeta=aZ(i)*G.ag(i)/2/betaZ(i);
                errzeta=1/2/betaZ(i)*sqrt(aZerr(i)^2+(aZ(i)*betaZerr(i)/betaZ(i))^2);
                fprintf(fid,'%s\t%f\t%f\n','zeta',zeta,errzeta);
                %fprintf(fid,'Absorption parameter(%d):\n',i);
                %fprintf(fid,'\t%s = %f +- %f\n','zeta = a/2/beta',zeta,errzeta);
               % fprintf(fid,'\t%s = %f +- %f\n','a',aZ(i)*G.ag(i),aZerr(i)*G.ag(i));
                %fprintf(fid,'Gauge volume width(%d):\n',i);
                fwhm=sqrt(4*log(2))/betaZ(i);
                fprintf(fid,'%s\t%f\t%f\n','fwhm[mm]',fwhm,betaZerr(i)/betaZ(i)*fwhm);
                %fprintf(fid,'\t%s = %f +- %f\n','fwhm',fwhm,betaZerr(i)/betaZ(i)*fwhm);
            end;         
            %n= size(G.fixed,2);
            %fprintf(fid,'Table:\n');
            %for i=1:n
            %    fprintf(fid,'%s\t%f\t%f\n',names{i},values(i),errors(i));
            %end;
        end
    end
end

