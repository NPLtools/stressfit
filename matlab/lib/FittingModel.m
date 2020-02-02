%%FittingModel -  basic class with a fitting model
%
% Implements the basic functions for a primitive model:
% y = a*x+b 

classdef FittingModel < matlab.mixin.Heterogeneous & handle
    properties(Constant)
%% constants
% degree
deg=pi/180;
%   
    end
    properties  


%% local variables
    a,b;    % physical parameters
    vars;   % array of all variable parameters
    fixed;  % flags 0|1 for free|fixed variables
    
    ndat;   % partitioning of x-values into multiple data sets
    nset;   % number of data sets
    nx;     % total number of x values (all data sets)
    xval;   % array with all x values
    xsgn;   % array of signs for the x-scale (one per data set)
            % ... xsgn allows to switch x-scale direction of the model; this feature
            % must be implemented in each descendant of the fitting model in modelFnc. 
    
    end
    
    methods(Static)
        function v = param2vars(a,b)            
           nax=size(a,1);
           nay=size(a,2);
           nbx=size(b,1);
           nby=size(b,2);
           c=reshape(a,1,nax*nay);
           if (nbx*nby>0)
               d=reshape(b,1,nbx*nby);
               v=cat(2,c,d);
           else
               v=c;
           end;
        end
                        
        function printString(txt, fid)
             fprintf('%s\n',txt);
             if (nargin==2);
                fprintf(fid,'%s\n',txt);   
             end;         
        end;
    end
    
    methods 
        function G = FittingModel() 
        end        
        
        function [a b] = vars2param(G,p)
           G.vars=p; 
           nax=size(G.a,1);
           nay=size(G.a,2);
           nbx=size(G.b,1);
           nby=size(G.b,2);
           n1=1;n2=nax*nay;
           a=reshape(p(n1:n2),nax,nay);
           n1=n2+1;n2=n2+nbx*nby;
           if (nbx*nby>0)
               b=reshape(p(n1:n2),nbx,nby);
           else
               b=[];
           end;
        end
        
        
        % set values for physical parameters
        % variant with only one parameter set
        function setParam1(G,a,afix)
           G.a = a;
           G.b = [];      
           G.vars=FittingModel.param2vars(a,[]);
           G.fixed=FittingModel.param2vars(afix,[]);
        end
        
        % set values for physical parameters       
        function setParam(G,a,b,afix,bfix)
           G.a = a;
           G.b = b;      
           G.vars=FittingModel.param2vars(a,b);
           G.fixed=FittingModel.param2vars(afix,bfix);
        end
        
        % get flags for fixed parameters       
        function f = getFixed(G)
           f=G.fixed;
        end   
        
        % get all variables as a single array      
        function v = getVariables(G)
           nax=size(G.a,1);
           nay=size(G.a,2);
           nbx=size(G.b,1);
           nby=size(G.b,2);
           c=reshape(G.a,1,nax*nay);
           if (nbx*nby>0)
               d=reshape(G.b,1,nbx*nby);
               v=cat(2,c,d);
           else
               v=c;
           end;
        end  
        
         % set all variables as a single array      
         function setVariables(G,p)           
           [G.a G.b] = G.vars2param(p);
         end   
        
         % select only non-fixed variables out of vars 
         function f = getFree(G)
            n=numel(G.vars);
            j=0;
            x=zeros(1,n);
            for i=1:n
                if (G.fixed(i)==0)
                    j=j+1;
                    x(j)=G.vars(i);
                end
            end
            f=x(1:j);
         end
         
         % set only non-fixed variables
         function  setFree(G,p)
            j=0;
            n=numel(G.vars);
            x=G.vars;
            for i=1:n
                if (G.fixed(i)==0)
                    j=j+1;
                    x(i)=p(j);
                end
            end
            G.setVariables(x);
         end
            
        % define x-values
        % npart is the partitioning of the x-array into multiple data sets
        % xsgn is an optional row-array with signs for x-scale of each dataset
        function setX(G,npart,x,xsgn)
            G.ndat=npart;
            G.nset=numel(npart);
            G.nx=min(size(x,1),sum(npart));
            G.xval=x;
            if (nargin<4)
                G.xsgn=ones(1,G.nset);
            else
                G.xsgn=xsgn;
            end;
        end  
        
         % return values of the fitting function for given x array
         % npart is the partitioning of the x-areray into multiple data
         % sets
        function y = fitFnc(G)  
            y=zeros(G.nx,1);
            n1=1;
            % loop through data sets
            for is=1:G.nset   
                n2=min(n1-1+G.ndat(is),G.nx);
                % loop through data points
                for ix=n1:n2 
                    y(ix)=G.modelFnc(is,ix);
                    %fprintf('fitFnc is=%d, ix=%d, y=%f\n',is,ix,y(ix));
                end;
                n1=n1+G.ndat(is);
            end;
            %display(y);
        end
        
        % return model value for is-th data set and ix-th data point
        function y = modelFnc(G,is,ix)
           y = G.xsgn(is)*G.a(is)*G.xval(ix)+G.b(is); 
        end
        
        function res = validated(G)
            res=(numel(G.a)+numel(G.b)==numel(G.vars));
            if (~res)
              fprintf('Model was not properly initialized\n');
              disp(G.a);
              if (numel(G.b)>0)
                disp(G.b);
              end;
            end
        end
        
        function nam = getNames(G)
          nam=cell(1,2*G.nset); 
          nn=0;
          for i=1:G.nset      
             for j=1:numel(G.a) 
                nn=nn+1; 
                nam{nn}=sprintf('a(%d,%d)',j,i);
             end;
          end
          if (numel(G.b)>0)
            for i=1:G.nset      
                for j=1:numel(G.b) 
                    nn=nn+1; 
                    nam{nn}=sprintf('b(%d,%d)',j,i);
                end;
            end
          end;
        end;
        

        
        function printReport(G, fid, header,names,values,errors, chi2)
            if (~isempty(header))
                fprintf(fid,'-------------------------------------\n');
                fprintf(fid,'%s\n',header);
            end                    
            if (chi2>0)
                fprintf(fid,'chi2=%f\n',chi2);
            end
            n= size(G.fixed,2);
            nf=0;
            nv=0;
            for i=1:n
                if (G.fixed(i))
                    nf=nf+1;
                else
                    nv=nv+1;
                end
            end
            if (nv>0) 
                fprintf(fid,'Fitted variables:\n');
                for i=1:n
                    if (~ G.fixed(i))
                        fprintf(fid,'%s\t%f\t%f\n',names{i},values(i),errors(i));
                    end;
                end;
            end;
            if (nf>0)
                fprintf(fid,'Fixed variables:\n');
                for i=1:n
                    if (G.fixed(i))
                        fprintf(fid,'%s\t%f\n',names{i},values(i));
                    end;
                end;
            end;            
        end;
        
        function reportResults(G,filename, log2scr, header,values,errors, chi2)        
            what=zeros(1,2);
            names=G.getNames();
            if (log2scr)
                what(1)=1;
            end
            if (~isempty(filename))
                what(2)=1;
            end
            % ------ start loop with log output ------
            for iw=1:2
                if (what(iw)>0)
                    if (iw==1)
                        G.printReport(1, header,names,values,errors, chi2);
                    else
                        fid=fopen(filename,'at'); 
                        G.printReport(fid, header,names,values,errors, chi2);
                        fclose(fid);
                    end
                end;
            end;   
        end
    end
end

