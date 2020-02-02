function [ afit,bfit,aerr,berr,fit,Cov ] = fit_stress_profile(npart,data,S,beta,zeta,dsam,a0,b0,fix,ni)
% Fit the model of stress profile:
% $$ \sigma_{i}(z) = (a_{i,1} + a_{i,2}z a_{i,3}z^2) exp(-2 b_i z) $$
% z .. depth under surface
% a .. 3x3 matrix with stress profile coeficients (in units of E/mm^j)
% b .. dim=3 vector with exponential decay coefficients in mm^-1
% data:   (eps_xx, error_xx, ... eps_zz, error_zz) table (6 columns)
% param:  [a,b]
% fix: flags for fixed parameters, array[12,1] of 0 or 1
% output: resulting parameters a,b
% std. deviation for fitted parameters
% fit   : fitted curves
% if ni==0, only evaluate the function

% initialization
afit=a0;
bfit=b0;
aerr=zeros(3,3);
berr=zeros(3,1);
Cov=0;

% check data format
nc=size(data,2);
np=size(data,1);
if ( (nc~=3) )
   fprintf('fit_stress_profile, wrong input data format:  ncol=%d',nc);
   fprintf('requires 3 columns with (x,y,error)');
   fit=zeros(np,1);
   Cov=0;
   return;
end

dx=data(:,1);
dy=data(:,2);
de=data(:,3);

par0=var2param(a0,b0);
p0=getfree(par0);
% if ni==0, only evaluate the function
if (ni==0) 
    fit=modelfnc(p0,dx);
    return;
end;

% data weighted by errors
yw = dy./de;

par=getparam(par0,p0);
[a,b]=param2var(par);
y=epsconv(npart,dx,S,a,b,beta,zeta,dsam);  

model = @(p,x) modelfnc(p,x)./de;
[p1,res,J,Cov] = nlinfit(dx,yw,model,p0);
[afit,bfit] = param2var(getparam(par0,p1));
[aerr,berr] = geterrors(p1,res,J);
fit=modelfnc(p1,dx);

    function y = modelfnc(p,x)
        par=getparam(par0,p);
        [a,b]=param2var(par);
        y=epsconv(npart,x,S,a,b,beta,zeta,dsam);        
    end

    function f = getfree(p)
       n=numel(p);
       j=0;
       x=zeros(1,n);
       for i=1:n
           if (fix(i)==0)
               j=j+1;
               x(j)=p(i);
           end
       end
       f=x(1:j);
    end

    function [ea, eb] = geterrors(p,res,J) 
        ytol = nlparci(p,res,'Jacobian',J);
        % convert 95% tolerance interval to std. dev.
        sd=abs(ytol(:,1)-ytol(:,2))/2;
        n=numel(par0);       
        err=getparam(zeros(n,1),sd);
        [ea, eb]=param2var(err);
    end

    function [par] = getparam(p,freepar)
       n=numel(p);
       j=0;
       par=p;
       for i=1:n
           if (fix(i)==0)
               j=j+1;
               par(i)=freepar(j);
           end;
       end;
    end
    
    function [a,b] = param2var(par)
       c=reshape(par,3,4);
       a=c(:,1:3);
       b=c(:,4);       
    end

    function [par] = var2param(a,b)
       c=reshape(a,1,9);
       par=cat(2,c,b);     
    end

end

