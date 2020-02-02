function [ strainErr ] = genStressErrors(d,a0,b0,fix,Cov)
%genStrainErrors Generate error bars from covariance matrix

par0=var2param(a0,b0);
p0=getfree(par0);
nx=numel(d);
strainErr1=zeros(nx,3);
strainErr2=zeros(nx,3);
nt=20;
[U,S,~] = svd(Cov);
dev=(diag(S)).^0.5;
for it=1:nt
    rx=normrnd(0.0,dev);
    ry=U*rx;
    p=p0+ry';
	tmp=modelfnc(p,d);
    strainErr1=strainErr1+tmp;
    strainErr2=strainErr2+tmp.^2;    
end;
z=strainErr2./nt-(strainErr1./nt).^2;
strainErr=(z.^0.5);

    function y = modelfnc(p,x)
        par=getparam(par0,p);
        [a,b]=param2var(par);
        n=numel(x);
        y=stress(n,x,a,b);
    end

% get the full set of parameters from the free subset
% substitute from p for the fixed ones.
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
% get the reduced set of parameters (exclude fixed)
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

    function [a,b] = param2var(par)
       c=reshape(par,3,4);
       a=c(:,1:3);
       b=c(:,4);       
    end

    function [par] = var2param(a,b)
       c=reshape(a,1,9);
       par=cat(2,c,b');     
    end

end

