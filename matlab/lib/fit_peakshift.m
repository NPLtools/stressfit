function [ output,sdev, fit ] = fit_peakshift( data,param,fix,ni)
% Fit the instensity scan to the analytical model (Gaussian approx.)
% model: $$ y = y0 + \D_{SE} / \beta (zeta + f(\beta (z-z_0)-\zeta, \beta d) $$
%  $$ \zeta = a / {2 \beta} $$
% d is the sample thickness [mm]
% data:   (x,y,error) table
% param:  [y0 DSE beta a z0 d]
% output: resulting parameters
% sdev  : errors for fitted parameters
% fit   : fitted curve
%
% if ni==0, only evaluate the function
%
% version: 1.2
% date: 4/4/2014
%

p0=getfree(param);
if (ni==0) 
    fit=modelfnc(p0,data(:,1));
    output=param;
    n= size(data,1);
    sdev=zeros(n,1);
    return;
end    
w=data(:,3);
yw = data(:,2)./w;
model = @(p,x) modelfnc(p,x)./w;
[p1,res,J] = nlinfit(data(:,1),yw,model,p0);
output=getparam(param,p1);
sdev=geterrors(p1,res,J);
fit=modelfnc(p1,data(:,1));

    function y = modelfnc(p,x)
        par=getparam(param,p);
        [y0 DSE beta zeta h0 d]=param2var(par);
        n= size(x,1);
        y=zeros(n,1);
        for i=1:n
            h=beta*x(i)-zeta;
            %fprintf('fit_peakshift: h=%f, x=%f, beta=%f, h0=%f, d=%f \n',h,x(i),beta,h0,d);
            %y(i) = y0 + DSE/beta*(zeta+f(h-h0,d));
            y(i) = y0 + DSE*(PSModel.gcentre( h-h0, d, 0.0, 0.0 )/beta-x(i));
        end;
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

    function err = geterrors(p,res,J) 
        ytol = nlparci(p,res,'Jacobian',J);
        % convert 95% tolerance interval to std. dev.
        sd=abs(ytol(:,1)-ytol(:,2))/2;
        n=numel(param);       
        err=getparam(zeros(n,1),sd);
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
    
    function [y0 DSE beta zeta h0 d] = param2var(par)
       y0 =  par(1);
       DSE = par(2);
       beta = par(3);
       zeta=par(4)/2/beta;
       h0=beta*par(5);
       d=beta*par(6);
    end

end

