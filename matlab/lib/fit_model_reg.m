function [ params,errors,fit,chi2,Cov ] = fit_model_reg(data,regP,model,ni)
% Fit a model function represented by given object "model" to data
% Version with regularization term (must be implemented by the model)
% It will run one fitting per regP element. 
% data:   3-columns with x,y and error values
% regP    regularization coefficients
% model:  a model class implementing the fit function, fixed and free
% parameters etc.
% It must implement:
% model.setX(npart,x)     ...  set x-values and initialize the model
% model.fitFnc()          ...  fitting function
% model.getVariables()    ...  return an array of variable parameters
% model.setFree(p)        ...  set new values of free variables
% model.getFree()         ...  get free variables
% model.getFixed()        ...  return an array of (0|1) for free|fixed parameters
% model.vars              ...  pointer to the variables array
% params,errors: fitted parameters and their errors
% fit : fitted curves
% Cov: covariance matrix
% if ni==0, only evaluate the fit function, return zero errors and Cov
opts = statset('nlinfit');  
opts.MaxIter=200;
opts.TolFun=1e-10;
opts.TolX=1e-10;
opts.Robust='on';
opts.RobustWgtFun='cauchy';

% initialization
params=model.getVariables();
errors=zeros(size(params,1),size(params,2));
Cov=0;
chi2=0;

% check data format
nc=size(data,2);
np=size(data,1);
if ( (nc~=3) )
   fprintf('fit_model, wrong input data format:  ncol=%d\n',nc);
   fprintf('requires 3 columns with (x,y,error)');
   fit=zeros(np,1);
   Cov=0;
   return;
end

% if ni==0, only evaluate the function
if (ni==0) 
    fit=model.fitFnc();
    return;
end;

% get the 3 data columns
dx=data(:,1);
dy=data(:,2);
de=data(:,3);

% "fixed" flags
fix=model.getFixed();

% select free variables
pfree=model.getFree();

% data weighted by errors
yw = dy./de;

% declare model function
fitmodel = @(p,x) fitfnc(p,x);

% fit data
%[p1,res,J,Cov] = nlinfit(dx,yw,fitmodel,pfree);
nreg=numel(regP);
for ireg=1:nreg;
    creg=regP(ireg);
    [p1,res,J,Cov] = nlinfit(dx,yw,fitmodel,pfree,opts);
    fit=modelfnc(p1,dx);
    chi2=getchi2(numel(params),dy,de,fit); 
    reg=model.regValue();
    fprintf('iter %d: creg=%g, chi2=%g, reg=%g\n',ireg,creg,chi2,reg);
end;
% get parameters, errors and fit
model.setFree(p1);
params=model.vars;
errors = geterrors(p1,res,J);
fit=modelfnc(p1,dx);
chi2=getchi2(numel(params),dy,de,fit);  

    function y = modelfnc(p,~)
        % fprintf('modelfnc x=%f\n',x(1));
        model.setFree(p);
        y=model.fitFnc(); 
        %display(y);
        %fprintf('np=%d\n',numel(y));
    end

    function y = fitfnc(p,~)
        % fprintf('modelfnc x=%f\n',x(1));
        model.setFree(p);
        tmp=model.fitFnc()./de; 
        y=cat(1,tmp,creg*model.regValue()); 
        %display(y);
        %fprintf('np=%d\n',numel(y));
    end

% get complete set of parameters (fixed + free)
% p contains all default values, freepar only free parameters
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

    function err = geterrors(p,res,J) 
        ytol = nlparci(p,res,'Jacobian',J);
        % convert 95% tolerance interval to std. dev.
        n=numel(model.vars);       
        sd=abs(ytol(:,1)-ytol(:,2))/2;
        err=getparam(zeros(n,1),sd);
    end


    function chi2 = getchi2(np,y,err,yfit)      
        sd=(y(:,1)-yfit(:,1))./abs(err(:,1));
        sd2=sd.^2;
        nfree=size(y,1)-np-1;
        chi2=sum(sd2)/nfree;
    end

end

