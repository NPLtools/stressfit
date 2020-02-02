function [beta,r,J,Sigma,mse] = nlinfit2(X,y,model,rmodel,beta,options)
%NLINFIT Nonlinear least-squares regression.
%   BETA = NLINFIT(X,Y,MODELFUN,BETA0) estimates the coefficients of a
%   nonlinear regression function, using least squares estimation.  Y is a
%   vector of response (dependent variable) values.  Typically, X is a
%   design matrix of predictor (independent variable) values, with one row
%   for each value in Y and one column for each coefficient.  However, X
%   may be any array that MODELFUN is prepared to accept.  MODELFUN is a
%   function, specified using @, that accepts two arguments, a coefficient
%   vector and the array X, and returns a vector of fitted Y values.  BETA0
%   is a vector containing initial values for the coefficients.
%
%   [BETA,R,J,COVB,MSE] = NLINFIT(X,Y,MODELFUN,BETA0) returns the fitted
%   coefficients BETA, the residuals R, the Jacobian J of MODELFUN, the
%   estimated covariance matrix COVB for the fitted coefficients, and an
%   estimate MSE of the variance of the error term.  You can use these
%   outputs with NLPREDCI to produce confidence intervals for predictions,
%   and with NLPARCI to produce confidence intervals for the estimated
%   coefficients.  If you use a robust option (see below), you must use
%   COVB and may need MSE as input to NLPREDCI or NLPARCI to insure that
%   the confidence intervals take the robust fit properly into account.
%
%   [...] = NLINFIT(X,Y,MODELFUN,BETA0,OPTIONS) specifies control parameters
%   for the algorithm used in NLINFIT.  OPTIONS is a structure that can be
%   created by a call to STATSET.  Applicable STATSET parameters are:
%
%      'MaxIter'     - Maximum number of iterations allowed.  Defaults to 100.
%      'TolFun'      - Termination tolerance on the residual sum of squares.
%                      Defaults to 1e-8.
%      'TolX'        - Termination tolerance on the estimated coefficients
%                      BETA.  Defaults to 1e-8.
%      'Display'     - Level of display output during estimation.  Choices
%                      are 'off' (the default), 'iter', or 'final'.
%      'DerivStep'   - Relative difference used in finite difference gradient
%                      calculation.  May be a scalar, or the same size as
%                      the parameter vector BETA.  Defaults to EPS^(1/3).
%      'FunValCheck' - Check for invalid values, such as NaN or Inf, from
%                      the objective function.  'off' or 'on' (default).
%      'Robust'      - Flag to invoke the robust fitting option.  'off' (the default)
%                      or 'on'.
%      'WgtFun'      - A weight function for robust fitting.  Valid only when Robust
%                      is 'on'.  'bisquare' (the default), 'andrews', 'cauchy',
%                      'fair', 'huber', 'logistic', 'talwar', or 'welsch'.  Can
%                      also be a function handle that accepts a normalized residual
%                      as input and returns the robust weights as output.
%      'Tune'        - The tuning constant used in robust fitting to normalize the
%                      residuals before applying the weight function.  A positive
%                      scalar.  The default value depends upon the weight function.
%                      This parameter is required if the weight function is
%                      specified as a function handle.
%
%   NLINFIT treats NaNs in Y or MODELFUN(BETA0,X) as missing data, and
%   ignores the corresponding observations.
%
%   Examples:
%
%      Use @ to specify MODELFUN:
%         load reaction;
%         beta = nlinfit(reactants,rate,@mymodel,beta);
%
%      where MYMODEL is a MATLAB function such as:
%         function yhat = mymodel(beta, x)
%         yhat = (beta(1)*x(:,2) - x(:,3)/beta(5)) ./ ...
%                        (1+beta(2)*x(:,1)+beta(3)*x(:,2)+beta(4)*x(:,3));
%   
%      For an example of weighted fitting, see the Statistics Toolbox demo
%      "Weighted Nonlinear Regression".
%
%   See also NLPARCI, NLPREDCI, NLMEFIT, NLINTOOL, STATSET.

%   References:
%      [1] Seber, G.A.F, and Wild, C.J. (1989) Nonlinear Regression, Wiley.

%   NLINFIT can be used to make a weighted fit with known weights:
%
%      load reaction;
%      w = [8 2 1 6 12 9 12 10 10 12 2 10 8]'; % some example known weights
%      ratew = sqrt(w).*rate;
%      mymodelw = @(beta,X) sqrt(w).*mymodel(beta,X);
% 
%      [betaw,residw,Jw] = nlinfit(reactants,ratew,mymodelw,beta);
%      betaciw = nlparci(betaw,residw,Jw);
%      [ratefitw, deltafitw] = nlpredci(@mymodel,reactants,betaw,residw,Jw);
%      rmse = norm(residw) / (length(w)-length(rate))
% 
%   Predict at the observed x values.  However, the prediction band
%   assumes a weight (measurement precision) of 1 at these points.
%
%      [ratepredw, deltapredw] = ...
%            nlpredci(@mymodel,reactants,betaw,residw,Jw,[],[],'observation');

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.7.2.1 $  $Date: 2011/06/20 21:23:52 $

if nargin < 4
    error(message('stats:nlinfit:TooFewInputs'));
elseif ~isvector(y)
    error(message('stats:nlinfit:NonVectorY'));
end
if nargin < 5
    options = statset('nlinfit');
else
    options = statset(statset('nlinfit'), options);
end

% Check sizes of the model function's outputs while initializing the fitted
% values, residuals, and SSE at the given starting coefficient values.
model = fcnchk(model);
try
    yfit = model(beta,X);
catch ME
    if isa(model, 'inline')
        m = message('stats:nlinfit:ModelInlineError');
        throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
    elseif strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
                && ~isempty(strfind(ME.message, func2str(model)))
        error(message('stats:nlinfit:ModelFunctionNotFound', func2str( model )));
    else
        m = message('stats:nlinfit:ModelFunctionError',func2str(model));
        throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
    end
end
if ~isequal(size(yfit), size(y))
    expSizeStr = sprintf('%-d-by-',size(y));
    actSizeStr = sprintf('%-d-by-',size(yfit));
    error(message('stats:nlinfit:WrongSizeFunOutput',expSizeStr(1:end-4),actSizeStr(1:end-4)));
end

% Find NaNs in either the responses or in the fitted values at the starting
% point.  Since X is allowed to be anything, we can't just check for rows
% with NaNs, so checking yhat is the appropriate thing.  Those positions in
% the fit will be ignored as missing values.  NaNs that show up anywhere
% else during iteration will be treated as bad values.
nans = (isnan(y(:)) | isnan(yfit(:))); % a col vector
r = y(:) - yfit(:);
r(nans) = [];
n = numel(r);
p = numel(beta);
sse = r'*r;

funValCheck = strcmp(options.FunValCheck, 'on');
if funValCheck && ~isfinite(sse), checkFunVals(r); end

% Set the level of display
switch options.Display
    case 'off',    verbose = 0;
    case 'notify', verbose = 1;
    case 'final',  verbose = 2;
    case 'iter',   verbose = 3;
end
maxiter = options.MaxIter;

if strcmp(options.Robust,'off')
    % display format for non-robust fit
    if verbose > 2 % iter
        disp(' ');
        disp('                                     Norm of         Norm of');
        disp('   Iteration             SSE        Gradient           Step ');
        disp('  -----------------------------------------------------------');
        disp(sprintf('      %6d    %12g',0,sse)); %#ok<DSPS>
    end
    [beta,J,~,cause] = LMfit(X,y, model,rmodel,beta,options,verbose,maxiter);
else
    % Do a preliminary fit just to get residuals and leverage from the
    % least squares coefficient.  We won't count this against the iteration
    % limit.
    [beta_ls,J] = LMfit(X,y, model,rmodel,beta,options,0,maxiter);
    res = y - model(beta_ls,X);
    res(isnan(res)) = [];
    ols_s = norm(res) / sqrt(max(1,length(res)-numel(beta)));

    % display format for robust fit
    % Please note there are two loops for robust fit. It would be very
    % confusing if we display all iteration results. Instead, only the last
    % step of each inner loop (LM fit) will be output.
    if verbose > 2 % iter
        disp(' ');
        disp('Displaying iterations that re-calculate the robust weights');
        disp(' ');
        disp('   Iteration             SSE ');
        disp('  -----------------------------');
        disp(sprintf('      %6d    %12g',0,sse)); %#ok<DSPS>
    end
    [beta,J,sig,cause] = nlrobustfit(X,y,beta,model,rmodel,J,ols_s,options,verbose,maxiter);  
end;

switch(cause)
    case 'maxiter'
        warning(message('stats:nlinfit:IterationLimitExceeded'));
    case 'tolx'
        if verbose > 1 % 'final' or 'iter'
            disp('Iterations terminated: relative norm of the current step is less than OPTIONS.TolX');
        end
    case 'tolfun'
        if verbose > 1 % 'final' or 'iter'
            disp('Iterations terminated: relative change in SSE less than OPTIONS.TolFun');
        end
    case 'stall'
        warning(message('stats:nlinfit:UnableToDecreaseSSE'));
end

% If the Jacobian is ill-conditioned, then it's likely that two parameters
% are aliased and the estimates will be highly correlated.  Prediction at
% new x values not in the same column space is dubious.  NLPARCI will have
% trouble computing CIs because the inverse of J'*J is difficult to get
% accurately. NLPREDCI will have the same difficulty, and in addition, 
% will in effect end up taking the difference of two very large, but nearly 
% equal, variance and covariance terms, lose precision, and so the
% prediction bands will be erratic. It may also be that the Jacobian has
% one or more columns of zeros, meaning model is constant with respect to
% one or more parameters.  This may be because those parameters are not
% even in the expression in the model function, or they are multiplied by
% another param that is estimated at exactly zero (or something similar),
% or because some part of the model function is underflowing, making it a
% constant zero.

[~,R] = qr(J,0);
if n <= p
     warning(message('stats:nlinfit:Overparameterized')); 
elseif condest(R) > 1/sqrt(eps(class(beta)))
     if any(all(abs(J)<sqrt(eps(norm(J,1))),1),2) % one or more columns of zeros
         warning(message('stats:nlinfit:ModelConstantWRTParam'));
     else % no columns of zeros
         warning(message('stats:nlinfit:IllConditionedJacobian'));
     end
end

if nargout > 1
    % Return residuals and Jacobian that have missing values where needed.
    yfit = model(beta,X);
    r = y - yfit;
    JJ(~nans,:) = J;
    JJ(nans,:) = NaN;
    J = JJ;  
end
if nargout > 3
    if strcmp(options.Robust,'off')
        % We could estimate the population variance and the covariance matrix
        % for beta here as
        mse = sum(abs(r(~nans)).^2)/(n-p);
    else
        mse = sig.^2;
    end
    Rinv = inv(R);
    Sigma = Rinv*Rinv'*mse; %#ok<MINV>
end
end % function nlinfit

%----------------------------------------------------------------------
function  [beta,J,iter,cause] = LMfit(X,y, model,rmodel,beta,options,verbose,maxiter) 
% Levenberg-Marquardt algorithm for nonlinear regression

% Set up convergence tolerances from options.
betatol = options.TolX;
rtol = options.TolFun;
fdiffstep = options.DerivStep;
if isscalar(fdiffstep)
    fdiffstep = repmat(fdiffstep, size(beta));
else
    % statset ensures fdiffstep is not a matrix.
    % Here, we ensure fdiffstep has the same shape as beta.
    fdiffstep = reshape(fdiffstep, size(beta));
end
funValCheck = strcmp(options.FunValCheck, 'on');

% Set initial weight for LM algorithm.
lambda = .01;

% Set the iteration step
sqrteps = sqrt(eps(class(beta)));

% treatment for nans
yfit = model(beta,X);
r = y(:) - yfit(:);
nans = (isnan(y(:)) | isnan(yfit(:))); % a col vector
r(nans) = [];

reg=rmodel.getReg(beta);
sse = r'*r + reg;

iter = 0;
breakOut = false;
cause = '';

while iter < maxiter
    iter = iter + 1;
    betaold = beta;
    sseold = sse;

    % Compute a finite difference approximation to the Jacobian
    J = getjacobian(beta,fdiffstep,model,X,yfit,nans);

    % Levenberg-Marquardt step: inv(J'*J+lambda*D)*J'*r
    H=J'*J;
    diagJtJ = diag(H);    
    if funValCheck && ~all(isfinite(diagJtJ)), checkFunVals(J(:)); end
    diagH=diag(diagJtJ);
    regMat=rmodel.getRegMat;
    Jplus = H + lambda*diagH + regMat;
    rplus = J'*r;
    step = Jplus \ rplus;
    beta(:) = beta(:) + step;

    % Evaluate the fitted values at the new coefficients and
    % compute the residuals and the SSE.
    yfit = model(beta,X);
    r = y(:) - yfit(:);
    r(nans) = [];
    reg=rmodel.getReg(beta);
    sse = r'*r + reg;
    if funValCheck && ~isfinite(sse), checkFunVals(r); end
    % If the LM step decreased the SSE, decrease lambda to downweight the
    % steepest descent direction.  Prevent underflowing to zero after many
    % successful steps; smaller than eps is effectively zero anyway.
    if sse < sseold
        lambda = max(0.1*lambda,eps);
        
    % If the LM step increased the SSE, repeatedly increase lambda to
    % upweight the steepest descent direction and decrease the step size
    % until we get a step that does decrease SSE.
    else
        while sse > sseold
            lambda = 10*lambda;
            if lambda > 1e16
                breakOut = true;
                break
            end
            Jplus = H + lambda*diagH + regMat;
            step = Jplus \ rplus;    
            beta(:) = betaold(:) + step;            
            yfit = model(beta,X);
            r = y(:) - yfit(:);
            r(nans) = [];
            reg=rmodel.getReg(beta);
            sse = r'*r + reg;
            if funValCheck && ~isfinite(sse), checkFunVals(r); end
        end
    end 
    if verbose > 2 % iter
        disp(sprintf('      %6d    %12g   %12g   %12g    %12g', ...
                     iter,sse,reg,norm(2*r'*J),norm(step))); %#ok<DSPS>
    end

    % Check step size and change in SSE for convergence.
    if norm(step) < betatol*(sqrteps+norm(beta))
        cause = 'tolx';
        break
    elseif abs(sse-sseold) <= rtol*sse
        cause = 'tolfun';
        break
    elseif breakOut
        cause = 'stall';
        break
    end
end
if (iter >= maxiter)
    cause = 'maxiter';
end
end % function LMfit

%--------------------------------------------------------------------------
function checkFunVals(v)
% check if the function has finite output
if any(~isfinite(v))
    error(message('stats:nlinfit:NonFiniteFunOutput'));
end
end % function checkFunVals

%--------------------------------------------------------------------------
function [beta,J,sig,cause]=nlrobustfit(x,y,beta,model,rmodel,J,ols_s,options,verbose,maxiter)
% nonlinear robust fit

tune = options.Tune;
WgtFun = options.WgtFun;
[WgtFun,tune] = statrobustwfun(WgtFun,tune);

yfit = model(beta,x);
fullr = y(:) - yfit(:);
ok = ~isnan(fullr);
r = fullr(ok);
Delta = sqrt(eps(class(beta)));

% Adjust residuals using leverage, as advised by DuMouchel & O'Brien
% Compute leverage based on X, the Jacobian
[Q,~]=qr(J,0);
h = min(.9999, sum(Q.*Q,2));

% Compute adjustment factor
adjfactor = 1 ./ sqrt(1-h);

radj = r .* adjfactor;

% If we get a perfect or near perfect fit, the whole idea of finding
% outliers by comparing them to the residual standard deviation becomes
% difficult.  We'll deal with that by never allowing our estimate of the
% standard deviation of the error term to get below a value that is a small
% fraction of the standard deviation of the raw response values.
tiny_s = 1e-6 * std(y);
if tiny_s==0
    tiny_s = 1;
end

% Main loop of repeated nonlinear fits, adjust weights each time
totiter = 0;
w = NaN(size(y));
while maxiter>0
    beta0=beta;
    s = madsigma(radj, numel(beta)); % robust estimate of sigma for residual

    % Compute robust weights based on current residuals
    w(ok) = feval(WgtFun, radj/(max(s,tiny_s)*tune));

    % this is the weighted nlinfit
    sw = sqrt(w);
    yw = y .* sw;
    modelw = @(b,x) sqrt(w).*model(b,x);
    [beta,~,lsiter,cause] = LMfit(x,yw,modelw,rmodel,beta0,options,0,maxiter); % 6th arg always silences display
    totiter = totiter + lsiter;
    maxiter = maxiter - lsiter;
    yfit = model(beta,x);
    fullr = y(:) - yfit(:);
    r = fullr(ok);
    radj = r .* adjfactor;
    
    % if there is no change in any coefficient, the iterations stop.
    if  all(abs(beta-beta0) < Delta*max(abs(beta),abs(beta0)))
        break;
    end
    
    if verbose > 2 % iter
        disp(sprintf('      %6d    %12g', ...
                     totiter, r'*r)); %#ok<DSPS>
    end
end

% this is a warning about the non-convergence
if maxiter<=0
    cause = 'maxiter';
end

% We need the Jacobian at the final coefficient estimates, but not the J1
% version returned by LMfit because it has robust weights included
fdiffstep = options.DerivStep;
if isscalar(fdiffstep)
    fdiffstep = repmat(fdiffstep,size(beta));
end
J = getjacobian(beta,fdiffstep,model,x,yfit,~ok);

% Compute MAD of adjusted residuals after dropping p-1 closest to 0
p = numel(beta);
n = length(radj);
mad_s = madsigma(radj, p);   

% Compute a robust scale estimate for the covariance matrix
sig = statrobustsigma(WgtFun,radj,p,mad_s,tune,h);

% Be conservative by not allowing this to be much bigger than the ols value
% if the sample size is not large compared to p^2
sig = max(sig, ...
          sqrt((ols_s^2 * p^2 + sig^2 * n) / (p^2 + n)));
end % function nlrobustfit

%----------------------- Robust estimate of sigma
function s = madsigma(r,p)
%MADSIGMA    Compute sigma estimate using MAD of residuals from 0
n = length(r);
rs = sort(abs(r));
s = median(rs(max(1,min(n,p)):end)) / 0.6745;
end % function madsigma

% ---------------------- Jacobian
function J = getjacobian(beta,fdiffstep,model,X,yfit,nans)
    function yplus = call_model_nested(betaNew)
        yplus = model(betaNew, X);
        yplus(nans) = [];
    end
J = statjacobian(@call_model_nested, beta, fdiffstep, yfit(~nans));
end % function getjacobian
