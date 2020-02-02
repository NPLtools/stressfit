%% Shift function
% Analytical function describing the difference between true and nominal gauge volume center of mass
%
% $$f(x,d)  =  - \pi ^{-{1/2}}\frac{  \exp(-x^2) - \exp(-(x-d )^2)} {erf(x) - erf(x-d )}  $$ 
% 
% $x$ = distance from surface (< 0 outside the sample), 
% $d$ = sample thickness 
%
% version: 1.2
% date: 4/4/2014
%

function [ val ] = f( x, d )
% approximation for large x
 if (x-d>5)
     val=x-d;
 elseif (x<-5)
     val=x;
 else        
    zjm= erf(x) - erf(x-d);
    val =  -1/sqrt(pi)*(exp(-x^2)-exp(-(x-d)^2))/zjm ;
 end;
end

