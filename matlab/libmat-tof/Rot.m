%% Elementary rotation matrix
% see equation (A2)
%
% $\gamma$ = rotation angle (applied to x,y only)
%

function [ R ] = Rot(gamma)
    ct=cos(gamma);
    st=sin(gamma);
    R=[ct, -st;st, ct];
end

