%% Elementary transport matrix
% Operator for rotation by gamma and subsequent flight by L
%
% L = distance [mm]
% gamma = rotation angle [rad]
% lambda = nominal wavelength [A]
%
function [ C ] = CT( L, gamma, lambda )
    R=Rot(gamma);
    t0=L*lambda/3.95603; % nominal flight time
    Lmat=[L,0.0;0.0,0.0];
    C=[R, Lmat, zeros(2,1); zeros(2,2), eye(2,2), zeros(2,1); zeros(1,3), t0, 1];
end

