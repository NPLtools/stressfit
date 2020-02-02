%% Elementary scattering  matrix
% Applies Bragg condition, see equation (A4)
%$$ \left( \begin{matrix} 1 & 2\tan(\theta ) \\ 0 & 1\end{matrix} \right) $$

%%
function [ SIG ] = Sigma( theta )
SIG=[1.0, 2*tan(theta); 0.0, 1.0];
end
