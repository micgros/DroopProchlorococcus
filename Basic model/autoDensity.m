function [dxa_dt ] = autoDensity( mu, m_star, x)
% Autotrophs density dX/dt function
% MG 12.7.12
% last modification: 15 Feb 16, MG

dxa_dt= (mu- m_star)*x;

end

