function [dON_dt] = organicN(m_star, b_n, gamma_n)
% External organic N concentaration dON/dt function
% MG 12.7.12
% last modification: 15 Feb 16, MG

dON_dt = (1-gamma_n) * m_star * b_n;

end

