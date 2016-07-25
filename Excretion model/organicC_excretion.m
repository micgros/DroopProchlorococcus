function [ dOC_dt ] = organicC_excretion(m_star, b_c, gamma_c, epsilon_c)
% External organic C concentaration dOC/dt function
% MG 12.7.12
% Last modification 15 Feb 16, MG

dOC_dt = (1 - gamma_c) * m_star * b_c + epsilon_c*b_c;

end

