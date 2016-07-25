function [ dN_dt ] = inorganicN(V_n_max, N, k_n, x, q_n, q_n_max, q_n_min, gamma_n, m_star, b_n)
% External inorganic N concentaration dN/dt function
% MG 12.7.12
% last modification: 15 Feb 16, MG

dN_dt=-V_n_max * (N / (N + k_n)) * x * ((q_n_max - q_n) / (q_n_max - q_n_min)) + ...
    gamma_n * m_star * b_n;

end

