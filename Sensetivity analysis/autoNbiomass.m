function [ dba_n_dt ] = autoNbiomass( V_n_max, N, k_n, x, q_n, q_n_max, q_n_min, m_star, b_n)
% Autotrophs N biomass dB/dt function
% MG 12.7.12
% last modification: 15 Feb 16, MG

 dba_n_dt = (V_n_max * N * x * (q_n_max - q_n)) / ((N + k_n) * (q_n_max - q_n_min)) - ...
     m_star * b_n;

end

