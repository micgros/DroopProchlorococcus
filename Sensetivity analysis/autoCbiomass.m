function [ dba_c_dt ] = autoCbiomass(V_c_max, C, k_c, x, q_c, q_c_max, q_c_min, m, b_c, r0, b, mu)
%Autotrophs C biomass dB/dt function
% MG 12.7.12
% last modification: 15 Feb 16, MG

dba_c_dt = (V_c_max * C * x * (q_c_max - q_c)) / ((C + k_c) * (q_c_max - q_c_min)) - ...
    m * b_c - (r0 + b * mu) * b_c;

end

