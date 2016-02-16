function [ dC_dt ] = inorganicC( V_c_max, C, k_c, x, m, b_c, q_c, q_c_max, q_c_min, mu, gamma_c, Csat, r0, b, h, Kg, B, CO2_dis )
% External inorganic C concentaration dC/dt function
% MG 16.4.13
% last modification: 15 Feb 16, MG

dC_dt = -(V_c_max * C * (q_c_max - q_c) * x) / ((C + k_c) * (q_c_max - q_c_min)) + ...
    gamma_c * m * b_c + (r0 + b * mu) * b_c - (C - Csat) / ((h * C) / (Kg * B * CO2_dis));

end


