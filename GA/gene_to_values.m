% gene_to_values(gene)
%   Given a binary array chromosome, returns real values for beta and delta
%
% Inputs:
%   gene  - A 1D binary array, the first half representing beta and
%           the second half representing delta.
% Output:
%   beta  - The value for beta
%   delta - The value for delta
%
function [gamma_N, mu_inf, m, V_N, V_C, q_N_min, q_C_min, q_N_max, q_C_max, r0, Kg] = gene_to_values(gene) 

mu_inf_guess=.864;
m_guess=.1;
V_N_guess=3.0958e-9; % N=3.0958e-9, P=1.2e-11 
V_C_guess=3.0958e-9;
q_N_min_guess=3.9389e-10; % N=3.9389e-10, P=2.0533e-11 
q_N_max_guess=1.9694e-9; % N=1.9694e-09, P=1.0266e-10
q_C_min_guess=2.1765e-10;  % real value=2.1765e-9
q_C_max_guess=6.5295e-9;
r0_guess=.18;
Kg_guess=0.4968;

% Set the min and max ranges for the exponents (base 10)
gamma_N_min_exponent = log10(0.00001); % 0<=gamma_N<=1
gamma_N_max_exponent = log10(1);
mu_inf_min_exponent = log10(mu_inf_guess*.1);
mu_inf_max_exponent = log10(mu_inf_guess*5);
m_min_exponent = log10(m_guess*.5);
m_max_exponent = log10(m_guess*5);
V_N_min_exponent = log10(V_N_guess*.1);
V_N_max_exponent = log10(V_N_guess*10);
V_C_min_exponent = log10(V_C_guess*.1);
V_C_max_exponent = log10(V_C_guess*10);
q_N_min_min_exponent = log10(q_N_min_guess*.25);
q_N_min_max_exponent = log10(q_N_min_guess*2);
q_N_max_min_exponent = log10(q_N_max_guess*.5);
q_N_max_max_exponent = log10(q_N_max_guess*2);
q_C_min_min_exponent = log10(q_C_min_guess*.025); % .025
q_C_min_max_exponent = log10(q_C_min_guess*2); % 2
q_C_max_min_exponent = log10(q_C_max_guess*.5);
q_C_max_max_exponent = log10(q_C_max_guess*2);
r0_min_exponent = log10(r0_guess*.1);
r0_max_exponent = log10(r0_guess*10);
Kg_min_exponent = log10(Kg_guess*.1);
Kg_max_exponent = log10(Kg_guess*10);

% The length of a gene within the chromosome is half the total length
gene_len = size(gene, 2)/11;

% Use bi2de to convert binary arrays to decimal numbers
gamma_N_decimal = bi2de(gene(1:gene_len));
mu_inf_decimal = bi2de(gene((gene_len+1):2*gene_len));
m_decimal = bi2de(gene((2*gene_len+1):3*gene_len));
V_N_decimal = bi2de(gene((3*gene_len+1):4*gene_len));
V_C_decimal = bi2de(gene((4*gene_len+1):5*gene_len));
q_N_min_decimal = bi2de(gene((5*gene_len+1):6*gene_len));
q_C_min_decimal = bi2de(gene((6*gene_len+1):7*gene_len));
q_N_max_decimal = bi2de(gene((7*gene_len+1):8*gene_len));
q_C_max_decimal = bi2de(gene((8*gene_len+1):9*gene_len));
r0_decimal = bi2de(gene((9*gene_len+1):10*gene_len));
Kg_decimal = bi2de(gene((10*gene_len+1):end));
    
% Convert the decimal values to real numbers inside the
% ranges established above.
gamma_N_exponent = gamma_N_min_exponent + ...
                  ((gamma_N_max_exponent-gamma_N_min_exponent) * ...
                   gamma_N_decimal / 2^gene_len);
mu_inf_exponent = mu_inf_min_exponent + ...
                  ((mu_inf_max_exponent-mu_inf_min_exponent) * ...
                   mu_inf_decimal / 2^gene_len);
m_exponent = m_min_exponent + ...
                  ((m_max_exponent-m_min_exponent) * ...
                   m_decimal / 2^gene_len);
V_N_exponent = V_N_min_exponent + ...
                  ((V_N_max_exponent-V_N_min_exponent) * ...
                   V_N_decimal / 2^gene_len); 
V_C_exponent = V_C_min_exponent + ...
                  ((V_C_max_exponent-V_C_min_exponent) * ...
                   V_C_decimal / 2^gene_len);                
q_N_min_exponent = q_N_min_min_exponent + ...
                  ((q_N_min_max_exponent-q_N_min_min_exponent) * ...
                   q_N_min_decimal / 2^gene_len);  
q_C_min_exponent = q_C_min_min_exponent + ...
                  ((q_C_min_max_exponent-q_C_min_min_exponent) * ...
                   q_C_min_decimal / 2^gene_len); 
q_N_max_exponent = q_N_max_min_exponent + ...
                  ((q_N_max_max_exponent-q_N_max_min_exponent) * ...
                   q_N_max_decimal / 2^gene_len);                
q_C_max_exponent = q_C_max_min_exponent + ...
                  ((q_C_max_max_exponent-q_C_max_min_exponent) * ...
                   q_C_max_decimal / 2^gene_len); 
r0_exponent = r0_min_exponent + ...
                  ((r0_max_exponent-r0_min_exponent) * ...
                   r0_decimal / 2^gene_len);
Kg_exponent = Kg_min_exponent + ...
                  ((Kg_max_exponent-Kg_min_exponent) * ...
                   Kg_decimal / 2^gene_len);               
% We are using base 10 exponents               
gamma_N = 10^gamma_N_exponent;
mu_inf = 10^mu_inf_exponent;
m = 10^m_exponent;
V_N = 10^V_N_exponent;
V_C = 10^V_C_exponent;
q_N_min = 10^q_N_min_exponent;
q_C_min = 10^q_C_min_exponent;
q_N_max = 10^q_N_max_exponent;
q_C_max = 10^q_C_max_exponent;
r0 = 10^r0_exponent;
Kg = 10^Kg_exponent;


