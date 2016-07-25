%% For sensetivity analysis (Fig. 4), the model
% Grossowicz et al. 2016
% created by MG, micgros@gmail.com
% last midification : 25 Jul 2016


% biomass (umol N l-1 )
b_n = 0.5148; % start at 0.5148 N per L
b_c = b_n * R_CN;

% start quotas (umol cell-1)
q_n = q_n_min;  
q_c = q_c_min;

% start cell density (cell L-1)
x = b_n / q_n;

% Nutrients (umol N l-1 )
N  = 800/8;
C  = 3000; % Csat + added 1 mmol HCO3- 
ON = 20; 
OC = ON*R_CN; 

time=linspace(0,md,1000); 

%output vectors
% Ba_n_v is the vector used for keeping all of the ba_n values
%output vectors
% Ba_n_v is the vector used for keeping all of the ba_n values
B_n_v = zeros(86400*md, 1); %preallocate
B_c_v = zeros(86400*md, 1); 
X_v   = zeros(86400*md, 1); 
N_v   = zeros(86400*md, 1); 
C_v   = zeros(86400*md, 1); 
ON_v  = zeros(86400*md, 1); 
OC_v  = zeros(86400*md, 1);
mu_v  = zeros(86400*md, 1);
q_n_v = zeros(86400*md, 1);
q_c_v = zeros(86400*md, 1);

% vectors for graphs
B_n_vv = zeros(1, 1000); 
B_c_vv = zeros(1, 1000); 
X_vv   = zeros(1, 1000); 
N_vv   = zeros(1, 1000); 
C_vv   = zeros(1, 1000); 
ON_vv  = zeros(1, 1000); 
OC_vv  = zeros(1, 1000);
mu_vv  = zeros(1, 1000);
q_n_vv = zeros(1, 1000);
q_c_vv = zeros(1, 1000);

% The LOOP 
for i=1:86400*md
    
    h = (.09/-27/86400)*i+.105; % the hight of the water columns in meters
   
    % alkalinity and Csat calculation
    ta      = (112.46/86400)*i  + 2617.6;
    ta      = ta * rhoref / 1e6;
    [csatd] = calc_csat(thetaK, salt, pCO2, pt, sit, ta);
    Csat    = csatd * 1e6 / rhoref;  

    % HCO3- additions (1 mM) on days: 5, 11, 18
    if i == 432000 
        C = C+1000;
    elseif i == 950400 
        C = C+1000;
    elseif i == 1555200
        C = C+1000;    
    end
    
   % growth rate
    limit   = max(q_n_min/q_n, q_c_min/q_c); 
    mu      = mu_inf*(1-limit);
    mu_v(i) = mu;

   % m star
   if q_n < q_n_min || q_c < q_c_min
        m_star = m + r0 + b*mu;
   else
        m_star = m;
   end
 
% biomass N auto
b_n = b_n + autoNbiomass(V_n_max, N, k_n, x, q_n, q_n_max, q_n_min, m_star, b_n);
B_n_v(i) = b_n;

% biomass C auto
b_c = b_c + autoCbiomass(V_c_max, C, k_c, x, q_c, q_c_max, q_c_min, m, b_c, r0, b, mu);
B_c_v(i) = b_c;
  
% Number density X auto
x = x + autoDensity(mu, m_star, x);
X_v(i) = x;

% Nutrients
% inorganic
N       = N + inorganicN(V_n_max, N, k_n, x, q_n, q_n_max, q_n_min, gamma_n, m_star, b_n);
N_v(i)  = N;
CO2_dis = 0.01 * C;
C       = C + inorganicC(V_c_max, C, k_c, x, m, b_c, q_c, q_c_max, q_c_min, mu, gamma_c, Csat, r0, b, h, Kg, B, CO2_dis);
C_v(i)  = C;

% organic nutrients
ON      = ON + organicN(m_star, b_n, gamma_n);
ON_v(i) = ON;
OC      = OC + organicC(m_star, b_c, gamma_c);
OC_v(i) = OC;

% quotas
% auto quotas
q_n      = b_n / x;
q_n_v(i) = q_n;
q_c      = b_c / x;
q_c_v(i) = q_c; 

end


% vectors for graphs loop
for k = 1:1000
B_n_vv(k) = B_n_v(round(86.4*md*k)); 
B_c_vv(k) = B_c_v(round(86.4*md*k)); 
X_vv(k)   = X_v(round(86.4*md*k)); 
N_vv(k)   = N_v(round(86.4*md*k)); 
C_vv(k)   = C_v(round(86.4*md*k)); 
ON_vv(k)  = ON_v(round(86.4*md*k)); 
OC_vv(k)  = OC_v(round(86.4*md*k));
mu_vv(k)  = mu_v(round(86.4*md*k));
q_n_vv(k) = q_n_v(round(86.4*md*k));
q_c_vv(k) = q_c_v(round(86.4*md*k));
end

save multiplotfrom9312_Low_N_Csat_h X_vv -append -ascii

