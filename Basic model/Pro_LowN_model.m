%% The basic model
% Grossowicz et al. 2016
% created by MG, micgros@gmail.com
% last modification : Droop-model ver 15 Feb 16
% MIT equations
% bicarbonate addition - 13.12.13
%%

% clear all variables
clear 
clc

% Number of days for run
md = 100;

% Redfield ratio
R_CN = 6.625;

% gama=fraction of heterotroph mort/resp to inorganic form
gamma_n = 0.5; 
gamma_c = 0.5; 

radius = 0.3628;  % "MED4" = 9312
vol = (4/3) * pi * radius^3;

% Nitrogen uptake (umol N cell-1 day-1)/86400 for per sec
% (9.1e-9)*vola^0.67 (umol N cell-1 day-1)
V_n_max = ((9.1e-9) / 86400) * vol^0.67; 
V_c_max = V_n_max; 

% growth rates. sec-1
mu_inf = 1e-5; 

% k - half saturation
k_n = 0.17 * vol^0.27; 
k_c = 0.17 * vol^0.27; 

% minimum N quota (umol N cell-1 )
q_n_min = (1.36e-9) * vol^0.77; 
q_c_min = q_n_min * R_CN; 

% maximum N quota (umol N cell-1 )
q_n_max=5 * q_n_min; 
q_c_max=3 * q_c_min; 

% mortality - s-1
m = .1/86400; 

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

% CO2 parameters
r0 = 0.18/86400;   % dark respiration, sec-1 = 0.18 d-1, Geider & Osborne 1989 
b  = 0.01;         % respiration coefficient, no units, Geider & Osborne 1989
Kg = 0.4968/86400; % m sec-1 = 0.7 cm h-1 from Cole & Caraco 1998
B  = 10;           % Revelle buffer factor, Mick's book
 
rhoref = 1024.5;        % reference density of seawater
thetaK = 273.15 + 24.0; % temperature 
salt   = 34.5;          % salinity 
pt     = 0.0;           % total phosphorus
sit    = 0.0;           % total silicon
pCO2   = 400e-6;

% time vector
time = linspace(0, md, 1000); 

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
for i = 1:86400*md
  
    h = (.09/-27/86400)*i+.105; % the hight of the water columns in meters
    
    % alkalinity and Csat calculation
    ta1     = (112.46/86400)*i + 2617.6;
    ta      = ta1 * rhoref / 1e6;
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

% end 

 %% plotting
 
% sample A results, Low N
xdot9312 = [0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27];
cell_data = [1419052.488   886758.691	2510054.533	4342450.579	6618183.367	9836826.858	13807174.51	22447085.89	40878834.36	55403033.4	70336826.86	NaN	101285105.7	110418967.3	96088786.64	96225630.54	NaN         87538854.81	83861877.98	75219154.74	71500937.29	NaN         66534417.34	56434688.35	47757181.57	39622764.23	34944173.44	29145799.46;
             1118204.628   1716336.772	2274463.34	3708809.59	6079871.759	8884165.04	11600641.2	18485503.21	25570393.09	40518817.95	NaN         NaN	90782408.7	112800111.5	99607053.25	105934650.2	77660181.2	97932422.4	81077528.59	96109906.43	87369523.24 84302539.73	78051091.64	61545522.06	53957077.08	45512995.69	34744244.76	30976533.49;
             1383908.046   1882758.621	2556475.096	3233409.962	5221302.682	7640459.77	NaN         18259003.83	25006896.55	39793256.7	NaN         NaN	91902222.22	116875862.1	119497777.8	118448582.4	NaN         101173842	85132970.03	87182425.07	80228201.63	73440871.93	61688692.1	49030517.71	34757765.67	32867574.93	28680926.43	23077588.56]*1000;
ydot9312 = nanmean(cell_data);
SD9312 = nanstd(cell_data);
 				

Nx = [0 2 4	6	8	9	11	13	14	15	16	23];
N_data = [209.2235294	132.3294118	119.8117647	116.2352941	103.7176471	41.12941176	27.65803922	0           0           68.66823529	0	0;
          194.3215686	80.47058824	114.4470588	125.1764706	78.68235294	91.2        43.87137255	73.91372549	34.33411765	82.61647059	0	0;
          156.172549	82.25882353	118.0235294	155.5764706	121.6       44.10980392	82.02039216	75.64235294	21.45882353	15.03905882	0	0];
Ny = nanmean(N_data);
N_SD = nanstd(N_data);

Cx = [0 7 13 19];
C_data = [2524	3492.7	3772.4	4300.4;
                2524	3244.5	3361.8	4181;
                2524	3697	3718.3	4274.7];
Cy = nanmean(C_data);
C_SD = nanstd(C_data);


figure
titles={'X','N','C','Q^N','Q^C','C:N ratio'};
subplot(2,3,1)
plot(time,X_vv)
hold on
errorbar(xdot9312,ydot9312,SD9312,'o','MarkerFaceColor','r')
title(titles{1})
ylabel('Density (cell L^{-1})')

subplot(2,3,2)
plot(time,N_vv)
hold on
errorbar(Nx,Ny,N_SD,'o','MarkerFaceColor','r')
title(titles{2})
ylabel('N (\mumol L^{-1})')

subplot(2,3,3)
plot(time,C_vv)
hold on
errorbar(Cx,Cy,C_SD,'o','MarkerFaceColor','r')
title(titles{3})
ylabel('C (\mumol L^{-1})')

subplot(2,3,4)
plot(time,q_n_vv)
title(titles{4})
ylabel('Q^N (\mumol cell^{-1})')
xlabel('Time (days)')

subplot(2,3,5)
plot(time,q_c_vv)
title(titles{5})
ylabel('Q^C (\mumol cell^{-1})')
xlabel('Time (days)')

subplot(2,3,6)
plot(time,q_c_vv./q_n_vv)
title(titles{6})
ylabel('C:N ratio')
xlabel('Time (days)')

