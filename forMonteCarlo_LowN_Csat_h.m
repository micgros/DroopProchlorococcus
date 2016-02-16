%% for Monte-Carlo ruffle
% created by MG, micgros@gmail.com
% last midification : 15 Feb 16

clc

% Number of days for run
md = 100;

% gama=fraction of heterotroph mort/resp to inorganic form
gamma_n = 0.5; 
l = find(X(1,:) == 1);  %if in X(1,:) there is 1, l is it location
if size(l) == 1        
gamma_n = X(2,l);
end
gamma_c=0.5; 

radius = 0.3628;  % "MED4" = 9312
vol = (4/3) * pi * radius^3;

% Nitrogen uptake (umol N cell-1 day-1)/86400 for per sec
V_n_max = ((9.1e-9) / 86400) * vol^0.67; 
l = find(X(1,:) == 4);  
if size(l) == 1        
V_n_max = X(2,l);
end

V_c_max = ((9.1e-9)/86400) * vol^0.67;
l = find(X(1,:) == 5);  
if size(l) == 1        
V_c_max = X(2,l);
end

% growth rates. sec-1
mu_inf=1e-5; 
l = find(X(1,:) == 2);  
if size(l) == 1        
mu_inf = X(2,l);
end

% k - half saturation
k_n = 0.17 * vol^0.27; 
k_c = 0.17 * vol^0.27; 

% minimum N quota (umol N cell-1 )
q_n_min = (1.36e-9) * vol^0.77; 
l = find(X(1,:) == 6);  
if size(l) == 1        
q_n_min = X(2,l);
end

q_c_min = (1.36e-9) * vol^0.77 * R_CN;
l = find(X(1,:) == 7);  
if size(l) == 1        
q_c_min = X(2,l);
end

% maximum N quota (umol N cell-1 )
q_n_max = 5 * q_n_min;
l = find(X(1,:) == 8);  
if size(l) == 1        
q_n_max = X(2,l);
end
q_c_max = 3 * q_c_min;
l = find(X(1,:) == 9);  
if size(l) == 1        
q_c_max = X(2,l);
end

% mortality - s-1
m = 7.1429e-7; 
l = find(X(1,:) == 3);  
if size(l)==1        
m=X(2,l);
end

% biomass (umol N l-1 )
b_n = 0.5148; % 9312
b_c = b_n * R_CN;
q_n = q_n_min;  
q_c = q_c_min;
x = b_n / q_n;

%nutrients
N = N0; % 800
C = C0; % 3000
ON = ON0; % 20
OC = OC0; % ON*R_CN

%CO2 parameters
r0 = 0.18/86400; % sec-1 = 0.18 d-1, Geider & Osborne 1989
l = find(X(1,:) == 10);  
if size(l) == 1        
r0 = X(2,l);
end

b = 0.01; % no units, Geider & Osborne 1989

Kg = 1.9444e-06; % m sec-1 = 0.7 cm h-1 from Cole & Caraco 1998
l = find(X(1,:) == 11);  
if size(l) == 1        
Kg = X(2,l);
end

B = 10; % Mick's book

rhoref = 1024.5 ;% reference density of seawater
thetaK = 273.15 + 24.0; % temperature 
salt = 34.5; % salinity 
pt = 0.0 ;% total phosphorus
sit = 0.0 ;% total silicon
pCO2 = 400e-6 ;

%time stepping
% dt=timestep (sec)
time=linspace(0,md,1000); 

%output vectors
% Ba_n_v is the vector used for keeping all of the ba_n values
B_n_v = zeros(86400*md, 1); 
B_c_v = zeros(86400*md, 1); 
X_v = zeros(86400*md, 1); 
N_v = zeros(86400*md, 1); 
C_v = zeros(86400*md, 1); 
ON_v = zeros(86400*md, 1); 
OC_v = zeros(86400*md, 1);
mu_v = zeros(86400*md, 1);
q_n_v = zeros(86400*md, 1);
q_c_v = zeros(86400*md, 1);

% vectors for graphs
B_n_vv=zeros(1, 1000); 
B_c_vv=zeros(1, 1000); 
X_vv=zeros(1, 1000); 
N_vv=zeros(1, 1000); 
C_vv=zeros(1, 1000); 
ON_vv=zeros(1, 1000); 
OC_vv=zeros(1, 1000);
mu_vv=zeros(1, 1000);
q_n_vv=zeros(1, 1000);
q_c_vv=zeros(1, 1000);

% The LOOP 
for i=1:86400*md
  h=(.09/-27/86400)*i+.105; % the hight of the water columns in meters
  
  ta = 112.46/86400 * i + 2617.6;
  ta = ta * rhoref / 1.0e6 ;
  [csatd] = calc_csat(thetaK, salt, pCO2, pt, sit, ta);
  Csat = csatd * 1.0e6 / rhoref; %umol kg-1

    % HCO3- additions (1 mM) on days: 0, 5, 11, 18
    if i==432000&&950400&&1555200
        C=C+1000;
    end

% growth rate
    limit = max(q_n_min/q_n, q_c_min/q_c); 
    mu = mu_inf*(1-limit);
    mu_v(i) = mu;

   % m star
if q_n < q_n_min || q_c < q_c_min
    m_star = m + r0 + b * mu;
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
N = N + inorganicN(V_n_max, N, k_n, x, q_n, q_n_max, q_n_min, gamma_n, m_star, b_n);
N_v(i) = N;
CO2_dis = 0.01 * C;
C = C + inorganicC(V_c_max, C, k_c, x, m, b_c, q_c, q_c_max, q_c_min, mu, gamma_c, Csat, r0, b, h, Kg, B, CO2_dis);
C_v(i) = C;

% organic nutrients
ON = ON + organicN(m_star, b_n, gamma_n);
ON_v(i) = ON;
OC = OC + organicC(m_star, b_c, gamma_c);
OC_v(i) = OC;

% quotas
% auto quotas
q_n = b_n / x;
q_n_v(i) = q_n;
q_c = b_c / x;
q_c_v(i) = q_c; 

end

% vectors for graphs loop
for k=1:1000
B_n_vv(k)=B_n_v(round(86.4*md*k)); 
B_c_vv(k)=B_c_v(round(86.4*md*k)); 
X_vv(k)=X_v(round(86.4*md*k)); 
N_vv(k)=N_v(round(86.4*md*k)); 
C_vv(k)=C_v(round(86.4*md*k)); 
ON_vv(k)=ON_v(round(86.4*md*k)); 
OC_vv(k)=OC_v(round(86.4*md*k));
mu_vv(k)=mu_v(round(86.4*md*k));
q_n_vv(k)=q_n_v(round(86.4*md*k));
q_c_vv(k)=q_c_v(round(86.4*md*k));
end


% save multiplotMonteCarlo_9312_LowN_Csat_h2 X_vv -append -ascii
% save N_9312_LowN_Csat_h2 N_vv -append -ascii
% save C_9312_LowN_Csat_h2 C_vv -append -ascii
% save QN_9312_LowN_Csat_h2 q_n_vv -append -ascii
% save QC_LowN_Csat_h2 q_c_vv -append -ascii

