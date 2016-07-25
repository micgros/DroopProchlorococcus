%% Sensetivity analysis (Fig. 4)
% Grossowicz et al. 2016
% created by MG, micgros@gmail.com
% last midification : 25 Jul 2016

md=100;
R_CN=106/16;
radius=0.3628;  %"MED4"
vol = (4/3)*pi*radius^3;

rhoref = 1024.5;        % reference density of seawater
thetaK = 273.15 + 24.0; % temperature 
salt   = 34.5;          % salinity 
pt     = 0.0;           % total phosphorus
sit    = 0.0;           % total silicon
pCO2   = 400e-6;

load('ParamMat.mat')

for n=ParamMat(1,:)

gamma_n = n; 
gamma_c = 0.5;

% Nitrogen uptake (umol N cell-1 day-1)/86400 for per sec
% (9.1e-9)*vola^0.67 (umol N cell-1 day-1)
V_n_max = ((9.1e-9) / 86400) * vol^0.67; 
V_c_max = V_n_max; % 2.92e-8 / 86400;  

% growth rates. sec-1
mu_inf = 1e-5; 

% k - half saturation
k_n = 0.17 * vol^0.27; 
k_c = 0.17 * vol^0.27; 

% minimum N quota (umol N cell-1 )
q_n_min = (1.36e-9) * vol^0.77; % 5.93e-10;   
q_c_min = q_n_min * R_CN; % 2.176500000000000e-09; 

% maximum N quota (umol N cell-1 )
q_n_max=5 * q_n_min; 
q_c_max=3 * q_c_min; % 6.529500000000000e-09; 

% mortality - s-1
m = .1/86400; % 0.08/86400;  

% CO2 parameters
r0 = 0.18/86400;   % dark respiration, sec-1 = 0.18 d-1, Geider & Osborne 1989 
b  = 0.01;         % respiration coefficient, no units, Geider & Osborne 1989
Kg = 0.4968/86400; % m sec-1 = 0.7 cm h-1 from Cole & Caraco 1998
B  = 10;           % Revelle buffer factor, Mick's book
 
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(3,:)
gamma_n=0.5;
gamma_c=n;
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(5,:)
gamma_c=0.5; 
mu_inf=n*1e-5;  
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(7,:)
mu_inf=1e-5;
m=n*.1/86400; 
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(9,:)
m=.1/86400;
V_n_max=n*((9.1e-9)/86400)*vol^0.67;  
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(10,:)
V_n_max=((9.1e-9)/86400)*vol^0.67;  
V_c_max=n*((9.1e-9)/86400)*vol^0.67;
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(14,:)
V_c_max=((9.1e-9)/86400)*vol^0.67;
k_c=n*0.17*vol^0.27; 
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(18,:)
k_c=0.17*vol^0.27; 
k_n=n*0.17*vol^0.27; 
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(19,:)
k_n=0.17*vol^0.27; 
q_n_min=n*(1.36e-9)*vol^0.77; 
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(22,:)
q_n_min=(1.36e-9)*vol^0.77; 
q_c_min=n*q_n_min*R_CN; 
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(23,:)
q_c_min=q_n_min*R_CN; 
q_n_max=n*5*q_n_min;
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(24,:)
q_n_max=5*q_n_min;
q_c_max=n*3*q_c_min;
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(27,:)
q_c_max=3*q_c_min;
r0=n*0.18/86400; 
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(29,:)
r0=0.18/86400; 
b=n*0.01; 
forPCA9312_LowN_axenic %standart
end

% for n=ParamMat(31,:)
% ba=0.01;
% Csat =n* csatd * 1.0e6 / rhoref; 
% forPCA9312_Pro99_axenic %standart
% end
% 
% for n=ParamMat(32,:)
% Csat=csatd * 1.0e6 / rhoref;
% B =n*10;  
% forPCA9312_Pro99_axenic %standart
% end

for n=ParamMat(32,:)
b=0.01;
B=n*10;  
forPCA9312_LowN_axenic %standart
end

for n=ParamMat(30,:)
B=10;  
Kg=n*0.4968/86400;
forPCA9312_LowN_axenic %standart
end

