%% Monte-Carlo ruffle
% Grossowicz et al. 2016
% created by MG, micgros@gmail.com
% last modification : 15 Feb 16

clear;
clc

rng('shuffle');

% parameter vals:
mu_inf = 1e-5;
m = .1/86400;
V_N = 3.0958e-9/86400;  
V_C = 3.0958e-9/86400;
q_N_min = 3.9389e-10 ; 
q_N_max = 1.9694e-9; 
q_C_min = 2.1765e-9;  
q_C_max = 6.5295e-9;
r0 = .18/86400;
Kg = 0.4968/86400;
paramRange=[0, 1;
            mu_inf*.1,      mu_inf*5;
            m*.5,                m*5;
            V_N*.1,           V_N*10;
            V_C*.1,           V_C*10;
            q_N_min*.25,   q_N_min*2;
            q_C_min*.025,  q_C_min*2;
            q_N_max*.5,    q_N_max*2;
            q_C_max*.5,    q_C_max*2;
            r0*.1,             r0*10;
            Kg*0.0001,         Kg*10];

for u=1:1000

p=randi([2,5]);                      % Select how many parameters you want to change

rows=zeros(1,p);                     % which param were chosen  
Val=zeros(1,p);                      % random valus of chosen param 

% for continues algorithm
row = randperm(11);
for i = 1:p    
   Val(i) = (paramRange(row(i), 2) - paramRange(row(i), 1)) * rand(1, 1) + paramRange(row(i), 1);
   rows(i) = row(i);
end
Y = [rows; Val];                                    % which param with what value
if numel(Y) == 4
    X = [Y zeros(2,3)];
elseif numel(Y) == 6
    X = [Y zeros(2)];
elseif numel(Y) == 8
    X = [Y zeros(2,1)];
else
    X = Y;
end
save ParamVals_9312_LowN_axenic_Csat_h2 X -append -ascii  


R_CN = 6.625; % Redfeild
% nutrient starts at:
N0=800/8; 
C0=3000; % Csat +HCO3-
ON0=20; 
OC0=ON0*R_CN; 

forMonteCarlo_LowN_Csat_h
fprintf('Round #%d\n finished', u);

end





