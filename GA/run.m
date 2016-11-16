% Grossowicz et al. 2016
% Last modify: 25 Jul 16
% Written by M. Grossowicz
% Run genetic algorithm. The output is the parameters values (on Workspace) and a plot 

clear

global gamma_N mu_inf m  V_N V_C q_N_min q_C_min q_N_max q_C_max r0 Kg xdot9312 ydot9312 Nx Ny Cx Cy

    
% Experiment results
xdot9312 = [0,1,2,3,4,5,7,8,9,12,13,14,16,17,18,19,21,22,23,24,25,26,27]';
ydot9312 = [1307055054,1828618028,2446997656.33333,3761556710.33333,5973119269.33333,8787150556,19730530976.6667,30485374666.6667,45238369350,94656578873.3333,113364980300,105064539230,106869621046.667,95548373070,83357458866.6667,86170495413.3333,79699554053.3333,68758067026.6667,55670242706.6667,45490674773.3333,39334444950,32789781543.3333,27733307170]';
SD9312   = [164490563.562652,97258986.6115434,151211060.778527,556398657.948174,704532456.717958,1101392733.72024,2355329794.34976,9005408673.97790,8810329518.17906,5767713715.04979,3265299402.11872,12622730754.4232,11140939152.7850,7123270830.43816,2074242102.48932,10482074115.6434,7947490560.44554,8404789891.68722,6292425313.52882,9798274309.74665,6327638781.94501,3559776766.32198,4134571755.78986]';
Nx       = [0     2     4     6     8     9    11    13    14    15    16    23];
Ny       = [100 98.3529  117.4275  132.3294  101.3333   58.8131   51.1833   49.8520   18.5976   55.4412   0  0];
NSD      = [0   29.4381  2.7316    20.6229   21.5579    28.0875   27.9090   43.1818   17.3450   35.6776   0  0];
Cx       = [0         7         13        19];
Cy       = [2.5240    3.4781    3.6175    4.2520]*1000;
CSD      = [0         226.6046  223.0887  62.8444];

[P, best] = ga(50, 110, 0.01, 0.5, 100); % run GA
[gamma_N, mu_inf, m, V_N, V_C, q_N_min, q_C_min, q_N_max, q_C_max, r0, Kg] = gene_to_values(best); % get values of best individual

% Run model with GA results
X0 = [0.5148, 0.5148*6.625, ydot9312(1), 100, 20, 3000, 20*6.625]; % b0 for MED4 - 11.0842, b0 for 9312 - 0.5148
[t, y] = ode15s('Pro_Csat', [0 5] , X0); % xdot9312(end)]
t_all=t; t_all(end)=[];
y_all=y; y_all(end,:)=[];

% HCO3- additions (1 mM) on days: 5, 11, 18 
X0=[y(end,1:5) y(end,6)+1000 y(end,7)];
[t,y]=ode15s('Pro_Csat', [5 11] , X0);
t_all=[t_all; t]; t_all(end)=[];
y_all=[y_all; y]; y_all(end,:)=[];

X0=[y(end,1:5) y(end,6)+1000 y(end,7)];
[t,y]=ode15s('Pro_Csat', [11 18] , X0);
t_all=[t_all; t]; t_all(end)=[];
y_all=[y_all; y]; y_all(end,:)=[];

X0=[y(end,1:5) y(end,6)+1000 y(end,7)];
[t,y]=ode15s('Pro_Csat', [18 xdot9312(end)] , X0);
t_all=[t_all; t];
y_all=[y_all; y]; 

% figure
titles={'B^N','B^C','X','N','ON','C','OC','Q^N','Q^C'};
figure
subplot(3,3,3)
plot(t_all, y_all(:,3))
hold on
errorbar(xdot9312, ydot9312, SD9312,'rs')
title(titles{3})
subplot(3,3,4)
plot(t_all, y_all(:,4))
hold on
errorbar(Nx,Ny,NSD,'og')
title(titles{4})
subplot(3,3,6)
plot(t_all, y_all(:,6))
hold on
title(titles{6})
errorbar(Cx,Cy,CSD,'or')
for j=[1 2 5 7]
    subplot(3,3,j)
    plot(t_all,y_all(:,j))
    title(titles{j})
end
Q_N_min=ones(length(t_all),1)*q_N_min; Q_N_max=ones(length(t_all),1)*q_N_max; Q_C_min=ones(length(t_all),1)*q_C_min; Q_C_max=ones(length(t_all),1)*q_C_max;
subplot(3,3,8)
plot(t_all,y_all(:,1)./y_all(:,3),t_all,Q_N_min,'k',t_all,Q_N_max,'k')
title(titles{8})
subplot(3,3,9)
plot(t_all,y_all(:,2)./y_all(:,3),t_all,Q_C_min,t_all,Q_C_max)
title(titles{9})
