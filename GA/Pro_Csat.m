% Grossowicz et al. 2016
% Last modify: 25 Jul 16
% Written by M. Grossowicz
% The model under ode15s solver


function [diff] = Pro_Csat(t,statVar)

global gamma_N mu_inf m V_N V_C q_N_min q_C_min q_N_max q_C_max r0 Kg 

diff=statVar;
% b_n=statVar(1); % N biomass
% b_c=statVar(2); % C biomass
% x=statVar(3); % X density
% N=statVar(4); % nitrogen conc.
% ON=statVar(5); % organic nitrogen conc.
% C=statVar(6); % carbon conc.
% OC=statVar(7); % organic carbon conc.

q_n=statVar(1)/statVar(3);
q_c=statVar(2)/statVar(3);
mu = mu_inf*(1-max(q_N_min/q_n,q_C_min/q_c));
CO2_dis=0.01*statVar(6);
h=(.09/-27)*t+.105; % the hight of the water columns in meters
rhoref=1024.5;
ta=112.46*t+2617.6;
ta=ta*rhoref/1e6;
[csatd]=calc_csat(273.15+24,34.5,400e-6,0,0,ta);
Csat=csatd*1e6/rhoref;

% MIT stars!!
if q_n<q_N_min || q_c<q_C_min
    m_star=m+r0+.01*mu;
else
    m_star=m;
end


% the differential equatons
diff(1)=(V_N*statVar(4)*statVar(3)*(q_N_max-q_n))/((statVar(4)+.1101)*(q_N_max-q_N_min))-m_star*statVar(1);  % k_N=.1101, k_P=2.1
diff(2)=(V_C*statVar(6)*statVar(3)*(q_C_max-q_c))/((statVar(6)+.1101)*(q_C_max-q_C_min))-m*statVar(2)-(r0+.01*mu)*statVar(2);
diff(3)=(mu- m_star)*statVar(3);
diff(4)=-V_N*(statVar(4)/(statVar(4)+.1101))*statVar(3)*((q_N_max-q_n)/(q_N_max-q_N_min))+gamma_N*m_star*statVar(1);
diff(5)=(1-gamma_N)*m_star*statVar(1);
diff(6)=-(V_C*statVar(6)*(q_C_max-q_c)*statVar(3))/((statVar(6)+.1101)*(q_C_max-q_C_min))+.5*m*statVar(2)+(r0+.01*mu)*statVar(2)-(statVar(6)-Csat)/((h*statVar(6))/(Kg*10*CO2_dis));
diff(7)=.5*m*statVar(2);

end
% ===========