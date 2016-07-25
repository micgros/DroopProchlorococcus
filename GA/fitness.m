% fitness(P)
%   Returns a vector with the fitness scores for each member of the
%   population, P
%
% Inputs:
%   P - A 2D array of the population genomes
%   gene_len - number of binary digits per number (1/2 chrom_len)
% Output:
%   fit = A 1D array of fitness scores

function fit = fitness(P)

% Use globals so we can send the values into our sir.m funciton
global gamma_N mu_inf m  V_N V_C q_N_min q_C_min q_N_max q_C_max r0 Kg xdot9312 ydot9312 Nx Ny Cx Cy

% I like to initialize any arrays I use
fit = zeros(size(P,1), 1);

% Loop through each individual in the population
for i=1:size(P,1)
    
    % Convert the current chromosome to real values.  I chose to interperet
    % the chromosome as two binary representations of numbers.  You may
    % also want to consider a Grey Code representation (or others).
   [gamma_N, mu_inf, m,  V_N, V_C, q_N_min, q_C_min, q_N_max, q_C_max, r0, Kg] = gene_to_values(P(i,:));
   
   % Get the model output for the given values of beta and delta
   % for LowN
%    X0 = [0.5148, 0.5148*6.625, ydot9312(1), 100, 20, 2000, 20*6.625]; % b0 for MED4 - 11.0842, b0 for 9312 - 0.5148
   % for  LowP
%    X0 = [0.0213, 0.0213*106, ydot9312(1), 50/8, 2000, 20*6.625/106, 20*6.625];
%    [t, y] = ode15s('Pro', [0 xdot9312(end)], X0); 

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

   
   % The fitness is related to the squared error between the model and
   % the data, so we use our squared_error.m function.
%    error1 = squared_error(xdot9312,ydot9312,t_all, y_all(:,3));
%    error2 = squared_error(Nx,Ny,t_all, y_all(:,4));
   error1 = squared_error(xdot9312,ydot9312,t_all, y_all(:,3));
   error2 = squared_error(Nx,Ny,t_all, y_all(:,4));
   error3 = squared_error(Cx,Cy,t_all, y_all(:,6));
   error=error1*error2*error3; %                                       
                      
    % The GA maximizes values, so we return the inverse of the error to
    % trick the GA into minimzing the error.
    fit(i) = error^-1;       
end