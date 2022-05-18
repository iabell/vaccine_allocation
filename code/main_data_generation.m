function main_data_generation

clear all
close all
hold off

%baseline model parameters
params = struct("contact_matrix",[2 2;1 1], "R0", 2, "theta", 1,...
    "gamma", 1,"pop_size", [5000;5000], "omega", [0.5;0.5],...
    "sigma", [0.5;0.5], "delta_A", 0.5, "I0",[50;50]);

params.alpha_1 = 0.5; % susceptibility 
params.alpha_2 = 0; % symptomatic infection
params.alpha_3 = 0; % severity 
params.alpha_4 = 0; % infectivity
params.max_vaccines = 5000;

%Calibrating beta under no vaccination
[beta_temp,temp] = r0_beta_calibration(params,[0;0],params.R0,-1);
params.beta = beta_temp;

params.dim = length(params.pop_size);

size = 100;

%discrete data options
contact_matrix_options = struct('opt1',[4 2; 2 1], 'opt2', [2 2; 1 1], 'opt3', [2 4; 1 2]);

omega_sigma_options = struct('opt1', [0.8; 0.2], 'opt2', [0.5; 0.5], 'opt3', [0.2; 0.8]);

delta_options = struct('opt1',0.1, 'opt2', 0.17, 'opt3', 0.5);

disease_example_data = true;
contact_matrix_data = false; 
omega_data = false; 
sigma_data = false; 
delta_A_data = true;
R0_data = true;


if disease_example_data
    illustrative_example(params) 
end

if contact_matrix_data
    params.alpha_1 = 0;
    params.alpha_2 = 0;
    params.alpha_3 = 0;
    params.alpha_4 = 0;
    
    %vaccine 1
    params.alpha_1 = 0.75;
    contact_matrix_1 = discrete_data_generation(params,'contact_matrix',contact_matrix_options,3);
    writematrix(contact_matrix_1,'contact_matrix_1.csv')
    
    %vaccine 2
    params.alpha_1 = 0;
    params.alpha_2 = 0.75;
    contact_matrix_2 = discrete_data_generation(params,'contact_matrix',contact_matrix_options,3);
    writematrix(contact_matrix_2,'contact_matrix_2.csv')
   
    %vaccine 3
    params.alpha_2 = 0;
    params.alpha_3 = 0.75;
    contact_matrix_3 = discrete_data_generation(params,'contact_matrix',contact_matrix_options,3);
    writematrix(contact_matrix_3,'contact_matrix_3.csv')

    %vaccine 4
    params.alpha_3 = 0;
    params.alpha_4 = 0.75;
    contact_matrix_4 = discrete_data_generation(params,'contact_matrix',contact_matrix_options,3);
    writematrix(contact_matrix_4,'contact_matrix_4.csv')
    
end

if omega_data 
    params.alpha_1 = 0;
    params.alpha_2 = 0;
    params.alpha_3 = 0;
    params.alpha_4 = 0;
    
    %vaccine 1
    params.alpha_1 = 0.75;
    omega_1 = discrete_data_generation(params,'omega',omega_sigma_options,3);
    writematrix(omega_1,'omega_data_1.csv')

    
    %vaccine 2
    params.alpha_1 = 0;
    params.alpha_2 = 0.75;
    omega_2 = discrete_data_generation(params,'omega',omega_sigma_options,3);
    writematrix(omega_2,'omega_data_2.csv')

    %vaccine 3
    params.alpha_2 = 0;
    params.alpha_3 = 0.75;
    omega_3 = discrete_data_generation(params,'omega',omega_sigma_options,3);
    writematrix(omega_3,'omega_data_3.csv')

    %vaccine 4
    params.alpha_3 = 0;
    params.alpha_4 = 0.75;
    omega_4 = discrete_data_generation(params,'omega',omega_sigma_options,3);
    writematrix(omega_4,'omega_data_4.csv')
end 

if sigma_data 
    params.alpha_1 = 0;
    params.alpha_2 = 0;
    params.alpha_3 = 0;
    params.alpha_4 = 0;
    
    %vaccine 1
    params.alpha_1 = 0.75;
    sigma_1 = discrete_data_generation(params,'sigma',omega_sigma_options,3);
    writematrix(sigma_1,'sigma_data_1.csv')

     
    %vaccine 2
    params.alpha_1 = 0;
    params.alpha_2 = 0.75;
    sigma_2 = discrete_data_generation(params,'sigma',omega_sigma_options,3);
    writematrix(sigma_2,'sigma_data_2.csv')

    %vaccine 3
    params.alpha_2 = 0;
    params.alpha_3 = 0.75;
    sigma_3 = discrete_data_generation(params,'sigma',omega_sigma_options,3);
    writematrix(sigma_3,'sigma_data_3_75.csv')

    %vaccine 4
    params.alpha_3 = 0;
    params.alpha_4 = 0.75;
    sigma_4 = discrete_data_generation(params,'sigma',omega_sigma_options,3);
    writematrix(sigma_4,'sigma_data_4.csv')
end 

if delta_A_data
    params.alpha_1 = 0;
    params.alpha_2 = 0;
    params.alpha_3 = 0;
    params.alpha_4 = 0;
    
    %vaccine 1
    params.alpha_1 = 0.75;
    delta_A_matrix_1 = continuous_data_generation(params,size,'delta_A',linspace(0,1,size));
    writematrix(delta_A_matrix_1, 'data_delta_A_1.csv')
    
    %vaccine 2
    params.alpha_1 = 0;
    params.alpha_2 = 0.75;
    delta_A_matrix_2 = continuous_data_generation(params,size,'delta_A',linspace(0,1,size));
    writematrix(delta_A_matrix_2, 'data_delta_A_2.csv')

    %vaccine 3
    params.alpha_2 = 0;
    params.alpha_3 = 0.75;
    delta_A_matrix_3 = continuous_data_generation(params,size,'delta_A',linspace(0,1,size));
    writematrix(delta_A_matrix_3, 'data_delta_A_3.csv')

    %vaccine 4
    params.alpha_3 = 0;
    params.alpha_4 = 0.75;
    delta_A_matrix_4 = continuous_data_generation(params,size,'delta_A',linspace(0,1,size));
    writematrix(delta_A_matrix_4, 'data_delta_A_4.csv')
end 

if R0_data 
    params.alpha_1 = 0;
    params.alpha_2 = 0;
    params.alpha_3 = 0;
    params.alpha_4 = 0;
    
    %vaccine 1
    params.alpha_1 = 0.75;
    R0_matrix_1 = continuous_data_generation(params,size,'R0',linspace(1,10,size));
    writematrix(R0_matrix_1,'data_R0_1.csv')

    %vaccine 2
    params.alpha_1 = 0;
    params.alpha_2 = 0.75;
    R0_matrix_2 = continuous_data_generation(params,size,'R0',linspace(1,10,size));
    writematrix(R0_matrix_2,'data_R0_2.csv')

    %vaccine 3
    params.alpha_2 = 0;
    params.alpha_3 = 0.75;
    R0_matrix_3 = continuous_data_generation(params,size,'R0',linspace(1,10,size));
    writematrix(R0_matrix_3,'data_R0_3.csv')

    %vaccine 4
    params.alpha_3 = 0;
    params.alpha_4 = 0.75;
    R0_matrix_4 = continuous_data_generation(params,size,'R0',linspace(1,10,size));
    writematrix(R0_matrix_4,'data_R0_4.csv')
end
end

%% Using specific parameters as example
function illustrative_example(params)
%define parameters for example disease 
params.contact_matrix = [2 1; 1 0.5];
params.R0 = 2;

%recalculate beta
[beta_temp,temp] = r0_beta_calibration(params,[0;0],params.R0,-1);
params.beta = beta_temp;

params.omega = [0.8; 0.2];
params.sigma = [0.8; 0.2];

size = 100;
ve_options = linspace(0,1,size);
ve_options = ve_options(1,2:end);
max_vaccines_options = struct('opt1',2000,'opt2', 5000, 'opt3', 8000);

%vaccine 1
params.alpha_1 = 0;
params.alpha_2 = 0;
params.alpha_3 = 0;
params.alpha_4 = 0;
ve_1 = continuous_data_generation(params,size-1,'alpha_1',ve_options);
writematrix(ve_1, 've_1.csv')

params.alpha_1 = 0.75;
max_vax_1 = discrete_data_generation(params,'max_vaccines',max_vaccines_options,3);
writematrix(max_vax_1, 'max_vax_1.csv')

%vaccine 2
params.alpha_1 = 0;
params.alpha_2 = 0;
params.alpha_3 = 0;
params.alpha_4 = 0;
ve_2 = continuous_data_generation(params,size-1,'alpha_2',ve_options);
writematrix(ve_2, 've_2.csv')

params.alpha_2 = 0.75;
max_vax_2 = discrete_data_generation(params,'max_vaccines',max_vaccines_options,3);
writematrix(max_vax_2, 'max_vax_2.csv')

%vaccine 3
params.alpha_1 = 0;
params.alpha_2 = 0;
params.alpha_3 = 0;
params.alpha_4 = 0;
ve_3 = continuous_data_generation(params,size-1,'alpha_3',ve_options);
writematrix(ve_3, 've_3.csv')

params.alpha_3 = 0.75;
max_vax_3 = discrete_data_generation(params,'max_vaccines',max_vaccines_options,3);
writematrix(max_vax_3, 'max_vax_3.csv')

%vaccine 4
params.alpha_1 = 0;
params.alpha_2 = 0;
params.alpha_3 = 0;
params.alpha_4 = 0;
ve_4 = continuous_data_generation(params,size-1,'alpha_4',ve_options);
writematrix(ve_4, 've_4.csv')

params.alpha_4 = 0.75;
max_vax_4 = discrete_data_generation(params,'max_vaccines',max_vaccines_options,3);
writematrix(max_vax_4, 'max_vax_4.csv')
end


%% Generate discrete data
function discrete_data = discrete_data_generation(params,variable_name,variable_options,options_num)
%preallocate data matrix
step = 250;

% Different steps when varying vaccine coverage
if params.max_vaccines == 2000
    vaccine_scenarios = [50:step:1950; 1950:-step:50];      
end 

if params.max_vaccines == 5000
    vaccine_scenarios = [50:step:4950; 4950:-step:50];
end 

if params.max_vaccines == 8000
    vaccine_scenarios = [3050:step:4950; 4950:-step:3050];
end

num_points = length(vaccine_scenarios);
%looping over discrete options
for j = 1:options_num   
    
    %need to adjust beta for contact matrix
    if strcmp(string(variable_name),'contact_matrix')
        [beta_temp,temp] = r0_beta_calibration(params,[0;0],params.R0,-1);
        params.beta = beta_temp;
    else 
        params.(string(variable_name)) = variable_options.(strcat('opt',string(j)));
    end
    
    
    if params.max_vaccines == 2000
        vaccine_scenarios = [50:step:1950; 1950:-step:50];      
    end 

    if params.max_vaccines == 5000
        vaccine_scenarios = [50:step:4950; 4950:-step:50];
    end 

    if params.max_vaccines == 8000
        vaccine_scenarios = [3050:step:4950; 4950:-step:3050];
    end

    num_points = length(vaccine_scenarios);
    
    
    
    %total infections
    discrete_data(1:2,j) = optimisation(params,@total_infections);
    
    count_1 = 1;
    for i = 3:num_points
       discrete_data(i,j) = total_infections(params,vaccine_scenarios(1,count_1));
       count_1 = count_1+1;
    end

    %total symp_infections
    discrete_data(num_points + 1:num_points + 2,j) = optimisation(params,@total_symp_infections);
    count_2 = 1;
    for i = num_points + 3:2*num_points
       discrete_data(i,j) = total_symp_infections(params,vaccine_scenarios(1,count_2));
       count_2 = count_2+1;
    end
    
    %total deaths
    discrete_data(2*num_points + 1: 2*num_points + 2,j) = optimisation(params,@total_deaths);
    count_3 = 1;
    for i = 2*num_points + 3: 3*num_points
       discrete_data(i,j) = total_deaths(params,vaccine_scenarios(1,count_3));
       count_3 = count_3 + 1;
    end
end

end



function cont_data = continuous_data_generation(params,size,variable_name,variable_options)

%preallocate data matrix
cont_data = zeros(16,size);

%iterate over options
for i = 1:size
    %change variable
    params.(string(variable_name)) = variable_options(i);
    
    %only recompute beta for R0
    if strcmp('R0',string(variable_name))
        [beta_temp,temp] = r0_beta_calibration(params,[0;0],params.R0,-1);
        params.beta = beta_temp;
    end

    %1st row: variable value
    cont_data(1,i) = variable_options(i);
    
    %2nd and 3rd rows: optimal strategy total_infections
    optimal_vax = optimisation(params, @total_infections);
    cont_data(2:3,i) = optimal_vax';
    
    %4th, 5th and 6th rows: total_infections for optimal, average, worst strats
    optimal_inf = total_infections(params,optimal_vax(1));
    worst_strat = flip(optimal_vax,2);
    worst_inf = total_infections(params,worst_strat(1));
    average_inf = total_infections(params,params.max_vaccines/2);
    cont_data(4:6,i) = [optimal_inf; worst_inf; average_inf];
    
    %7th and 8th rows: optimal strategy total_symp_infections
    optimal_vax = optimisation(params, @total_symp_infections);
    cont_data(7:8,i) = optimal_vax';
    
    %9th, 10th and 11th rows: total_symp_infections for optimal, average, worst strats
    optimal_inf = total_symp_infections(params,optimal_vax(1));
    worst_strat = flip(optimal_vax,2);
    worst_inf = total_symp_infections(params,worst_strat(1));
    average_inf = total_symp_infections(params,params.max_vaccines/2);
    cont_data(9:11,i) = [optimal_inf; worst_inf; average_inf];
    
    %12th and 13th rows: optimal strategy total_deaths
    optimal_vax = optimisation(params, @total_deaths);
    cont_data(12:13,i) = optimal_vax';
    
    %14th, 15th and 16th rows: total_deaths for optimal, average, worst strats
    optimal_inf = total_deaths(params,optimal_vax(1));
    worst_strat = flip(optimal_vax,2);
    worst_inf = total_deaths(params,worst_strat(1));
    average_inf = total_deaths(params,params.max_vaccines/2);
    cont_data(14:16,i) = [optimal_inf; worst_inf; average_inf];
    
    
end

end 

%% Optimisation 
function optimal_strat = optimisation(params,objective_function)
%defining optimisation conditions
LB = zeros(params.dim,1);
UB = params.pop_size - params.I0;
A = ones(1,params.dim -1);

%starting points: all vaccines allocated in each compartment and average
start = zeros(params.dim +1, params.dim);
average = floor(params.max_vaccines/params.dim);
start(1,1:params.dim) = average*ones(1,params.dim);
start(2:end,1:params.dim) = params.max_vaccines*eye(params.dim);

%iterate over different starting allocations 
vax_allocations = zeros(params.dim + 1, params.dim);
for i = (1:params.dim+1)
        vax_allocations(i,1:params.dim-1) = fmincon(@(x) ...
            objective_function(params,x),...
            start(i,1:end-1)',A,params.max_vaccines,[],[],LB,UB);
end

vax_allocations(:,end) = params.max_vaccines - sum(vax_allocations,2);

% finding overall optimal strategy
vax_allocations_objective = zeros(params.dim+1,1);
for i = (1:params.dim+1)
        vax_allocations_objective(i,1) = objective_function(params,vax_allocations(i,1:params.dim-1));
end

[min_infections, min_index] = min(vax_allocations_objective);

optimal_strat = vax_allocations(min_index,:);

end


%% Objective functions
function total_i = total_infections(params,vax)
all_vax = [vax;params.max_vaccines - vax];
[temp_tout, temp_yout, total_i,total_i_1,total_symp_i,total_symp_i_1,total_d,total_d_1] = model(params,all_vax);
end

function total_i_1 = total_infections_group_1(params,vax)
all_vax = [vax;params.max_vaccines - vax];
[temp_tout, temp_yout, total_i,total_i_1,total_symp_i,total_symp_i_1,total_d,total_d_1] = model(params,all_vax);
end

function total_symp_i = total_symp_infections(params,vax)
all_vax = [vax;params.max_vaccines - vax];
[temp_tout, temp_yout, total_i,total_i_1,total_symp_i,total_symp_i_1,total_d,total_d_1] = model(params,all_vax);
end

function total_symp_i_1 = total_symp_infections_group_1(params,vax)
all_vax = [vax;params.max_vaccines - vax];
[temp_tout, temp_yout, total_i,total_i_1,total_symp_i,total_symp_i_1,total_d,total_d_1] = model(params,all_vax);
end

function total_d = total_deaths(params,vax)
all_vax = [vax;params.max_vaccines - vax];
[temp_tout, temp_yout,total_i,total_i_1,total_symp_i,total_symp_i_1,total_d,total_d_1] = model(params,all_vax);
end

function total_d = total_deaths_1(params,vax)
all_vax = [vax;params.max_vaccines - vax];
[temp_tout, temp_yout,total_i,total_i_1,total_symp_i,total_symp_i_1,total_d,total_d_1] = model(params,all_vax);
end

%% Modelling
function [tout,yout,total_i,total_i_1,total_symp_i,total_symp_i_1,total_d,total_d_1] = ...
    model(params,vax)

%time
T = 50;

%initial condition
%assumptions: initial infected people split evenly between symptomatic and
%               asymptomatic 
%[S1;S2;V1;V2;ES1;ES2;EV1;EV2;ISS1;ISS2;ISA1;ISA2;IVS1;IVS2;IVA1;IVA2;
%         RSS1;RSS2;DSS1;DSS2;RSA1;RVS1;RVS2;DVS1;DVS2;RVA1;RSA2;RVA2]
initial = [params.pop_size - vax - params.I0; vax; zeros(2*params.dim,1);...
    0.5*params.I0; 0.5*params.I0; zeros(8*params.dim,1)];

%ensuring against overallocating vaccines
initial(initial<0) = 0;

%simulating model 
[tout,yout] = ode45(@deriv,[0,T],initial,[],params);

%calculating objective functions
RSS1 = 17; RSS2 = 18;
DSS1 = 19; DSS2 = 20;
RVS1 = 21; RVS2 = 22;
DVS1 = 23; DVS2 = 24;
RSA1 = 25; RSA2 = 26;
RVA1 = 27; RVA2 = 28;

total_i = sum(yout(end,RSS1:RVA2));
total_i_1 = yout(end,RSS1) + yout(end,DSS1) + yout(end,RVS1) +...
            yout(end,DVS1) + yout(end,RSA1) + yout(end,RVA1);
total_symp_i = sum(yout(end,RSS1:DVS2));
total_symp_i_1 = yout(end,RSS1) + yout(end,DSS1) + yout(end,RVS1)+...
            yout(end,DVS1);
total_d = sum(yout(end,DSS1:DSS2)) + sum(yout(end,DVS1:DVS2));
total_d_1 = yout(end,DSS1) + yout(end,DVS1);


function dYdt = deriv(t,Y,params)
    dim = params.dim;
    
    S = Y(1:dim,1);
    V = Y((dim + 1):2*dim, 1);
    ES = Y((2*dim + 1):3*dim,1);
    EV = Y((3*dim + 1):4*dim,1);
    ISS = Y((4*dim + 1):5*dim,1);
    ISA = Y((5*dim + 1):6*dim,1);
    IVS = Y((6*dim + 1):7*dim,1);
    IVA = Y((7*dim + 1):8*dim,1);
    RSS = Y((8*dim + 1):9*dim,1);
    DSS = Y((9*dim + 1):10*dim,1);
    RVS = Y((10*dim + 1):11*dim,1);
    DVS = Y((11*dim + 1):12*dim,1);
    RSA = Y((12*dim + 1):13*dim,1);
    RVA = Y((13*dim + 1):14*dim,1);
    
    LP = params.pop_size - DSS - DVS;
    
    inverse_pop_size = LP.^-1;
    beta = params.beta * params.contact_matrix * diag(inverse_pop_size);
  
    lambda = beta * diag(ISS) +...
        params.delta_A* beta * diag(ISA)+...
        (1-params.alpha_4) * (beta * diag(IVS) + ...
        params.delta_A * beta * diag(IVA));
    
    lambdaS = diag(lambda * [S';S']);
    lambdaV = diag(lambda * [V';V']);
    
    dSdt = - lambdaS;
    dVdt = - (1-params.alpha_1) * lambdaV;
    dESdt = lambdaS - params.theta.*ES;
    dEVdt = (1-params.alpha_1) * lambdaV - params.theta.*EV;
    dISSdt = params.omega*params.theta.*ES - params.gamma.*ISS;
    dISAdt = (1-params.omega)*params.theta.*ES - params.gamma.*ISA;
    dIVSdt = (1-params.alpha_2)*params.omega*params.theta.*EV...
        - params.gamma.*IVS;
    dIVAdt = (1-(1-params.alpha_2)*params.omega)*params.theta.*EV...
        - params.gamma.*IVA;
    dRSSdt = (1-params.sigma)*params.gamma.*ISS;
    dDSSdt = params.sigma.*params.gamma.*ISS;
    dRVSdt = (1-(1-params.alpha_3)*params.sigma)*params.gamma.*IVS;
    dDVSdt = (1-params.alpha_3)*params.sigma*params.gamma.*IVS;
    dRSAdt = params.gamma * ISA;
    dRVAdt = params.gamma * IVA;
    
    dYdt = [dSdt; dVdt; dESdt; dEVdt; dISSdt; dISAdt; dIVSdt; dIVAdt;...
        dRSSdt; dDSSdt; dRVSdt; dDVSdt; dRSAdt; dRVAdt];
    
end
end

%% Finding optimal vaccination strategy 
function [optimal_vax_i, optimal_vax_symp_i, optimal_vax_d] = vax_optimisation(params)

%defining optimisation conditions
LB = 0;
UB = params.pop_size(1) - params.I0(1);
A = 1;

%multiple starting points 
average = floor(params.max_vaccines/params.dim);
start = [params.max_vaccines; average; 0];

%iterate over different starting allocations
vax_allocations_i = zeros(params.dim + 1, params.dim);
vax_allocations_symp_i = zeros(params.dim + 1, params.dim);
vax_allocations_d = zeros(params.dim + 1, params.dim);
for i = (1:params.dim+1)
    %objective function: total infections
    vax_allocations_i(i,1) = fmincon(@(x) total_infections(params,x),...
        start(i,1)',A,params.max_vaccines,[],[],LB,UB);
    
    %objective function: total symptomatic infections
    vax_allocations_symp_i(i,1) = fmincon(@(x) total_symp_infections(params,x),...
        start(i,1),A,params.max_vaccines,[],[],LB,UB);
    
    %objective function: total deaths
    vax_allocations_d(i,1) = fmincon(@(x) total_deaths(params,x),...
        start(i,1)',A,params.max_vaccines,[],[],LB,UB);
    
end

%allocating vaccines in second group
vax_allocations_i(:,params.dim) = params.max_vaccines - vax_allocations_i(:,1);
vax_allocations_symp_i(:,params.dim) = params.max_vaccines - vax_allocations_symp_i(:,1);
vax_allocations_d(:,params.dim) = params.max_vaccines - vax_allocations_d(:,1);

%finding optimal strategy
objective_i = zeros(params.dim+1,1);
objective_symp_i = zeros(params.dim+1,1);
objective_d = zeros(params.dim+1,1);

for i = (1:params.dim+1)
    objective_i(i,1) = total_infections(params,vax_allocations_i(i,1));
    objective_symp_i(i,1) = total_symp_infections(params,vax_allocations_symp_i(i,1));
    objective_d(i,1) = total_deaths(params,vax_allocations_d(i,1));
end

[temp, min_index_i] = min(objective_i);
[temp, min_index_symp_i] = min(objective_symp_i);
[temp, min_index_d] = min(objective_d);


optimal_vax_i = vax_allocations_i(min_index_i,:)
optimal_vax_symp_i = vax_allocations_symp_i(min_index_symp_i,:)
optimal_vax_d = vax_allocations_d(min_index_d,:)

end






