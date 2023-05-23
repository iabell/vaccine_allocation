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

% R0 sensitivity analysis
R0_sensitivity_analysis(params)

disease_example_data = true;
omega_data = true; 
sigma_data = true; 
delta_A_data = true;
R0_data = true;


if disease_example_data
    illustrative_example(params) 
end

if omega_data 
    params.alpha_1 = 0;
    params.alpha_2 = 0;
    params.alpha_3 = 0;
    params.alpha_4 = 0;
    
    %vaccine 1
    params.alpha_1 = 0.75;
    omega_1 = continuous_data_generation(params,size,'omega',[linspace(0,1,size)', 0.5*ones(size,1)]');
    writematrix(omega_1,'omega_data_1_cont.csv')

    
    %vaccine 2
    params.alpha_1 = 0;
    params.alpha_2 = 0.75;
    omega_2 = continuous_data_generation(params,size,'omega',[linspace(0,1,size)', 0.5*ones(size,1)]');
    writematrix(omega_2,'omega_data_2_cont.csv')

    %vaccine 3
    params.alpha_2 = 0;
    params.alpha_3 = 0.75;
    omega_3 = continuous_data_generation(params,size,'omega',[linspace(0,1,size)', 0.5*ones(size,1)]');
    writematrix(omega_3,'omega_data_3_cont.csv')

    %vaccine 4
    params.alpha_3 = 0;
    params.alpha_4 = 0.75;
    omega_4 = continuous_data_generation(params,size,'omega',[linspace(0,1,size)', 0.5*ones(size,1)]');
    writematrix(omega_4,'omega_data_4_cont.csv')
end 

if sigma_data
    params.alpha_1 = 0;
    params.alpha_2 = 0;
    params.alpha_3 = 0;
    params.alpha_4 = 0;
    
    %vaccine 1
    params.alpha_1 = 0.75;
    sigma_1 = continuous_data_generation(params,size,'sigma',[linspace(0,1,size)', 0.5*ones(size,1)]');
    writematrix(sigma_1,'sigma_data_1.csv')
     
    %vaccine 2
    params.alpha_1 = 0;
    params.alpha_2 = 0.75;
    sigma_2 = continuous_data_generation(params,size,'sigma',[linspace(0,1,size)', 0.5*ones(size,1)]');
    writematrix(sigma_2,'sigma_data_2.csv')

    %vaccine 3
    params.alpha_2 = 0;
    params.alpha_3 = 0.75;
    sigma_3 = continuous_data_generation(params,size,'sigma',[linspace(0,1,size)', 0.5*ones(size,1)]');
    writematrix(sigma_3,'sigma_data_3_75.csv')

    %vaccine 4
    params.alpha_3 = 0;
    params.alpha_4 = 0.75;
    sigma_4 = continuous_data_generation(params,size,'sigma',[linspace(0,1,size)', 0.5*ones(size,1)]');
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
total_population = sum(params.pop_size);
cont_max_vaccines_options = 0:100:total_population - sum(params.I0);
size_max_vaccines = length(cont_max_vaccines_options);

%vaccine 1
params.alpha_1 = 0;
params.alpha_2 = 0;
params.alpha_3 = 0;
params.alpha_4 = 0;
ve_1 = continuous_data_generation(params,size-1,'alpha_1',ve_options);
writematrix(ve_1, 've_1.csv')

params.alpha_1 = 0.75;
cont_max_vax_1 = continuous_data_generation(params,size_max_vaccines,'max_vaccines',cont_max_vaccines_options);
writematrix(cont_max_vax_1, 'cont_max_vax_1.csv')

%vaccine 2
params.alpha_1 = 0;
params.alpha_2 = 0;
params.alpha_3 = 0;
params.alpha_4 = 0;
ve_2 = continuous_data_generation(params,size-1,'alpha_2',ve_options);
writematrix(ve_2, 've_2.csv')

params.alpha_2 = 0.75;
cont_max_vax_2 = continuous_data_generation(params,size_max_vaccines,'max_vaccines',cont_max_vaccines_options);
writematrix(cont_max_vax_2, 'cont_max_vax_2.csv')

%vaccine 3
params.alpha_1 = 0;
params.alpha_2 = 0;
params.alpha_3 = 0;
params.alpha_4 = 0;
ve_3 = continuous_data_generation(params,size-1,'alpha_3',ve_options);
writematrix(ve_3, 've_3.csv')

params.alpha_3 = 0.75;
cont_max_vax_3 = continuous_data_generation(params,size_max_vaccines,'max_vaccines',cont_max_vaccines_options);
writematrix(cont_max_vax_3, 'cont_max_vax_3.csv')

%vaccine 4
params.alpha_1 = 0;
params.alpha_2 = 0;
params.alpha_3 = 0;
params.alpha_4 = 0;
ve_4 = continuous_data_generation(params,size-1,'alpha_4',ve_options);
writematrix(ve_4, 've_4.csv')

params.alpha_4 = 0.75;
cont_max_vax_4 = continuous_data_generation(params,size_max_vaccines,'max_vaccines',cont_max_vaccines_options);
writematrix(cont_max_vax_4, 'cont_max_vax_4.csv')
end


function cont_data = continuous_data_generation(params,size,variable_name,variable_options)

%preallocate data matrix
cont_data = zeros(16,size);

%iterate over options
for i = 1:size
    
    %change variable
    if strcmp('contact_matrix', string(variable_name))
        cmat = params.contact_matrix;
        params.contact_matrix = ...
        [variable_options(i), variable_options(i); 1, 1].*cmat;
        [beta_temp,temp] = r0_beta_calibration(params,[0;0],params.R0,-1);
        params.beta = beta_temp;
    
    else
        params.(string(variable_name)) = variable_options(:,i);

    end
    
    %recompute beta for R0
    if strcmp('R0',string(variable_name))
        [beta_temp,temp] = r0_beta_calibration(params,[0;0],params.R0,-1);
        params.beta = beta_temp;
    end

    

    %1st row: variable value
    cont_data(1,i) = variable_options(1,i);
    
    %2nd and 3rd rows: optimal strategy total_infections
    optimal_vax = optimisation(params, @total_infections);
    cont_data(2:3,i) = optimal_vax';
    
    %4th, 5th and 6th rows: total_infections for optimal, average, worst strats
    optimal_inf = total_infections(params,optimal_vax(1));
    worst_strat = flip(optimal_vax,2);
    worst_inf = total_infections(params,worst_strat(1));
    average_strat = params.max_vaccines/2;
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

function [total_i_1, total_i_2] = total_infections_by_group(params,vax)
all_vax = [vax;params.max_vaccines - vax];
[temp_tout, temp_yout, total_i,total_i_1,total_symp_i,total_symp_i_1,total_d,total_d_1] = model(params,all_vax);
total_i_2 = total_i - total_i_1;
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


function R0_sensitivity_analysis(params) 
    params.alpha_1 = 0.75;
    params.alpha_2 = 0;
    params.alpha_3 = 0;
    params.alpha_4 = 0;
    
    %varying R0
    R0_options = 1:0.2:10;
    
    %sensitivity analysis model 1
    vax_allocations = zeros(2,length(R0_options));
    total_infections_group_1 = zeros(1, length(R0_options));
    total_infections_group_2 = zeros(1, length(R0_options));
    total_infections = zeros(1, length(R0_options));
    total_infections_vax_gp_1 = zeros(3, length(R0_options));
    total_infections_vax_gp_2 = zeros(3, length(R0_options));
    total_infections_no_vax = zeros(3, length(R0_options));
    count = 1;
    for i = R0_options
        params.max_vaccines = 5000;
        params.R0 = i;
        [beta_temp, temp] = r0_beta_calibration(params, [0;0], params.R0, -1);
        params.beta = beta_temp;
%         total infections from vaccinating group 1
        [total_1_vax_1, total_2_vax_1] = total_infections_by_group(params, 4950);      
        total_infections_vax_gp_1(1, count) = total_1_vax_1;
        total_infections_vax_gp_1(2, count) = total_2_vax_1;
        total_infections_vax_gp_1(3, count) = total_1_vax_1 + total_2_vax_1;

%         total infections from vaccinating group 2
        [total_1_vax_2, total_2_vax_2] = total_infections_by_group(params, 50);
        total_infections_vax_gp_2(1, count) = total_1_vax_2;
        total_infections_vax_gp_2(2, count) = total_2_vax_2;
        total_infections_vax_gp_2(3, count) = total_1_vax_2 + total_2_vax_2;

%         total infections from no vaccination
        params.max_vaccines = 0;
        [total_1_no_vax, total_2_no_vax] = total_infections_by_group(params, 0);     
        total_infections_no_vax(1,count) = total_1_no_vax;
        total_infections_no_vax(2,count) = total_2_no_vax;
        total_infections_no_vax(3,count) = total_1_no_vax + total_2_no_vax;
        count = count+1;
    end

    writematrix(R0_options, 'R0_sense_R0_vals.csv')
    writematrix(total_infections_vax_gp_1,'R0_sens_vaccinate_gp_1_infections.csv')
    writematrix(total_infections_vax_gp_2,'R0_sens_vaccinate_gp_2_infections.csv')
    writematrix(total_infections_no_vax, 'R0_sens_no_vax_infections.csv')

end



