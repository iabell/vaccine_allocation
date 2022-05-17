function [beta,R0] = r0_beta_calibration(params,v,R0_in,beta_in)
%R0_in: R0 for calculating beta (input -1 if calculating R0)
%beta_in: beta for calculating R0 (input -1 if calculating beta)
N1 = params.pop_size(1);
N2 = params.pop_size(2);

%Disease free equilibrium (DFE)
v1 = v(1);
v2 = v(2);
S1 = N1 - v1;
S2 = N2 - v2;

%Living population (frequency dependent transmission)
LP1 = N1;
LP2 = N2;


%Contact matrix values
B_11 = @(b) b*params.contact_matrix(1,1)/LP1;
B_12 = @(b) b*params.contact_matrix(1,2)/LP2;
B_21 = @(b) b*params.contact_matrix(2,1)/LP1;
B_22 = @(b) b*params.contact_matrix(2,2)/LP2;

%Transmission matrix
F_left = zeros(4,4);
F_right_row_1 = @(b) [params.delta_A*B_11(b)*S1, B_11(b)*S1...
    params.alpha_4*params.delta_A*B_11(b)*S1, params.alpha_4*B_11(b)*S1...
    params.delta_A*B_12(b)*S1, B_12(b)*S1...
    params.alpha_4*params.delta_A*B_12(b)*S1, params.alpha_4*B_12(b)*S1];

F_right_row_2 = @(b) params.alpha_1 * F_right_row_1(b);

F_right_row_3 = @(b) [params.delta_A*B_21(b)*S2, B_21(b)*S2...
    params.alpha_4*params.delta_A*B_21(b)*S2, params.alpha_4*B_21(b)*S2...
    params.delta_A*B_22(b)*S2, B_22(b)*S2...
    params.alpha_4*params.delta_A*B_22(b)*S2, params.alpha_4*B_22(b)*S2];

F_right_row_4 = @(b) params.alpha_1*F_right_row_3(b);

F_right = @(b) [F_right_row_1(b); F_right_row_2(b); F_right_row_3(b);... 
    F_right_row_4(b)];

F_bottom = zeros(8,12);

F_mat = @(b) [F_left, F_right(b); F_bottom];

%Transition matrix
V_left = params.theta*eye(4);
V_right = zeros(4,8);

V_bottom_left_diag_vect = -[(1-params.omega(1))*params.theta;... 
    params.omega(1)*params.theta;...
    (1-params.alpha_2*params.omega(1))*params.theta;...
    params.alpha_2*params.omega(1)*params.theta;...
    (1-params.omega(2))*params.theta;... 
    params.omega(2)*params.theta;...
    (1-params.alpha_2*params.omega(2))*params.theta;...
    params.alpha_2*params.omega(2)*params.theta];

V_bottom_left = [diag(V_bottom_left_diag_vect) zeros(8,4)];

V_bottom_right = [zeros(8,4) params.gamma*eye(8)];

V_bottom = V_bottom_left + V_bottom_right;

V_mat = [V_left V_right; V_bottom];

V_mat_inverse = inv(V_mat);

%Next generation matrix
NGM = @(b) F_mat(b) * V_mat_inverse;
    
    
%Calculating beta given R0
if beta_in<0
    R0_func = @(b) eigs(NGM(b),1) - R0_in;
    beta = fzero(R0_func,0);
    R0 = R0_in;
end
    
%Calculating R0 given beta
if R0_in <0
    R0 = eigs(NGM(beta_in),1);
    beta = beta_in;    
end

end