clc;
clear;
close all;

%% === Load the predefined parameters ===
syms s;
syms k_tilde_11_bar k_tilde_12_bar k_tilde_13_bar k_tilde_14_bar k_tilde_15_bar k_tilde_16_bar k_tilde_21_bar k_tilde_22_bar k_tilde_23_bar k_tilde_24_bar k_tilde_25_bar k_tilde_26_bar k_tilde_31_bar k_tilde_32_bar k_tilde_33_bar k_tilde_34_bar k_tilde_35_bar k_tilde_36_bar;
output = parameter_init();
A = output{1};
B = output{2};
C = output{3};
x_0 = output{5};

%% === Controllability check ===
W_c = [B, A*B, A^2*B, A^3*B, A^4*B, A^5*B]; % Controllability matrix
r_Wc = rank(W_c); % Check if controllability matrix is full rank
if r_Wc == 6
    disp("The system is controllable.");
else
    fprintf("The rank of controllability matrix is %d. The system is NOT controllable.", r_Wc);
end

%% === Observability check ===
W_o = [C; C*A; C*A^2; C*A^3; C*A^4; C*A^5]; % Observability matrix
r_Wo = rank(W_o); % Check if obervability matrix is full rank
if r_Wo == 6
    disp("The system is observable.");
else
    fprintf("The rank of observability matrix is %d. The system is NOT observable.", r_Wo);
end

%% === LQR design and solve for the state feedback gain K ===
Q = diag([30, 20, 5, 10, 10, 10]); % State weighting matrix Q: penalize deviations in systems states
R = diag([5, 5]); % Control input weighting matrix R: penalize control energy

% Generalized Hamiltonian matrix Gamma used for solving the Algebraic Riccati Equation
Gamma = [A, -B/R*B';
        -Q, -A'];

[eigen_vector, eigen_value] = eig(Gamma);

% Find the eigenvectors associated with eigenvalues that have negative real parts
stable_idx = find(real(diag(eigen_value))<0);
vector_stable = eigen_vector(:, stable_idx);

n_state = size(A, 1);
v_vector = vector_stable(1:n_state, :); % Upper block v
mu_vector = vector_stable(n_state + 1:end, :); % Lower block mu
P = real(mu_vector/v_vector);
K = R\B'*P; % Get the optimal LQR feedback gain matrix K

%% === Controllability check ===
A_tilde = A';
B_tilde = C';

W_est_c = [B_tilde, A_tilde*B_tilde, A_tilde^2*B_tilde, A_tilde^3*B_tilde, A_tilde^4*B_tilde, A_tilde^5*B_tilde]; % Controllability matrix
r_W_est_c = rank(W_est_c); % Check if controllability matrix is full rank
if r_W_est_c == 6
    disp("The dual system is controllable.");
else
    fprintf("The rank of controllability matrix is %d. The dual system is NOT controllable.", r_W_est_c);
end

%%  === Define the structured controllability basis matrix and transformation matrix ===
C_matrix = [B_tilde(:,1), A_tilde*B_tilde(:,1), B_tilde(:,2), A_tilde*B_tilde(:,2), B_tilde(:,3), A_tilde*B_tilde(:,3)]; % Select 6 linearly-independent column vectors to form a new C square matrix
r_Cm = rank(C_matrix); % Check if C_matrix is full rank
if r_Cm == 6
    disp("The structured controllability basis matrix is legal.");
else
    fprintf("The rank of structured controllability basis matrix is %d. It's NOT legal.", r_Cm);
end

C_inv = inv(C_matrix);
d_1 = 2; 
d_2 = 2; 
d_3 = 2; % The number of vectors associated with the 3 second-order subsystems
q_1 = C_inv(d_1, :);
q_2 = C_inv(d_1 + d_2, :);
q_3 = C_inv(d_1 + d_2 + d_3, :); % q_1, q_2 and q_3 are rows taken from C_inv to form T
T = [q_1; q_1*A_tilde; q_2; q_2*A_tilde; q_3; q_3*A_tilde]; % Transformation matrix

%% === Solve for the state feedback gain L via symbolic pole placement in controllable canonical form ===
desired_poles = [-10, -10, -20, -20, -30, -30];
poly_1 = (s - desired_poles(1))*(s - desired_poles(2)); % 2 poles for the first subsystem
poly_2 = (s - desired_poles(3))*(s - desired_poles(4)); % 2 poles for the second subsystem
poly_3 = (s - desired_poles(5))*(s - desired_poles(6)); % 2 poles for the third subsystem
coef_1 = double(coeffs(poly_1)); % Get the coefficient [a2,a1,1]
coef_2 = double(coeffs(poly_2)); % Get the coefficient [a4,a3,1]
coef_3 = double(coeffs(poly_3)); % Get the coefficient [a6,a5,1]

% Transform the system from its original coordinates to the controllable canonical form coordinates
A_tilde_bar = T*A_tilde/T;
B_tilde_bar = T*B_tilde;

% Controllable canonical form gain matrix K_tilde_bar  
K_tilde_bar = [k_tilde_11_bar, k_tilde_12_bar, k_tilde_13_bar, k_tilde_14_bar, k_tilde_15_bar, k_tilde_16_bar;
               k_tilde_21_bar, k_tilde_22_bar, k_tilde_23_bar, k_tilde_24_bar, k_tilde_25_bar, k_tilde_26_bar;
               k_tilde_31_bar, k_tilde_32_bar, k_tilde_33_bar, k_tilde_34_bar, k_tilde_35_bar, k_tilde_36_bar];

A_tilde_d = A_tilde_bar - B_tilde_bar*K_tilde_bar; % Controllable canonical form state matrix A_tilde_d

target = [ A_tilde_d(2,:) == -[coef_1(1:2),0, 0, 0, 0], A_tilde_d(4,:) == -[0, 0, coef_2(1:2), 0, 0], A_tilde_d(6,:) == -[0, 0, 0, 0, coef_3(1:2)] ];
vars = [k_tilde_11_bar k_tilde_12_bar k_tilde_13_bar k_tilde_14_bar k_tilde_15_bar k_tilde_16_bar k_tilde_21_bar k_tilde_22_bar k_tilde_23_bar k_tilde_24_bar k_tilde_25_bar k_tilde_26_bar k_tilde_31_bar k_tilde_32_bar k_tilde_33_bar k_tilde_34_bar k_tilde_35_bar k_tilde_36_bar];
sol = solve(target, vars);
K_tilde_bar_sol = double([sol.k_tilde_11_bar sol.k_tilde_12_bar sol.k_tilde_13_bar sol.k_tilde_14_bar sol.k_tilde_15_bar sol.k_tilde_16_bar;
                          sol.k_tilde_21_bar sol.k_tilde_22_bar sol.k_tilde_23_bar sol.k_tilde_24_bar sol.k_tilde_25_bar sol.k_tilde_26_bar;
                          sol.k_tilde_31_bar sol.k_tilde_32_bar sol.k_tilde_33_bar sol.k_tilde_34_bar sol.k_tilde_35_bar sol.k_tilde_36_bar;]);

K_tilde = K_tilde_bar_sol*T; % The original gain matrix K_tilde
L = K_tilde'; % Get the observer gain