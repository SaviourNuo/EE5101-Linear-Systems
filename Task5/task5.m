clc;
clear;
close all;

%% === Load the predefined parameters ===
output = parameter_init();
A = output{1};
B = output{2};
C = output{3};
y_sp = output{4};
x_0 = output{5};
w = [-1;1]; % Disturbance (step disturbance at 10s)

%% === Augmented system design ===
% Preparations for the augamented system
A_bar = [A, zeros(6, 3);
        -C, zeros(3, 3)];
B_bar = [B; zeros(3,2)];
B_wbar = [B; zeros(3,2)];
B_rbar = [zeros(6,3); eye(3)];
C_bar = [C, zeros(3,3)];

%% === Controllability check ===
% Check controllability of original system
W_c = [B, A*B, A^2*B, A^3*B, A^4*B, A^5*B]; % Controllability matrix
r_Wc = rank(W_c); % Check if controllability matrix is full rank
if r_Wc == 6
    disp("The system is controllable.");
else
    fprintf("The rank of controllability matrix is %d. The system is NOT controllable.", r_Wc);
end

% Check controllability of the augmented system
W_qc = [A, B;
        C, zeros(3,2)];
r_Wqc = rank(W_qc); % Check if controllability matrix is full rank
if r_Wqc == 8
    disp("The augmented system is controllable.");
else
    fprintf("The rank of controllability matrix is %d. The augmented system is NOT controllable.", r_Wqc);
end

%% === LQR design and solve for the state feedback gain K ===
Q = diag([30, 20, 5, 10, 10, 10, 10, 60, 80]); % State weighting matrix Q: penalize deviations in systems states
R = diag([5, 5]); % Control input weighting matrix R: penalize control energy

% Generalized Hamiltonian matrix Gamma used for solving the Algebraic Riccati Equation
Gamma = [A_bar, -B_bar/R*B_bar';
        -Q, -A_bar'];

[eigen_vector, eigen_value] = eig(Gamma);

% Find the eigenvectors associated with eigenvalues that have negative real parts
stable_idx = find(real(diag(eigen_value))<0);
vector_stable = eigen_vector(:, stable_idx);

n_state = size(A_bar, 1);
v_vector = vector_stable(1:n_state, :); % Upper block v
mu_vector = vector_stable(n_state + 1:end, :); % Lower block mu
P = real(mu_vector/v_vector);
K = R\B_bar'*P; % Get the optimal LQR feedback gain matrix K

%% === LQR design and solve for the state feedback gain L ===
A_tilde = A';
B_tilde = C';
Q_tilde = eye(6);
R_tilde = eye(3);

% Generalized Hamiltonian matrix Gamma used for solving the Algebraic Riccati Equation
Gamma_tilde = [A_tilde, -B_tilde/R_tilde*B_tilde';
        -Q_tilde, -A_tilde'];

[eigen_vector_tilde, eigen_value_tilde] = eig(Gamma_tilde);

% Find the eigenvectors associated with eigenvalues that have negative real parts
stable_idx_tilde = find(real(diag(eigen_value_tilde))<0);
vector_stable_tilde = eigen_vector_tilde(:, stable_idx_tilde);

n_state_tilde = size(A, 1);
v_vector_tilde = vector_stable_tilde(1:n_state_tilde, :); % Upper block v
mu_vector_tilde = vector_stable_tilde(n_state_tilde + 1:end, :); % Lower block mu
P_tilde = real(mu_vector_tilde/v_vector_tilde);
K_tilde = R_tilde\B_tilde'*P_tilde; % Get the optimal LQR feedback gain matrix K_tilde
L = K_tilde';

disp(L);