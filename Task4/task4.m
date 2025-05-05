clc;
clear;
close all;

%% === Load the predefined parameters ===
syms s;
output = parameter_init();
A = output{1};
B = output{2};
C = output{3};
x_0 = output{5};
C_2 = [1, 0, 0, 0, 0, 0;
       0, 0, 1, 0, 0, 0]; % Only 2 ouputs are considered

%% === Controllability check ===
W_c = [B, A*B, A^2*B, A^3*B, A^4*B, A^5*B]; % Controllability matrix
r_Wc = rank(W_c); % Check if controllability matrix is full rank
if r_Wc == 6
    disp("The original system is controllable.");
else
    fprintf("The rank of controllability matrix is %d. The original system is NOT controllable.", r_Wc);
end

%% === Find the relative degree ===
degree = zeros(2, 1); % Initialize the relative degree for each output channel
for k = 1:2 % Find the smallest i that makes C_2(k,:)*A^(i-1)*B not zero
    degree(k) = find(arrayfun(@(i) norm(C_2(k,:)*A^(i-1)*B, 'fro') > 1e-10, 1:6), 1);
end

degree_1 = degree(1); % First relative degree sigma_1
degree_2 = degree(2); % Second relative degree sigma_1

B_star = [C_2(1,:)*A^(degree_1-1)*B;C_2(2,:)*A^(degree_2-1)*B]; % B* for the decoupled system

phi_poles_1 = [-6, -8]; % Pole placement for the first output channel
phi_poles_2 = [-10, -12]; % Pole placement for the second output channel

I = eye(size(A, 1)); % Identity matrix I
phi_f1 = (A - phi_poles_1(1)*I)*(A - phi_poles_1(2)*I);
phi_f2 = (A - phi_poles_2(1)*I)*(A - phi_poles_2(2)*I); % Stable characteristic polynomial

C_star = [C_2(1,:)*phi_f1; C_2(2,:)*phi_f2]; % C* for the decoupled system
F = inv(B_star); % Feedforward compensation matrix F
K = B_star\C_star; % Feedback gain matrix K

%% === Step response (zero initial condition, step input) ===
Af = A - B * K; % Closed-loop matrix Af
Bf = B * F; % Input matrix Bf
Cf = C_2; % Output matrix Cf
sys_decouple = ss(Af, Bf, Cf, 0); % Create state-space system

t = 0:0.01:10; % Time span
u_0 = zeros(length(t), 2); % Initial input set to zero for free response
u_1 = [ones(length(t),1), zeros(length(t),1)]; % r1=[1,0] step function
u_2 = [zeros(length(t),1), ones(length(t),1)]; % r2=[0,1] step function

[y1, t_out1, ~] = lsim(sys_decouple, u_1, t);
figure;
plot(t_out1, y1, 'LineWidth', 1.2);
title('Step Response to r = [1; 0]');
xlabel('Time (s)');
ylabel('Output');
legend('d(t)', '\psi(t)');
grid on;

[y2, t_out2, ~] = lsim(sys_decouple, u_2, t);
figure;
plot(t_out2, y2, 'LineWidth', 1.2);
title('Step Response to r = [0; 1]');
xlabel('Time (s)');
ylabel('Output');
legend('d(t)', '\psi(t)');
grid on;

%% === Free response (non-zero initial condition, zero input) ===
% Simulate the free reponse of the system with initial state as x_0
[y_free, t_free, x_free] = lsim(sys_decouple, u_0, t, x_0);
u_free = -K * x_free';
u_free = u_free';

figure;
plot(t_free, x_free, 'LineWidth', 1.2);
title('State Response: x(0) ≠ 0, r = 0');
xlabel('Time (s)');
ylabel('States');
legend('x_1 = d(t)', 'x_2 = φ(t)', 'x_3 = ψ(t)', 'x_4 = d''(t)', 'x_5 = φ''(t)', 'x_6 = ψ''(t)', 'Location', 'southwest');
grid on;

figure;
plot(t_free, y_free, 'LineWidth', 1.2);
title('Output Response: x(0) ≠ 0, r = 0');
xlabel('Time (s)');
ylabel('Output');
legend('d(t)', '\psi(t)');
grid on;

figure;
plot(t_free, u_free, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Control Input u(t)');
title('Control Inputs during Free Response');
legend('u_c(t)', 'u_h(t)');
grid on;