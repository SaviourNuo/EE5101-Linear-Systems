clc;
clear;
close all;

%% === Load the predefined parameters ===
syms s;
syms k_11_bar k_12_bar k_13_bar k_14_bar k_15_bar k_16_bar k_21_bar k_22_bar k_23_bar k_24_bar k_25_bar k_26_bar;
output = parameter_init();
A = output{1};
B = output{2};
C = output{3};
x_0 = output{5};

%% === Design desired poles ===
zeta = 0.8;
omega = 1.5; % (ζ = 0.8, ωn = 1.5)
p_dominant = [-zeta*omega + omega*sqrt(1-zeta^2)*1j, -zeta*omega - omega*sqrt(1-zeta^2)*1j]; % Dominant poles
p_extra = [-4, -4, -6, -6]; % Extra non-dominant poles (2~5 times for faster convergence)
desired_poles = [p_dominant, p_extra]; % Get all the poles
% K_ideal = place(A, B, desired_poles); % Get the ideal feedback gain matrix K to place the desired poles on closed-loop system

max_ost = 0.1; % The overshoot is less than 10%
max_st = 5; % The 2% settling time is less than 5 seconds
M_p = exp(-pi*zeta/sqrt(1-zeta^2)); % Ideal overshoot calculation
t_s = 4/(zeta*omega); % Ideal settling time calculation

fprintf("According to the ideal second-order system formula\nThe overshoot is %.2f%%.\nThe 2%% settling time is %.3f seconds.\n", M_p*100, t_s);

%% === Controllability check ===
W_c = [B, A*B, A^2*B, A^3*B, A^4*B, A^5*B]; % Controllability matrix
r_Wc = rank(W_c); % Check if controllability matrix is full rank
if r_Wc == 6
    disp("The system is controllable.");
else
    fprintf("The rank of controllability matrix is %d. The system is NOT controllable.", r_Wc);
end

%%  === Define the structured controllability basis matrix and transformation matrix ===
C_matrix = [B(:,1), A*B(:,1), A^2*B(:,1), B(:,2), A*B(:,2), A^2*B(:,2)]; % Select 6 linearly-independent column vectors to form a new C square matrix
r_Cm = rank(C_matrix); % Check if C_matrix is full rank
if r_Cm == 6
    disp("The structured controllability basis matrix is legal.");
else
    fprintf("The rank of structured controllability basis matrix is %d. It's NOT legal.", r_Cm);
end

C_inv = inv(C_matrix);
d_1 = 3; % The number of vectors associated with u_c
d_2 = 3; % The number of vectors associated with u_h
q_1 = C_inv(d_1, :);
q_2 = C_inv(d_1 + d_2, :); % q_1 and q_2 are rows taken from C_inv to form T
T = [q_1; q_1*A; q_1*A^2; q_2; q_2*A; q_2*A^2;]; % Transformation matrix

%% === Solve for the state feedback gain K via symbolic pole placement in controllable canonical form ===
% Construct the characteristic polynomial of the subsystem
poly_1 = (s - desired_poles(1))*(s - desired_poles(2))*(s - desired_poles(3)); % 3 poles for the first subsystem
poly_2 = (s - desired_poles(4))*(s - desired_poles(5))*(s - desired_poles(6)); % 3 poles for the second subsystem
coef_1 = double(coeffs(poly_1)); % Get the coefficient [a3,a2,a1,1]
coef_2 = double(coeffs(poly_2)); % Get the coefficient [a6,a5,a4,1]

% Transform the system from its original coordinates to the controllable canonical form coordinates
A_bar = T*A/T;
B_bar = T*B;

% Controllable canonical form gain matrix K_bar  
K_bar=[k_11_bar, k_12_bar, k_13_bar, k_14_bar, k_15_bar, k_16_bar; % Affect u_c
       k_21_bar, k_22_bar, k_23_bar, k_24_bar, k_25_bar, k_26_bar]; % Affect u_h
A_d = A_bar - B_bar*K_bar; % Controllable canonical form state matrix A_d

% Match the 3rd and 6th rows of A_d to the desired controllable form 
target = [A_d(3,:) == -[coef_1(1:3), 0, 0, 0], A_d(6,:) == -[0, 0, 0, coef_2(1:3)]];
vars = [k_11_bar k_12_bar k_13_bar k_14_bar k_15_bar k_16_bar k_21_bar k_22_bar k_23_bar k_24_bar k_25_bar k_26_bar];
sol = solve(target, vars);
K_bar_sol = double([sol.k_11_bar sol.k_12_bar sol.k_13_bar sol.k_14_bar sol.k_15_bar sol.k_16_bar;
                    sol.k_21_bar sol.k_22_bar sol.k_23_bar sol.k_24_bar sol.k_25_bar sol.k_26_bar]);

K = K_bar_sol*T; % The original gain matrix K

%% === Step response (zero initial condition, step input) ===
A_cl = A - B*K; % The orignal closed-loop matrix A_cl
sys = ss(A_cl, B, C, 0); % Construct a state-space system

t = 0:0.01:10; % Time span
x0 = zeros(6,1); % Initial state set to zero for step response 
u_0 = zeros(length(t), 2); % Initial input set to zero for free response
u_1 = [ones(length(t),1), zeros(length(t),1)]; % r1=[1,0] step function
u_2 = [zeros(length(t),1), ones(length(t),1)]; % r2=[0,1] step function

% Simulate the reponse of the system to u_1 and u_2
[y1, t_out1] = lsim(sys, u_1, t, x0);
[y2, t_out2] = lsim(sys, u_2, t, x0);
figure; plot(t_out1, y1, 'LineWidth', 1.2); title('Step Response: [1 0] Input'); xlabel('Time (s)'); ylabel('Output y'); grid on; legend('d(t)', '\phi(t)', '\psi(t)');
figure; plot(t_out2, y2, 'LineWidth', 1.2); title('Step Response: [0 1] Input'); xlabel('Time (s)'); ylabel('Output y'); grid on; legend('d(t)', '\phi(t)', '\psi(t)');

%% === Free response (non-zero initial condition, zero input) ===
% Simulate the free reponse of the system with initial state as x_0
[y_free, t_free, x_free] = lsim(sys, u_0, t, x_0);
u_free = -K * x_free';
u_free = u_free';

figure;
plot(t_free, x_free, 'LineWidth', 1.3);
xlabel('Time (s)');
ylabel('States');
title('Free Response with x_0 ≠ 0 and u = 0');
legend('x_1 = d(t)', 'x_2 = φ(t)', 'x_3 = ψ(t)', 'x_4 = d''(t)', 'x_5 = φ''(t)', 'x_6 = ψ''(t)');
grid on;

figure;
plot(t_free, y_free, 'LineWidth', 1.3);
xlabel('Time (s)');
ylabel('Output y = Cx');
title('Output Response with x_0 ≠ 0 and u = 0');
legend('y_1 = d(t)', 'y_2 = φ(t)', 'y_3 = ψ(t)');
grid on;

figure;
plot(t_free, u_free, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Control Input u(t)');
title('Control Inputs during Free Response');
legend('u_c(t)', 'u_h(t)');
grid on;

fprintf('Max control input during free response:\n');
fprintf('u_c: %.2f\n', max(abs(u_free(:,1))));
fprintf('u_h: %.2f\n', max(abs(u_free(:,2))));

%% === Design specifications check ===
evaluate_response(y1, t_out1, '[1, 0]');
evaluate_response(y2, t_out2, '[0, 1]');

function evaluate_response(y, t, label)
    fprintf('\n=== Performance Report for Input %s ===\n', label);
    state_names = {'d(t)', 'phi(t)', 'psi(t)'};
    
    for i = 1:3
        s = y(:, i);   
        final = mean(s(end - floor(length(s)*0.1):end)); % Use average value to calculate the steady-state value, increase robustness
        peak = max(abs(s)); % Calculate the maximum offset value (positive/negative)
        overshoot = (peak - abs(final)) / abs(final) * 100; % Calculate the overshoot

        band = 0.02 * abs(final); % 2% interval
        settled_idx = find(abs(s - final) <= band); % Get the indices of those who fall into +-2% interval
        if isempty(settled_idx)
            T_s = NaN; % Not converge
        else
            d_idx = find(diff(settled_idx) > 1, 1, 'last'); % Find the moment when the system is finally stable
            if isempty(d_idx)
                T_s = t(settled_idx(1));
            else
                T_s = t(settled_idx(d_idx + 1));
            end
        end

        fprintf('Output: %-6s | Overshoot: %7.2f%% | Settling Time: %5.2fs\n', ...
            state_names{i}, overshoot, T_s);
    end
end
