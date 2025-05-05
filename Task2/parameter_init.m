function output = parameter_init()
    % Metriculation Number: A0313771H
    a = 3; b = 7; c = 7; d = 1;

    g = 9.8; % Gravitational acceleration
    
    M_f = 2.14 + c/20; % Mass of the front wheel
    M_r = 5.91 - b/10; % Mass of the rear wheel
    M_c = 1.74; % Mass of the cart system
    
    L_F = 0.133; % Horizontal length from a front wheel rotation axis to the steering axis
    L_Ff = 0.05; % Horizontal length from a front wheel rotation axis to a center-of-gravity of part of front wheel
    L_R = 0.308 + (a - d)/100; % Horizontal length from a rear wheel rotation axis to the steering axis
    L_r = 0.128; % Horizontal length from a rear wheel rotation axis to a center-of-gravity of part of rear wheel
    L_c = 0.259; % Horizontal length from a rear wheel rotation axis to a center-of-gravity of the cart system
    
    J_x = 0.5 + (c - d)/100; % Moment of inertia around center-of-gravity x axially
    
    H_f = 0.18; % Vertical length from a floor to a center-of-gravity of the front wheel
    H_r = 0.161; % Vertical length from a floor to a center-of-gravity of the rear wheel
    H_c = 0.098; % Vertical length from a floor to a center-of-gravity of the front wheel
    
    mu_x = 3.33 - b/20 + a*c/60; % Viscous coefficient around x axis
    
    alpha = 15.5 - a/3 + b/2;
    gamma = 11.5+ (a - c)/(b + d + 3);
    beta = 27.5 - d/2;
    delta = 60 + (a - d)*c/10;
    
    %% Parameters in the matrices
    den = M_f*H_f^2 + M_r*H_r^2 + M_c*H_c^2 + J_x;
    a_51 = -M_c*g/den;
    a_52 = (M_f*H_f + M_r*H_r + M_c*H_c)*g/den;
    a_53 = (M_r*L_r*L_F + M_c*L_c*L_F + M_f*L_Ff*L_R)*g/((L_R+L_F)*den);
    a_54 = -M_c*H_c*alpha/den;
    a_55 = -mu_x/den;
    a_56 = M_f*H_f*L_Ff*gamma/den;
    b_51 = M_c*H_c*beta/den;
    b_52 = -M_f*H_f*L_Ff*delta/den;
    
    A = [0, 0, 0, 1, 0, 0;
         0, 0, 0, 0, 1, 0;
         0, 0, 0, 0, 0, 1;
         0, 6.5, -10, -alpha, 0, 0;
         a_51, a_52, a_53, a_54, a_55, a_56;
         5, -3.6, 0, 0, 0, -gamma];

    B = [0, 0;
         0, 0;
         0, 0;
         beta, 11.2;
         b_51, b_52;
         40, delta];

    C = [1, 0, 0, 0, 0, 0;
         0, 1, 0, 0, 0, 0;
         0, 0, 1, 0, 0, 0];

    y_sp = -0.1*C/A*B*[-0.5 + (a - b)/20; 0.1 + (b - c)/(a + d + 10)];

    x_0 = [0.2; -0.1; 0.15; -1; 0.8; 0];
    
    output = {A, B, C, y_sp, x_0};
end