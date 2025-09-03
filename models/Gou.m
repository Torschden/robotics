%% Liquid Sloshing Model in Cylindrical Vessel
% Based on the paper: "Synchronous Hopf-Bifurcation and Damping Osmosis
% Phenomena of Liquid-Spacecraft Coupled System" by Gou et al. (2001)

clear all; close all; clc;

%% Physical Parameters
% Tank geometry
d = 3.1;        % Tank diameter [dimensionless]
h = 3.0;        % Liquid depth [dimensionless]
a = d/2;        % Tank radius

% Liquid properties
rho = 1000;     % Liquid density [kg/m³] (water)
g = 9.81;       % Gravitational acceleration [m/s²]
sigma = 0.073;  % Surface tension coefficient [N/m] (water)

% System parameters
M = 1;          % Dry mass [dimensionless]
mu = 0.16;      % Mass ratio (liquid mass / dry mass)
k = 1;          % Spring stiffness [dimensionless]

% Damping parameters
zeta_x = 0.05;  % Structural damping ratio
zeta_1 = 0.0348; % Primary mode damping ratio
zeta_2 = 0.0348; % Secondary mode damping ratio
zeta_3 = 0.035;  % Axisymmetric mode damping ratio

% Modal parameters (from Table 1 in paper)
% Primary modes (order ε¹)
k_11 = 1.841/a; % First antisymmetric mode wave number
lambda_11 = k_11; % Assuming small contact angle hysteresis

% Secondary modes (order ε²)
k_01 = 3.832/a; % Axisymmetric mode wave number
k_21 = 3.054/a; % Second antisymmetric mode wave number

% Natural frequencies
nu_1 = sqrt(g * k_11 * tanh(k_11 * h)) / sqrt(k/M); % Primary mode frequency ratio
nu_2 = nu_1; % Assuming symmetric system
nu_3 = sqrt(g * k_01 * tanh(k_01 * h)) / sqrt(k/M); % Axisymmetric mode frequency ratio

%% Excitation Parameters
E_ex = 0.01;    % Excitation amplitude [dimensionless]
nu_x = 0.90;    % Excitation frequency ratio
omega_x = nu_x * sqrt(k/M); % Excitation frequency [rad/s]

%% Modal Coefficients (simplified - key coupling terms)
% These would normally be computed from complex modal integrals
% Using representative values based on the paper's numerical examples

% Coupling coefficients
alpha_1 = 0.5;  % X-q1 coupling
alpha_2 = 0.5;  % Y-q2 coupling

% Nonlinear coefficients (representative values)
beta_111 = 0.1;
beta_1111 = 0.05;
beta_1122 = 0.03;
gamma_113 = 0.02;
gamma_1111 = 0.01;

%% Time Integration Setup
t_end = 200;    % Simulation time
dt = 0.01;      % Time step
t = 0:dt:t_end;
N = length(t);

% State vector: [X, X_dot, q1, q1_dot, q2, q2_dot, q3, q3_dot]
% X: structural displacement
% q1, q2: primary sloshing modes (planar and nonplanar)
% q3: axisymmetric mode

n_states = 8;
y = zeros(N, n_states);

% Initial conditions (small perturbations as in the paper)
y(1, :) = [0, 0, 0, 0, 0, 0, 0, 0]; % Small nonplanar mode perturbation

%% Main Integration Loop
fprintf('Running sloshing simulation...\n');
tic;

for i = 1:N-1
    % Current state
    X = y(i, 1);
    X_dot = y(i, 2);
    q1 = y(i, 3);
    q1_dot = y(i, 4);
    q2 = y(i, 5);
    q2_dot = y(i, 6);
    q3 = y(i, 7);
    q3_dot = y(i, 8);
    
    % Excitation force
    F_ex = E_ex * cos(nu_x * sqrt(k/M) * t(i));
    
    % Mass matrix (time-varying due to liquid motion)
    M11 = 1 + mu;
    M33 = mu * (1 + 0.1 * q1^2); % Simplified nonlinear mass coupling
    M55 = mu * (1 + 0.1 * q2^2);
    M77 = mu * (1 + 0.05 * (q1^2 + q2^2));
    
    % Equations of motion (simplified from Eq. 22 in paper)
    
    % Structural equation
    X_ddot = (F_ex - 2*zeta_x*X_dot - X - alpha_1*q1_dot) / M11;
    
    % Primary planar mode (q1)
    q1_ddot = (-alpha_1*X_ddot - 2*zeta_1*nu_1*q1_dot - nu_1^2*q1 ...
               - beta_111*q1*q1_dot^2 - beta_1111*q1^3 ...
               - beta_1122*q1*q2^2 - gamma_113*q1*q3) / M33;
    
    % Primary nonplanar mode (q2) 
    q2_ddot = (-2*zeta_2*nu_2*q2_dot - nu_2^2*q2 ...
               - beta_111*q2*q2_dot^2 - beta_1111*q2^3 ...
               - beta_1122*q2*q1^2 - gamma_113*q2*q3) / M55;
    
    % Axisymmetric mode (q3)
    q3_ddot = (-2*zeta_3*nu_3*q3_dot - nu_3^2*q3 ...
               - 0.5*gamma_113*(q1^2 + q2^2)) / M77;
    
    % Update state vector using 4th-order Runge-Kutta
    k1 = dt * dynamics_rhs([X; X_dot; q1; q1_dot; q2; q2_dot; q3; q3_dot], ...
                          F_ex, alpha_1, alpha_2, beta_111, beta_1111, beta_1122, ...
                          gamma_113, zeta_x, zeta_1, zeta_2, zeta_3, ...
                          nu_1, nu_2, nu_3, mu);
    
    k2 = dt * dynamics_rhs(y(i,:)' + k1/2, F_ex, alpha_1, alpha_2, ...
                          beta_111, beta_1111, beta_1122, gamma_113, ...
                          zeta_x, zeta_1, zeta_2, zeta_3, nu_1, nu_2, nu_3, mu);
    
    k3 = dt * dynamics_rhs(y(i,:)' + k2/2, F_ex, alpha_1, alpha_2, ...
                          beta_111, beta_1111, beta_1122, gamma_113, ...
                          zeta_x, zeta_1, zeta_2, zeta_3, nu_1, nu_2, nu_3, mu);
    
    k4 = dt * dynamics_rhs(y(i,:)' + k3, F_ex, alpha_1, alpha_2, ...
                          beta_111, beta_1111, beta_1122, gamma_113, ...
                          zeta_x, zeta_1, zeta_2, zeta_3, nu_1, nu_2, nu_3, mu);
    
    y(i+1, :) = y(i, :) + (k1 + 2*k2 + 2*k3 + k4)'/6;
end

elapsed_time = toc;
fprintf('Simulation completed in %.2f seconds.\n', elapsed_time);

%% Post-Processing and Visualization

% Extract time histories
X_history = y(:, 1);
q1_history = y(:, 3);
q2_history = y(:, 5);
q3_history = y(:, 7);

% Calculate wave height at a point (r=a, theta=0)
% Based on Eq. (7) and Table 1 from the paper
r_point = a;
theta_point = 0;

% Bessel functions for modal shapes
J1_k11a = besselj(1, k_11*r_point);
J0_k01a = besselj(0, k_01*r_point);
J2_k21a = besselj(2, k_21*r_point);

% Normalization coefficients (simplified)
B1 = 1; B2 = 1; B3 = 1;

% Wave height calculation
eta_history = B1 * J1_k11a * cos(theta_point) * q1_history + ...
              B2 * J1_k11a * sin(theta_point) * q2_history + ...
              B3 * J0_k01a * q3_history;

%% Plotting Results
figure('Position', [100, 100, 1200, 800]);

% Time histories
subplot(3, 2, 1);
plot(t, X_history, 'b-', 'LineWidth', 1.5);
xlabel('Time [dimensionless]');
ylabel('X [dimensionless]');
title('Structural Displacement');
grid on;

subplot(3, 2, 2);
plot(t, q1_history, 'r-', 'LineWidth', 1.5);
xlabel('Time [dimensionless]');
ylabel('q_1 [dimensionless]');
title('Primary Planar Mode');
grid on;

subplot(3, 2, 3);
plot(t, q2_history, 'g-', 'LineWidth', 1.5);
xlabel('Time [dimensionless]');
ylabel('q_2 [dimensionless]');
title('Primary Nonplanar Mode');
grid on;

subplot(3, 2, 4);
plot(t, q3_history, 'm-', 'LineWidth', 1.5);
xlabel('Time [dimensionless]');
ylabel('q_3 [dimensionless]');
title('Axisymmetric Mode');
grid on;

subplot(3, 2, 5);
plot(t, eta_history, 'k-', 'LineWidth', 1.5);
xlabel('Time [dimensionless]');
ylabel('\eta [dimensionless]');
title('Wave Height at (r=a, \theta=0)');
grid on;

% Phase portrait of nonplanar mode
subplot(3, 2, 6);
plot(q2_history, y(:, 6), 'g-', 'LineWidth', 1.5);
xlabel('q_2 [dimensionless]');
ylabel('dq_2/dt [dimensionless]');
title('Nonplanar Mode Phase Portrait');
grid on;
axis equal;

sgtitle('Liquid Sloshing Dynamics in Cylindrical Vessel', 'FontSize', 14, 'FontWeight', 'bold');

%% Frequency Analysis
figure('Position', [1350, 100, 600, 800]);

% FFT of primary modes
Fs = 1/dt;
L = length(t);
f = Fs*(0:(L/2))/L;

% Skip initial transient
skip_samples = round(0.2 * N);
q1_steady = q1_history(skip_samples:end);
q2_steady = q2_history(skip_samples:end);

Y1 = fft(q1_steady);
Y2 = fft(q2_steady);
P1 = abs(Y1/length(q1_steady));
P2 = abs(Y2/length(q2_steady));
P1 = P1(1:floor(length(q1_steady)/2)+1);
P2 = P2(1:floor(length(q2_steady)/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
P2(2:end-1) = 2*P2(2:end-1);

f_plot = f(1:length(P1));

subplot(2, 1, 1);
semilogy(f_plot, P1, 'r-', 'LineWidth', 1.5);
xlabel('Frequency [Hz]');
ylabel('|FFT(q_1)|');
title('Frequency Spectrum - Primary Planar Mode');
grid on;
xlim([0, 2]);

subplot(2, 1, 2);
semilogy(f_plot, P2, 'g-', 'LineWidth', 1.5);
xlabel('Frequency [Hz]');
ylabel('|FFT(q_2)|');
title('Frequency Spectrum - Primary Nonplanar Mode');
grid on;
xlim([0, 2]);

%% Display System Information
fprintf('\n=== System Parameters ===\n');
fprintf('Tank diameter: %.2f\n', d);
fprintf('Liquid depth: %.2f\n', h);
fprintf('Mass ratio (μ): %.4f\n', mu);
fprintf('Primary mode frequency ratio (ν₁): %.4f\n', nu_1);
fprintf('Excitation frequency ratio (νₓ): %.4f\n', nu_x);
fprintf('Excitation amplitude: %.4f\n', E_ex);

%% Function Definitions

function dydt = dynamics_rhs(y, F_ex, alpha_1, alpha_2, beta_111, beta_1111, ...
                            beta_1122, gamma_113, zeta_x, zeta_1, zeta_2, ...
                            zeta_3, nu_1, nu_2, nu_3, mu)
    % Right-hand side of the dynamics equations
    % State vector: [X, X_dot, q1, q1_dot, q2, q2_dot, q3, q3_dot]
    
    X = y(1); X_dot = y(2);
    q1 = y(3); q1_dot = y(4);
    q2 = y(5); q2_dot = y(6);
    q3 = y(7); q3_dot = y(8);
    
    % Mass matrix terms (simplified)
    M11 = 1 + mu;
    M33 = mu * (1 + 0.1 * q1^2);
    M55 = mu * (1 + 0.1 * q2^2);
    M77 = mu * (1 + 0.05 * (q1^2 + q2^2));
    
    % Equations of motion
    X_ddot = (F_ex - 2*zeta_x*X_dot - X - alpha_1*q1_dot) / M11;
    
    q1_ddot = (-alpha_1*X_ddot - 2*zeta_1*nu_1*q1_dot - nu_1^2*q1 ...
               - beta_111*q1*q1_dot^2 - beta_1111*q1^3 ...
               - beta_1122*q1*q2^2 - gamma_113*q1*q3) / M33;
    
    q2_ddot = (-2*zeta_2*nu_2*q2_dot - nu_2^2*q2 ...
               - beta_111*q2*q2_dot^2 - beta_1111*q2^3 ...
               - beta_1122*q2*q1^2 - gamma_113*q2*q3) / M55;
    
    q3_ddot = (-2*zeta_3*nu_3*q3_dot - nu_3^2*q3 ...
               - 0.5*gamma_113*(q1^2 + q2^2)) / M77;
    
    dydt = [X_dot; X_ddot; q1_dot; q1_ddot; q2_dot; q2_ddot; q3_dot; q3_ddot];
end