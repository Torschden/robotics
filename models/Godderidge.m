% Sloshing Tank Simulation Model
% Governing equations:
% θ₁̇ = θ₂
% θ₂̇ = -δ₃θ₁̇³ - δθ₁̇ + (1/l)[-gR(θ₁) + A(t)] - [a(θ/θ_crit)^b + c(θ/θ_crit)^(2d)θ₁̇^m]
% A(t) = -x₀ω²sin(ωt)

clear; clc; close all;

%% Parameters Definition
% System parameters (adjust these based on your specific tank)
params.g = 9.81;        % Gravitational acceleration [m/s²]
params.l = 0.323;         % Characteristic length [m]
params.delta = 0.032;     % Linear damping coefficient
params.delta3 = 0.0032;   % Cubic damping coefficient

% Nonlinear restoring force parameters
params.theta_crit = pi/6;  % Critical angle [rad]
params.a = 1.0;            % Nonlinear coefficient a
params.b = 2.0;            % Nonlinear exponent b
params.c = 0.5;            % Nonlinear coefficient c
params.d = 1.0;            % Nonlinear exponent d (note: 2d in equation)
params.m = 1.0;            % Velocity exponent m

% Forcing parameters
params.x0 = 0;        % Forcing amplitude [m]
params.omega = 2.0;     % Forcing frequency [rad/s]

% Simulation parameters
t_span = [0, 50];       % Time span [s]
dt = 0.01;              % Time step [s]
t = t_span(1):dt:t_span(2);

% Initial conditions [theta1, theta2]
y0 = [pi/2; 0];          % Initial angle and angular velocity

%% Solve the differential equation
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_ode, y] = ode45(@(t, y) sloshing_ode(t, y, params), t, y0, options);

%% Extract results
theta1 = y(:, 1);       % Angular displacement
theta2 = y(:, 2);       % Angular velocity

%% Calculate additional quantities
A_t = -params.x0 * params.omega^2 * sin(params.omega * t_ode);  % Forcing function
theta1_dot = theta2;    % Angular velocity
theta2_dot = zeros(size(t_ode));

% Calculate acceleration for each time point
for i = 1:length(t_ode)
    theta2_dot(i) = calculate_theta2_dot(t_ode(i), [theta1(i); theta2(i)], params);
end

%% Plotting
figure('Position', [100, 100, 1200, 800]);

% Time series plots
subplot(2, 3, 1);
plot(t_ode, theta1 * 180/pi, 'b-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Angular Displacement [deg]');
title('Angular Displacement vs Time');
grid on;

subplot(2, 3, 2);
plot(t_ode, theta2, 'r-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Angular Velocity [rad/s]');
title('Angular Velocity vs Time');
grid on;

subplot(2, 3, 3);
plot(t_ode, theta2_dot, 'g-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Angular Acceleration [rad/s²]');
title('Angular Acceleration vs Time');
grid on;

% Phase portrait
subplot(2, 3, 4);
plot(theta1 * 180/pi, theta2, 'b-', 'LineWidth', 1.5);
xlabel('Angular Displacement [deg]');
ylabel('Angular Velocity [rad/s]');
title('Phase Portrait');
grid on;

% Forcing function
subplot(2, 3, 5);
plot(t_ode, A_t, 'k-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Forcing A(t) [m/s²]');
title('External Forcing');
grid on;

% Energy analysis (approximate)
KE = 0.5 * params.l^2 * theta2.^2;  % Kinetic energy (simplified)
PE = params.g * params.l * (1 - cos(theta1));  % Potential energy (simplified)
Total_E = KE + PE;

subplot(2, 3, 6);
plot(t_ode, KE, 'r-', 'LineWidth', 1.5); hold on;
plot(t_ode, PE, 'b-', 'LineWidth', 1.5);
plot(t_ode, Total_E, 'k--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Energy');
title('Energy Components');
legend('Kinetic', 'Potential', 'Total', 'Location', 'best');
grid on;

sgtitle('Sloshing Tank Simulation Results');

%% Display some statistics
fprintf('Simulation Results Summary:\n');
fprintf('Max angular displacement: %.3f deg\n', max(abs(theta1)) * 180/pi);
fprintf('Max angular velocity: %.3f rad/s\n', max(abs(theta2)));
fprintf('RMS angular displacement: %.3f deg\n', rms(theta1) * 180/pi);
fprintf('RMS angular velocity: %.3f rad/s\n', rms(theta2));

%% Function definitions
function dydt = sloshing_ode(t, y, params)
    % Extract state variables
    theta1 = y(1);
    theta2 = y(2);  % theta1_dot
    
    % Calculate theta2_dot
    theta2_dot = calculate_theta2_dot(t, y, params);
    
    % State derivatives
    dydt = [theta2; theta2_dot];
end

function theta2_dot = calculate_theta2_dot(t, y, params)
    % Extract parameters
    g = params.g;
    l = params.l;
    delta = params.delta;
    delta3 = params.delta3;
    theta_crit = params.theta_crit;
    a = params.a;
    b = params.b;
    c = params.c;
    d = params.d;
    m = params.m;
    x0 = params.x0;
    omega = params.omega;
    
    % Extract state variables
    theta1 = y(1);
    theta1_dot = y(2);
    
    % External forcing
    A_t = -x0 * omega^2 * sin(omega * t);
    
    % Restoring force R(theta1) - assuming linear for now
    % You may need to modify this based on your specific R(theta1)
    R_theta1 = sin(theta1);  % Common assumption for small angles: R(θ) ≈ sin(θ)
    
    % Nonlinear damping term
    theta_ratio = theta1 / theta_crit;
    nonlinear_damping = a * (abs(theta_ratio))^b + c * (abs(theta_ratio))^(2*d) * (abs(theta1_dot))^m;
    
    % Apply sign to nonlinear damping
    if theta1_dot ~= 0
        nonlinear_damping = nonlinear_damping * sign(theta1_dot);
    end
    
    % Calculate theta2_dot according to the governing equation
    theta2_dot = -delta3 * theta1_dot^3 - delta * theta1_dot + ...
                 (1/l) * (-g * R_theta1 + A_t) - nonlinear_damping;
end

%% Parameter sensitivity analysis function (optional)
function parameter_sensitivity_analysis()
    % This function can be called separately to analyze parameter effects
    fprintf('Parameter Sensitivity Analysis:\n');
    fprintf('Modify parameters in the main script and observe changes in:\n');
    fprintf('1. Maximum displacement amplitude\n');
    fprintf('2. Resonance frequency\n');
    fprintf('3. Damping characteristics\n');
    fprintf('4. Nonlinear behavior onset\n');
end