% Corrected MATLAB code for free swinging pendulum simulation
% Fixed: Missing x0 definition and variable scope issues

% Parameters (same as vessel 1)
g = 9.81;                   % gravity [m/s^2]
mass = 0.5;                 % liquid mass [kg]
h_glass = 0.12;             % glass height [m]
l = h_glass/2;              % pendulum length [m]
I = mass * l^2;             % moment of inertia
k = mass * g / l;           % spring constant (gravity restoring force)
c = 0.03;                   % damping coefficient

% Initial conditions
theta0 = deg2rad(10);       % initial deflection [rad]
dtheta0 = 0;                % initial angular velocity [rad/s]
x0 = [theta0; dtheta0];     % FIXED: Define initial state vector

Tmax = 6;                   % simulation time [s]

% ODE function for free pendulum
% FIXED: Pass parameters to avoid scope issues
function dx = free_pendulum(t, x, mass, g, l, I, k, c)
    theta = x(1);
    dtheta = x(2);
    
    % Nonlinear pendulum equation
    ddtheta = (-c*dtheta - k*theta + mass*g*l*sin(theta)) / I;
    
    dx = zeros(2,1);
    dx(1) = dtheta;
    dx(2) = ddtheta;
end

% Solve ODE with parameter passing
% FIXED: Use anonymous function to pass parameters
ode_func = @(t,x) free_pendulum(t, x, mass, g, l, I, k, c);
[t, sol] = ode45(ode_func, [0 Tmax], x0);

% Plot results
figure;
plot(t, rad2deg(sol(:,1)), 'b-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Pendulum Angle [deg]');
title('Free Swinging Pendulum Simulation');
grid on;

% Display system properties
omega_n = sqrt(g/l1);        % Natural frequency
zeta = c / (2 * sqrt(k * I)); % Damping ratio
omega_d = omega_n * sqrt(1 - zeta^2); % Damped frequency (if underdamped)

fprintf('System Properties:\n');
fprintf('Natural frequency: %.2f rad/s (%.2f Hz)\n', omega_n, omega_n/(2*pi));
fprintf('Damping ratio: %.4f\n', zeta);
if zeta < 1
    fprintf('Damped frequency: %.2f rad/s (%.2f Hz)\n', omega_d, omega_d/(2*pi));
    fprintf('System is underdamped\n');
elseif zeta == 1
    fprintf('System is critically damped\n');
else
    fprintf('System is overdamped\n');
end