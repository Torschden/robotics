%% Sloshing Dynamics Simulation
% Implementation of the sloshing model: M(q)q_ddot + C(q,q_dot)q_dot + G(q) + D = Q*u
% State vector: q = [x_n, y_n, x_0, y_0]'
% Input vector: u = [u_x, u_y]'

clear; clc; close all;

%% System Parameters
params = struct();
params.omega_n = 2.5;      % Natural frequency [rad/s]
params.zeta_n = 0.05;      % Damping ratio
params.g = 9.81;           % Gravitational acceleration [m/s^2]
params.R = 0.5;            % Tank radius [m]
params.alpha_n = 1.2;      % Nonlinear parameter

% Damping matrix D (assuming diagonal damping)
params.D = diag([0.1, 0.1, 0.05, 0.05]);

%% State and Input Dimensions
% State vector q = [x_n, y_n, x_0, y_0]' (4 states)
% State derivative q_dot = [x_n_dot, y_n_dot, x_0_dot, y_0_dot]' (4 states)
% Total state vector for simulation: [q; q_dot] (8 states)
% Input vector u = [u_x, u_y]' (2 inputs)

n_states = 8;  % Total states [q; q_dot]
n_inputs = 2;  % Number of inputs
n_q = 4;       % Number of generalized coordinates

fprintf('System Dimensions:\n');
fprintf('- Total states: %d (q and q_dot)\n', n_states);
fprintf('- Generalized coordinates: %d\n', n_q);
fprintf('- Inputs: %d\n', n_inputs);

%% Mass Matrix M(q)
function M = mass_matrix(q, params)
    x_n = q(1); y_n = q(2);
    omega_n = params.omega_n; g = params.g;
    
    % Frequency-based coupling terms
    coupling_factor = omega_n^4 / g^2;
    
    M = [1 + coupling_factor * x_n^2,    coupling_factor * x_n * y_n,    0, 0;
         coupling_factor * x_n * y_n,    1 + coupling_factor * y_n^2,    0, 0;
         0,                              0,                              1, 0;
         0,                              0,                              0, 1];
end

%% Coriolis Matrix C(q, q_dot)
function C = coriolis_matrix(q, q_dot, params)
    x_n = q(1); y_n = q(2);
    x_n_dot = q_dot(1); y_n_dot = q_dot(2);
    omega_n = params.omega_n; zeta_n = params.zeta_n; g = params.g;
    
    coupling_factor = omega_n^4 / g^2;
    damping_base = 2 * omega_n * zeta_n;
    
    C = zeros(4, 4);
    
    % First row
    C(1,1) = damping_base * (1 + coupling_factor * x_n^2) + coupling_factor * x_n * x_n_dot;
    C(1,2) = damping_base * coupling_factor * x_n * y_n + coupling_factor * x_n * y_n_dot;
    
    % Second row  
    C(2,1) = damping_base * coupling_factor * x_n * y_n + coupling_factor * y_n * x_n_dot;
    C(2,2) = damping_base * (1 + coupling_factor * y_n^2) + coupling_factor * y_n * y_n_dot;
    
    % Remaining elements are zero (tank motion has no velocity coupling)
end

%% Gravitational/Restoring Force G(q)
function G = gravity_vector(q, params)
    x_n = q(1); y_n = q(2);
    omega_n = params.omega_n; alpha_n = params.alpha_n; R = params.R;
    
    % Nonlinear restoring force
    nonlinear_term = 1 + (alpha_n / R^2) * (x_n^2 + y_n^2);
    
    G = [omega_n^2 * x_n * nonlinear_term;
         omega_n^2 * y_n * nonlinear_term;
         0;
         0];
end

%% Input Matrix Q
function Q = input_matrix()
    Q = [-1,  0;
          0, -1;
          1,  0;
          0,  1];
end

%% Right-Hand Side Function for ODE Integration
function x_dot = sloshing_rhs(t, x, u_func, params)
    % State vector: x = [q; q_dot] where q = [x_n, y_n, x_0, y_0]
    
    % Extract states
    q = x(1:4);
    q_dot = x(5:8);
    
    % Get input at time t
    u = u_func(t);
    
    % Compute system matrices
    M = mass_matrix(q, params);
    C = coriolis_matrix(q, q_dot, params);
    G = gravity_vector(q, params);
    Q = input_matrix();
    
    % Compute acceleration: M(q)q_ddot = -C(q,q_dot)q_dot - G(q) - D*q_dot + Q*u
    q_ddot = M \ (-C * q_dot - G - params.D * q_dot + Q * u);
    
    % State derivative: x_dot = [q_dot; q_ddot]
    x_dot = [q_dot; q_ddot];
end

%% Advanced Trajectory Generation Functions
function profile = generate_trajectory_profile(waypoints, constraints)
    % Generate acceleration profile for point-to-point trajectory
    % waypoints: [time, x_target, y_target] - positions to reach
    % constraints: struct with max_acc, max_vel fields
    
    if nargin < 2
        constraints.max_acc = 0.2;  % m/s²
        constraints.max_vel = 0.5;  % m/s
    end
    
    profile = [];
    current_pos = [0; 0];
    current_vel = [0; 0];
    current_time = 0;
    
    for i = 1:size(waypoints, 1)
        target_time = waypoints(i, 1);
        target_pos = waypoints(i, 2:3)';
        
        % Generate trajectory segment
        segment = plan_trajectory_segment(current_time, target_time, ...
                                        current_pos, target_pos, ...
                                        current_vel, constraints);
        profile = [profile; segment];
        
        % Update current state (simplified - assumes we reach target)
        current_pos = target_pos;
        current_vel = [0; 0];  % Assume we stop at waypoints
        current_time = target_time;
    end
    
    % Add final zero acceleration
    profile = [profile; inf, 0.0, 0.0];
end

function segment = plan_trajectory_segment(t_start, t_end, pos_start, pos_end, vel_start, constraints)
    % Simple trapezoidal velocity profile for now
    % This is a placeholder for more sophisticated trajectory planning
    
    dt = t_end - t_start;
    displacement = pos_end - pos_start;
    
    % Simple constant acceleration approach
    if dt > 0
        required_acc = 2 * displacement / dt^2;
        
        % Limit acceleration
        acc_x = max(-constraints.max_acc, min(constraints.max_acc, required_acc(1)));
        acc_y = max(-constraints.max_acc, min(constraints.max_acc, required_acc(2)));
    else
        acc_x = 0; acc_y = 0;
    end
    
    segment = [t_end, acc_x, acc_y];
end

%% Alternative Input Functions

% Function 1: Step-wise profile (current implementation)
function u = step_profile_input(t)
    profile = [
        2,   0.0,  0.0;    % 0 to 2s: no acceleration
        5,   0.1,  0.0;    % 2 to 5s: +0.1 m/s² in x
        10,  0.0,  0.0;    % 5 to 10s: no acceleration  
        13, -0.1,  0.0;    % 10 to 13s: -0.1 m/s² in x
        inf, 0.0,  0.0     % 13s onwards: no acceleration
    ];
    u = get_acceleration_from_profile(t, profile);
end

% Function 2: Point-to-point trajectory
function u = trajectory_input(t)
    % Define waypoints: [time, x_position, y_position]
    waypoints = [
        10,   0.1,  0.0;    % Move to (0.5, 0) by t=5s
        15,  0.2,  0;    % Move to (0.5, 0.3) by t=10s
        % 15,  0.0,  0.0     % Return to origin by t=15s
    ];
    
    constraints.max_acc = 0.15;  % m/s²
    constraints.max_vel = 0.15;   % m/s
    
    persistent profile_generated profile_data
    if isempty(profile_generated)
        profile_data = generate_trajectory_profile(waypoints, constraints);
        profile_generated = true;
        fprintf('Generated trajectory profile with %d segments\n', size(profile_data, 1));
    end
    
    u = get_acceleration_from_profile(t, profile_data);
end

% Function 3: Custom smooth transitions
function u = smooth_profile_input(t)
    % Smooth transitions between acceleration levels
    if t < 1
        acc_x = 0;
    elseif t < 2
        % Smooth ramp up
        acc_x = 0.1 * (1 - cos(pi * (t - 1))) / 2;
    elseif t < 5
        acc_x = 0.1;
    elseif t < 6
        % Smooth ramp down
        acc_x = 0.1 * (1 + cos(pi * (t - 5))) / 2;
    elseif t < 9
        acc_x = 0;
    elseif t < 10
        % Smooth ramp to negative
        acc_x = -0.1 * (1 - cos(pi * (t - 9))) / 2;
    elseif t < 12
        acc_x = -0.1;
    elseif t < 13
        % Smooth ramp to zero
        acc_x = -0.1 * (1 + cos(pi * (t - 12))) / 2;
    else
        acc_x = 0;
    end
    
    u = [acc_x; 0];
end

%% Main Input Function (switch between different profiles)
function u = input_function(t)
    % Select which input profile to use:
    profile_type = 'smooth';  % Options: 'step', 'trajectory', 'smooth'
    
    switch profile_type
        case 'step'
            u = step_profile_input(t);
        case 'trajectory' 
            u = trajectory_input(t);
        case 't'
            u = trajectory(t);
        case 'smooth'
            u = smooth_profile_input(t);
        otherwise
            u = step_profile_input(t);
    end
end

function acc = get_acceleration_from_profile(t, profile)
    % Extract acceleration based on time and profile
    acc = [0; 0];  % Default acceleration
    
    for i = 1:size(profile, 1)
        if t <= profile(i, 1)
            acc = [profile(i, 2); profile(i, 3)];
            break;
        end
    end
end

%% Simulation Parameters
tspan = [0, 25];  % Time span [s] - extended to cover full profile
dt = 0.01;        % Time step [s]

% Initial conditions: [x_n, y_n, x_0, y_0, x_n_dot, y_n_dot, x_0_dot, y_0_dot]
x0 = [0; 0; 0; 0; 0; 0; 0; 0];  % Small initial displacement

fprintf('\nInitial Conditions:\n');
fprintf('- Fluid displacement (x_n, y_n): [%.3f, %.3f] m\n', x0(1), x0(2));
fprintf('- Tank position (x_0, y_0): [%.3f, %.3f] m\n', x0(3), x0(4));

%% Numerical Integration
fprintf('\nStarting simulation...\n');
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Create anonymous function for ODE solver
rhs_func = @(t, x) sloshing_rhs(t, x, @input_function, params);

% Solve ODE
[t, x] = ode45(rhs_func, tspan, x0, options);

fprintf('Simulation completed. %d time points computed.\n', length(t));

%% Extract Results
q = x(:, 1:4);           % Generalized coordinates
q_dot = x(:, 5:8);       % Generalized velocities

x_n = q(:, 1);           % Fluid x-displacement
y_n = q(:, 2);           % Fluid y-displacement  
x_0 = q(:, 3);           % Tank x-position
y_0 = q(:, 4);           % Tank y-position

%% Compute Input History
u_history = zeros(length(t), 2);
for i = 1:length(t)
    u_history(i, :) = input_function(t(i));
end

%% Trajectory Analysis Functions
function analyze_trajectory(t, x, u_history)
    fprintf('\nTrajectory Analysis:\n');
    
    % Extract positions and velocities
    x_0 = x(:, 3);  % Tank x-position
    y_0 = x(:, 4);  % Tank y-position  
    x_0_dot = x(:, 7);  % Tank x-velocity
    y_0_dot = x(:, 8);  % Tank y-velocity
    
    % Calculate actual trajectory metrics
    final_pos = [x_0(end), y_0(end)];
    max_vel = [max(abs(x_0_dot)), max(abs(y_0_dot))];
    max_acc = [max(abs(u_history(:, 1))), max(abs(u_history(:, 2)))];
    
    % Calculate distance traveled
    distance = 0;
    for i = 2:length(t)
        dx = x_0(i) - x_0(i-1);
        dy = y_0(i) - y_0(i-1);
        distance = distance + sqrt(dx^2 + dy^2);
    end
    
    fprintf('- Final position: [%.3f, %.3f] m\n', final_pos);
    fprintf('- Distance traveled: %.3f m\n', distance);
    fprintf('- Max velocity: [%.3f, %.3f] m/s\n', max_vel);
    fprintf('- Max acceleration: [%.3f, %.3f] m/s²\n', max_acc);
    fprintf('- Final velocity: [%.4f, %.4f] m/s\n', x_0_dot(end), y_0_dot(end));
    
    % Sloshing analysis
    x_n = x(:, 1); y_n = x(:, 2);
    max_slosh = [max(abs(x_n)), max(abs(y_n))];
    final_slosh = [abs(x_n(end)), abs(y_n(end))];
    
    fprintf('- Max sloshing: [%.4f, %.4f] m\n', max_slosh);
    fprintf('- Final sloshing: [%.4f, %.4f] m\n', final_slosh);
    
    % Energy analysis
    KE = 0.5 * sum(x(:, 5:8).^2, 2);
    fprintf('- Final kinetic energy: %.6f J\n', KE(end));
end

%% Analysis and Visualization
analyze_trajectory(t, x, u_history);
plot_trajectory_analysis(t, x, u_history,params);

% Keep the original visualization as well
figure('Position', [100, 100, 1200, 800]); 
function plot_trajectory_analysis(t, x, u_history, params)
    figure('Position', [100, 100, 1400, 900]);
    
    % Extract key variables
    x_n = x(:, 1); y_n = x(:, 2);
    x_0 = x(:, 3); y_0 = x(:, 4);
    x_0_dot = x(:, 7); y_0_dot = x(:, 8);
    
    % Subplot 1: Tank trajectory (X-Y plot)
    subplot(2, 4, 1);
    plot(x_0, y_0, 'b-', 'LineWidth', 2); hold on;
    plot(x_0(1), y_0(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    plot(x_0(end), y_0(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    xlabel('X Position [m]'); ylabel('Y Position [m]');
    title('Tank Trajectory');
    legend('Path', 'Start', 'End', 'Location', 'best');
    grid on; axis equal;
    
    % Subplot 2: Tank position vs time
    subplot(2, 4, 2);
    plot(t, x_0, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, y_0, 'r--', 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Position [m]');
    title('Tank Position vs Time');
    legend('x_0', 'y_0', 'Location', 'best');
    grid on;
    
    % Subplot 3: Tank velocity vs time
    subplot(2, 4, 3);
    plot(t, x_0_dot, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, y_0_dot, 'r--', 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Velocity [m/s]');
    title('Tank Velocity vs Time');
    legend('ẋ_0', 'ẏ_0', 'Location', 'best');
    grid on;
    
    % Subplot 4: Input acceleration
    subplot(2, 4, 4);
    plot(t, u_history(:, 1), 'k-', 'LineWidth', 2); hold on;
    plot(t, u_history(:, 2), 'c--', 'LineWidth', 2);
    xlabel('Time [s]'); ylabel('Acceleration [m/s²]');
    title('Control Input');
    legend('u_x', 'u_y', 'Location', 'best');
    grid on;
    
    % Subplot 5: Fluid sloshing displacements  
    subplot(2, 4, 5);
    plot(t, x_n, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, y_n, 'r--', 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Displacement [m]');
    title('Fluid Sloshing');
    legend('x_n', 'y_n', 'Location', 'best');
    grid on;
    
    % Subplot 6: Sloshing phase portrait
    subplot(2, 4, 6);
    plot(x_n, y_n, 'b-', 'LineWidth', 1); hold on;
    plot(x_n(1), y_n(1), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    plot(x_n(end), y_n(end), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    xlabel('x_n [m]'); ylabel('y_n [m]');
    title('Sloshing Phase Portrait');
    grid on; axis equal;
    
    % Subplot 7: Velocity profile magnitude
    subplot(2, 4, 7);
    vel_mag = sqrt(x_0_dot.^2 + y_0_dot.^2);
    plot(t, vel_mag, 'g-', 'LineWidth', 2);
    xlabel('Time [s]'); ylabel('Speed [m/s]');
    title('Tank Speed');
    grid on;
    
    % Subplot 8: Energy evolution
    subplot(2, 4, 8);
    KE = 0.5 * sum(x(:, 5:8).^2, 2);
    PE_slosh = 0.5 * params.omega_n^2 * (x_n.^2 + y_n.^2);  % Approximate
    Total_E = KE + PE_slosh;
    
    plot(t, KE, 'r-', 'LineWidth', 1); hold on;
    plot(t, PE_slosh, 'b-', 'LineWidth', 1);
    plot(t, Total_E, 'k--', 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('Energy');
    title('System Energy');
    legend('Kinetic', 'Potential', 'Total', 'Location', 'best');
    grid on;
    
    sgtitle('Comprehensive Trajectory Analysis');
end

% Subplot 1: Fluid displacements
subplot(2, 3, 1);
plot(t, x_n, 'b-', 'LineWidth', 1.5); hold on;
plot(t, y_n, 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Displacement [m]');
title('Fluid Sloshing Displacements');
legend('x_n', 'y_n', 'Location', 'best');
grid on;

% Subplot 2: Tank positions
subplot(2, 3, 2);
plot(t, x_0, 'g-', 'LineWidth', 1.5); hold on;
plot(t, y_0, 'm--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Position [m]');
title('Tank Positions');
legend('x_0', 'y_0', 'Location', 'best');
grid on;

% Subplot 3: Phase portrait (fluid motion)
subplot(2, 3, 3);
plot(x_n, y_n, 'b-', 'LineWidth', 1);
xlabel('x_n [m]'); ylabel('y_n [m]');
title('Fluid Motion Phase Portrait');
axis equal; grid on;

% Subplot 4: Input forces
subplot(2, 3, 4);
plot(t, u_history(:, 1), 'k-', 'LineWidth', 1.5); hold on;
plot(t, u_history(:, 2), 'c--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Input [m/s²]');
title('Tank Excitation Forces');
legend('u_x', 'u_y', 'Location', 'best');
grid on;

% Subplot 5: Velocities
subplot(2, 3, 5);
plot(t, q_dot(:, 1), 'b:', 'LineWidth', 1); hold on;
plot(t, q_dot(:, 2), 'r:', 'LineWidth', 1);
xlabel('Time [s]'); ylabel('Velocity [m/s]');
title('Fluid Velocities');
legend('ẋ_n', 'ẏ_n', 'Location', 'best');
grid on;

% Subplot 6: Energy analysis
subplot(2, 3, 6);
% Kinetic energy approximation
KE = 0.5 * sum(q_dot.^2, 2);
% Potential energy approximation  
PE = 0.5 * params.omega_n^2 * (x_n.^2 + y_n.^2);
Total_E = KE + PE;

plot(t, KE, 'r-', 'LineWidth', 1); hold on;
plot(t, PE, 'b-', 'LineWidth', 1);
plot(t, Total_E, 'k--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Energy');
title('System Energy');
legend('Kinetic', 'Potential', 'Total', 'Location', 'best');
grid on;

sgtitle('Sloshing Dynamics Simulation Results');

%% Summary Statistics
fprintf('\nSimulation Results Summary:\n');
fprintf('- Max fluid displacement |x_n|: %.4f m\n', max(abs(x_n)));
fprintf('- Max fluid displacement |y_n|: %.4f m\n', max(abs(y_n)));
fprintf('- Max tank displacement |x_0|: %.4f m\n', max(abs(x_0)));
fprintf('- Max tank displacement |y_0|: %.4f m\n', max(abs(y_0)));
fprintf('- Final total energy: %.6f J\n', Total_E(end));

%% System Analysis Function
function analyze_system(params)
    fprintf('\nSystem Analysis:\n');
    fprintf('- Natural frequency ω_n: %.2f rad/s (%.2f Hz)\n', ...
            params.omega_n, params.omega_n/(2*pi));
    fprintf('- Damping ratio ζ_n: %.3f\n', params.zeta_n);
    fprintf('- Coupling strength (ω_n⁴/g²): %.4f\n', params.omega_n^4/params.g^2);
    
    % Linearized natural frequency
    omega_linear = sqrt(params.omega_n^2);
    fprintf('- Linearized frequency: %.2f rad/s\n', omega_linear);
end

analyze_system(params);