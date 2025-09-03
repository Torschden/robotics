%% Sloshing Suppression Simulation
% Based on "A Dynamic Optimization Approach for Sloshing Free Transport 
% of Liquid Filled Containers using an Industrial Robot"
% This simulation demonstrates the difference between classical trajectory
% planning and optimized sloshing-free trajectory planning

clear; clc; close all;

%% Parameters from the paper (Table I)
% Liquid parameters
m = 0.5;           % Mass of liquid [kg]
l = 0.021;         % Pendulum length [m]
d = 1.51;          % Damping constant [N*s/m]
g = 9.81;          % Gravitational acceleration [m/s^2]

% Derived parameters
omega = sqrt(g/l); % Natural frequency [rad/s]
zeta = d/(2*m*omega); % Damping ratio

% Simulation parameters
dt = 0.001;        % Time step [s]
t_total = 3;       % Total simulation time [s]
t = 0:dt:t_total;  % Time vector
N = length(t);

% Motion parameters
y_start = -0.27;   % Initial position [m]
y_end = 0.27;      % Final position [m]
distance = y_end - y_start; % Total distance [m]

%% Classical Trajectory (Constant Acceleration/Deceleration)
fprintf('Simulating Classical Trajectory...\n');

% Classical trajectory parameters
t_accel = 0.5;     % Acceleration time [s]
t_decel = 0.5;     % Deceleration time [s]
t_const = 0.5;     % Constant velocity time [s]
t_start = 0.5;     % Start time [s]

% Calculate required acceleration
v_max = distance / (t_accel/2 + t_const + t_decel/2);
a_max = v_max / t_accel;

% Initialize arrays
y_classic = zeros(size(t));
v_classic = zeros(size(t));
a_classic = zeros(size(t));
phi_classic = zeros(size(t));
phi_dot_classic = zeros(size(t));

% Generate classical trajectory
for i = 1:N
    if t(i) < t_start
        % Before motion starts
        y_classic(i) = y_start;
        v_classic(i) = 0;
        a_classic(i) = 0;
    elseif t(i) < t_start + t_accel
        % Acceleration phase
        t_rel = t(i) - t_start;
        y_classic(i) = y_start + 0.5 * a_max * t_rel^2;
        v_classic(i) = a_max * t_rel;
        a_classic(i) = a_max;
    elseif t(i) < t_start + t_accel + t_const
        % Constant velocity phase
        t_rel = t(i) - t_start - t_accel;
        y_classic(i) = y_start + 0.5 * a_max * t_accel^2 + v_max * t_rel;
        v_classic(i) = v_max;
        a_classic(i) = 0;
    elseif t(i) < t_start + t_accel + t_const + t_decel
        % Deceleration phase
        t_rel = t(i) - t_start - t_accel - t_const;
        y_classic(i) = y_start + 0.5 * a_max * t_accel^2 + v_max * t_const + ...
                      v_max * t_rel - 0.5 * a_max * t_rel^2;
        v_classic(i) = v_max - a_max * t_rel;
        a_classic(i) = -a_max;
    else
        % After motion ends
        y_classic(i) = y_end;
        v_classic(i) = 0;
        a_classic(i) = 0;
    end
end

% Simulate sloshing dynamics for classical trajectory
for i = 2:N
    % Nonlinear pendulum dynamics (equation 2 from paper)
    phi_ddot = (a_classic(i)/l) * cos(phi_classic(i-1)) - ...
               (g/l) * sin(phi_classic(i-1)) - ...
               (d/m) * phi_dot_classic(i-1) * cos(phi_classic(i-1))^2;
    
    % Euler integration
    phi_dot_classic(i) = phi_dot_classic(i-1) + phi_ddot * dt;
    phi_classic(i) = phi_classic(i-1) + phi_dot_classic(i-1) * dt;
end

%% Optimized Trajectory (Sloshing Suppression)
fprintf('Simulating Optimized Trajectory...\n');

% For the optimized trajectory, we'll create a trajectory that compensates
% for sloshing by tilting the container appropriately
y_opt = zeros(size(t));
v_opt = zeros(size(t));
a_opt = zeros(size(t));
phi_container = zeros(size(t)); % Container tilt angle
phi_liquid_opt = zeros(size(t)); % Liquid sloshing angle (should be minimal)
phi_dot_liquid_opt = zeros(size(t));

% Generate smooth trajectory using polynomial
t_motion_start = 0.5;
t_motion_end = 2.5;
motion_duration = t_motion_end - t_motion_start;

% 5th order polynomial coefficients for smooth trajectory
% Boundary conditions: y(0) = y_start, y(T) = y_end, v(0) = v(T) = 0, a(0) = a(T) = 0
A = [0, 0, 0, 0, 0, 1;
     motion_duration^5, motion_duration^4, motion_duration^3, motion_duration^2, motion_duration, 1;
     0, 0, 0, 0, 1, 0;
     5*motion_duration^4, 4*motion_duration^3, 3*motion_duration^2, 2*motion_duration, 1, 0;
     0, 0, 0, 2, 0, 0;
     20*motion_duration^3, 12*motion_duration^2, 6*motion_duration, 2, 0, 0];

b = [y_start; y_end; 0; 0; 0; 0];
coeffs = A \ b;

% Generate optimized trajectory
for i = 1:N
    if t(i) < t_motion_start
        y_opt(i) = y_start;
        v_opt(i) = 0;
        a_opt(i) = 0;
    elseif t(i) <= t_motion_end
        t_rel = t(i) - t_motion_start;
        y_opt(i) = coeffs(1)*t_rel^5 + coeffs(2)*t_rel^4 + coeffs(3)*t_rel^3 + ...
                   coeffs(4)*t_rel^2 + coeffs(5)*t_rel + coeffs(6);
        v_opt(i) = 5*coeffs(1)*t_rel^4 + 4*coeffs(2)*t_rel^3 + 3*coeffs(3)*t_rel^2 + ...
                   2*coeffs(4)*t_rel + coeffs(5);
        a_opt(i) = 20*coeffs(1)*t_rel^3 + 12*coeffs(2)*t_rel^2 + 6*coeffs(3)*t_rel + ...
                   2*coeffs(4);
    else
        y_opt(i) = y_end;
        v_opt(i) = 0;
        a_opt(i) = 0;
    end
end

% Calculate optimal container tilt to suppress sloshing
% From the paper: tilt the container so that liquid surface remains level
for i = 1:N
    if abs(a_opt(i)) > 0.01  % Avoid division by small numbers
        % Calculate required tilt angle to keep liquid surface horizontal
        phi_container(i) = -atan(a_opt(i)/g);
    else
        phi_container(i) = 0;
    end
end

% Simulate liquid dynamics with optimized trajectory
% The key insight: with proper container tilting, relative acceleration is minimized
for i = 2:N
    % Effective acceleration felt by liquid (reduced due to container tilt)
    a_effective = a_opt(i) * cos(phi_container(i)) + g * sin(phi_container(i));
    
    % Pendulum dynamics with effective acceleration
    phi_ddot = (a_effective/l) * cos(phi_liquid_opt(i-1)) - ...
               (g/l) * sin(phi_liquid_opt(i-1)) - ...
               (d/m) * phi_dot_liquid_opt(i-1) * cos(phi_liquid_opt(i-1))^2;
    
    % Euler integration
    phi_dot_liquid_opt(i) = phi_dot_liquid_opt(i-1) + phi_ddot * dt;
    phi_liquid_opt(i) = phi_liquid_opt(i-1) + phi_dot_liquid_opt(i-1) * dt;
end

%% Visualization
fprintf('Creating visualizations...\n');

%% Figure 6 Recreation - Four-plot layout matching the paper
figure('Position', [100, 100, 1000, 800]);

% Plot 1: y position (liquid position)
subplot(4,1,1);
plot(t, y_opt, 'b-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('y^l_0 [m]');
title('Optimal State Trajectory (Figure 6 Recreation)');
grid on; xlim([0, 3]);
ylim([y_start-0.05, y_end+0.05]);

% Plot 2: y velocity (liquid velocity) 
subplot(4,1,2);
plot(t, v_opt, 'b-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('ẏ^l_0 [m s^{-1}]');
grid on; xlim([0, 3]);

% Plot 3: phi angle (liquid surface tilt angle)
subplot(4,1,3);
plot(t, rad2deg(phi_liquid_opt), 'b-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('φ^l_0 [°]');
grid on; xlim([0, 3]);
ylim([-5, 5]);

% Plot 4: phi_dot (angular velocity of liquid surface)
subplot(4,1,4);
plot(t, rad2deg(phi_dot_liquid_opt), 'b-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('φ̇^l_0 [° s^{-1}]');
grid on; xlim([0, 3]);

% Add vertical lines to show motion phases
for subplot_idx = 1:4
    subplot(4,1,subplot_idx);
    hold on;
    xline(t_motion_start, 'k--', 'LineWidth', 1, 'Alpha', 0.5);
    xline(t_motion_end, 'k--', 'LineWidth', 1, 'Alpha', 0.5);
end

% Figure: Sloshing Comparison
figure('Position', [200, 200, 1200, 600]);

subplot(2,1,1);
plot(t, rad2deg(phi_classic), 'r-', 'LineWidth', 2); hold on;
plot(t, rad2deg(phi_liquid_opt), 'b-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Sloshing Angle [°]');
title('Liquid Sloshing Angle Comparison');
legend('Classical Trajectory', 'Optimized Trajectory', 'Location', 'best');
grid on; ylim([-25, 25]);

subplot(2,1,2);
plot(t, rad2deg(phi_container), 'g--', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Container Tilt [°]');
title('Container Tilt Angle (Optimized Trajectory)');
legend('Container Tilt', 'Location', 'best');
grid on;

% Figure 3: Animation Setup
figure('Position', [300, 300, 800, 600]);
axis equal; axis([-0.4, 0.4, -0.4, 0.4]);
xlabel('X Position [m]'); ylabel('Y Position [m]');
title('Container Motion Animation');
hold on; grid on;

% Animation parameters
container_width = 0.1;
container_height = 0.15;
liquid_height = 0.1;

% Run animation
fprintf('Running animation (close figure to continue)...\n');
for i = 1:10:N  % Animate every 10th frame for speed
    clf;
    axis equal; axis([-0.4, 0.4, -0.4, 0.4]);
    xlabel('X Position [m]'); ylabel('Y Position [m]');
    title(sprintf('Container Motion Animation - t = %.2f s', t(i)));
    hold on; grid on;
    
    % Classical trajectory (left side)
    container_x_classic = -0.2;
    container_y_classic = y_classic(i);
    
    % Draw classical container (red)
    rectangle('Position', [container_x_classic - container_width/2, ...
                          container_y_classic - container_height/2, ...
                          container_width, container_height], ...
             'FaceColor', [1, 0.8, 0.8], 'EdgeColor', 'r', 'LineWidth', 2);
    
    % Draw liquid in classical container (tilted due to sloshing)
    liquid_tilt = phi_classic(i);
    liquid_x = [container_x_classic - container_width/2 + 0.01, ...
                container_x_classic + container_width/2 - 0.01];
    liquid_y = [container_y_classic - container_height/2 + 0.02, ...
                container_y_classic - container_height/2 + 0.02] + ...
               [liquid_height * sin(liquid_tilt), -liquid_height * sin(liquid_tilt)];
    fill([liquid_x, fliplr(liquid_x)], ...
         [liquid_y, liquid_y - liquid_height], ...
         'b', 'FaceAlpha', 0.6);
    
    % Optimized trajectory (right side)
    container_x_opt = 0.2;
    container_y_opt = y_opt(i);
    
    % Container rotation matrix
    R = [cos(phi_container(i)), -sin(phi_container(i));
         sin(phi_container(i)), cos(phi_container(i))];
    
    % Draw optimized container (blue, tilted)
    corners = [-container_width/2, -container_height/2;
               container_width/2, -container_height/2;
               container_width/2, container_height/2;
               -container_width/2, container_height/2];
    rotated_corners = (R * corners')';
    rotated_corners(:,1) = rotated_corners(:,1) + container_x_opt;
    rotated_corners(:,2) = rotated_corners(:,2) + container_y_opt;
    
    fill(rotated_corners(:,1), rotated_corners(:,2), ...
         [0.8, 0.8, 1], 'EdgeColor', 'b', 'LineWidth', 2);
    
    % Draw liquid in optimized container (should remain level)
    liquid_corners = [-container_width/2 + 0.01, -container_height/2 + 0.02;
                      container_width/2 - 0.01, -container_height/2 + 0.02;
                      container_width/2 - 0.01, -container_height/2 + 0.02 + liquid_height;
                      -container_width/2 + 0.01, -container_height/2 + 0.02 + liquid_height];
    
    % Apply small sloshing perturbation
    liquid_corners(3:4, 2) = liquid_corners(3:4, 2) + ...
        [liquid_height * sin(phi_liquid_opt(i)), -liquid_height * sin(phi_liquid_opt(i))];
    
    rotated_liquid = (R * liquid_corners')';
    rotated_liquid(:,1) = rotated_liquid(:,1) + container_x_opt;
    rotated_liquid(:,2) = rotated_liquid(:,2) + container_y_opt;
    
    fill(rotated_liquid(:,1), rotated_liquid(:,2), 'c', 'FaceAlpha', 0.8);
    
    % Add labels
    text(-0.2, -0.35, 'Classical', 'HorizontalAlignment', 'center', 'FontSize', 12);
    text(0.2, -0.35, 'Optimized', 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    drawnow;
    pause(0.05);
end

%% Results Summary
fprintf('\n=== SIMULATION RESULTS ===\n');
fprintf('Classical Trajectory:\n');
fprintf('  Max sloshing angle: %.2f°\n', max(abs(rad2deg(phi_classic))));
fprintf('  RMS sloshing angle: %.2f°\n', rms(rad2deg(phi_classic)));

fprintf('\nOptimized Trajectory:\n');
fprintf('  Max sloshing angle: %.2f°\n', max(abs(rad2deg(phi_liquid_opt))));
fprintf('  RMS sloshing angle: %.2f°\n', rms(rad2deg(phi_liquid_opt)));
fprintf('  Max container tilt: %.2f°\n', max(abs(rad2deg(phi_container))));

fprintf('\nSloshing Reduction: %.1f%%\n', ...
    100 * (1 - rms(rad2deg(phi_liquid_opt)) / rms(rad2deg(phi_classic))));

fprintf('\nSimulation completed successfully!\n');