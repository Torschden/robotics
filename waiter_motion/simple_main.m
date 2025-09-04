%% simple_main.m - Simplified main script that actually works
%
% This is a simplified version that focuses on the core functionality
% without complex optimization that might fail due to toolbox dependencies

clear all;
close all;
clc;

fprintf('=== Simplified Sloshing Controller Demo ===\n\n');

%% Create controller instance
controller = SloshingController();

%% Define a simple trajectory scenario
fprintf('Defining trajectory scenario...\n');

% Simple start and end positions
q_start = [0, 0, 0, 0, 0, 0]';
q_end = [pi/4, pi/6, -pi/6, 0, 0, pi/4]';

fprintf('Start configuration: [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f] rad\n', q_start);
fprintf('End configuration:   [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f] rad\n\n', q_end);

%% Generate simple trajectory without optimization
fprintf('Generating simple trajectory (no optimization)...\n');

% Time parameters
t_total = 3.0; % Total trajectory time
N_points = 50; % Number of trajectory points
time_vec = linspace(0, t_total, N_points);

% Generate simple polynomial trajectory
q_traj = generateSimpleTrajectory(q_start, q_end, time_vec);

fprintf('Trajectory generated with %d points over %.1f seconds\n\n', N_points, t_total);

%% Compute plate trajectory
fprintf('Computing plate trajectory...\n');
plate_traj = computePlateTrajectory(q_traj, time_vec);

%% Analyze sloshing
fprintf('Analyzing sloshing behavior...\n');
sloshing_analysis = analyzeSimpleSloshing(plate_traj, controller);

%% Display results
fprintf('\n=== RESULTS ===\n');
fprintf('Max plate acceleration: %.2f m/s²\n', max(vecnorm(plate_traj.acceleration, 2, 2)));
fprintf('Max sloshing amplitude: %.1f mm\n', sloshing_analysis.max_amplitude * 1000);
fprintf('Safety margin: %.1f mm\n', sloshing_analysis.safety_margin * 1000);

if sloshing_analysis.spill_risk
    fprintf('⚠️  WARNING: Spill risk detected!\n');
else
    fprintf('✅ Safe operation predicted\n');
end

%% Visualize results
fprintf('\nGenerating plots...\n');
visualizeSimpleResults(q_traj, plate_traj, sloshing_analysis, time_vec);

fprintf('\n=== Demo Complete ===\n');
fprintf('Check the generated plots for trajectory and sloshing analysis.\n');

%% Helper Functions

function q_traj = generateSimpleTrajectory(q_start, q_end, time_vec)
    % Generate smooth trajectory using polynomial interpolation
    
    N = length(time_vec);
    n_joints = length(q_start);
    q_traj = zeros(N, n_joints);
    
    % Normalize time to [0,1]
    t_norm = (time_vec - time_vec(1)) / (time_vec(end) - time_vec(1));
    
    % Use 5th order polynomial for smooth trajectory
    for i = 1:n_joints
        % Boundary conditions: position, velocity, acceleration = 0 at start/end
        % p(t) = 6*t^5 - 15*t^4 + 10*t^3 (smooth S-curve)
        s = 6*t_norm.^5 - 15*t_norm.^4 + 10*t_norm.^3;
        q_traj(:, i) = q_start(i) + (q_end(i) - q_start(i)) * s';
    end
end

function plate_traj = computePlateTrajectory(q_traj, time_vec)
    % Compute end-effector trajectory from joint trajectory
    
    N = length(time_vec);
    plate_traj = struct();
    plate_traj.time = time_vec';
    plate_traj.position = zeros(N, 3);
    plate_traj.orientation = zeros(N, 3);
    plate_traj.velocity = zeros(N, 3);
    plate_traj.acceleration = zeros(N, 3);
    
    % Simple forward kinematics
    L1 = 0.5; L2 = 0.4; L3 = 0.3; % Link lengths
    
    for i = 1:N
        q = q_traj(i, :);
        
        % Position (simplified 3-link planar robot extended to 6D)
        x = L1*cos(q(1)) + L2*cos(q(1)+q(2)) + L3*cos(q(1)+q(2)+q(3));
        y = L1*sin(q(1)) + L2*sin(q(1)+q(2)) + L3*sin(q(1)+q(2)+q(3));
        z = 0.8 + 0.2*sin(q(3)); % Height variation
        
        plate_traj.position(i, :) = [x, y, z];
        plate_traj.orientation(i, :) = [q(4), q(5), q(6)];
    end
    
    % Numerical differentiation for velocities and accelerations
    dt = time_vec(2) - time_vec(1);
    
    for i = 1:3
        plate_traj.velocity(:, i) = gradient(plate_traj.position(:, i), dt);
        plate_traj.acceleration(:, i) = gradient(plate_traj.velocity(:, i), dt);
    end
end

function analysis = analyzeSimpleSloshing(plate_traj, controller)
    % Simple sloshing analysis
    
    % Horizontal acceleration magnitude
    a_horizontal = sqrt(plate_traj.acceleration(:,1).^2 + plate_traj.acceleration(:,2).^2);
    
    % Simple sloshing model: amplitude proportional to acceleration
    sloshing_gain = 0.005; % 5mm per m/s²
    sloshing_amplitude = a_horizontal * sloshing_gain;
    
    % Analysis
    analysis = struct();
    analysis.sloshing_amplitude = sloshing_amplitude;
    analysis.max_amplitude = max(sloshing_amplitude);
    analysis.mean_amplitude = mean(sloshing_amplitude);
    analysis.safety_margin = controller.container_height - controller.filling_height - analysis.max_amplitude;
    analysis.spill_risk = analysis.safety_margin < 0.005; % 5mm threshold
    
    [~, max_idx] = max(sloshing_amplitude);
    analysis.critical_time = plate_traj.time(max_idx);
end

function visualizeSimpleResults(q_traj, plate_traj, sloshing_analysis, time_vec)
    % Create visualization plots
    
    % Joint trajectories
    figure('Name', 'Joint Trajectories', 'Position', [100, 100, 800, 600]);
    
    subplot(3,1,1);
    plot(time_vec, q_traj);
    xlabel('Time [s]');
    ylabel('Joint Position [rad]');
    title('Joint Position Trajectories');
    legend('J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'Location', 'best');
    grid on;
    
    % Joint velocities
    subplot(3,1,2);
    qd_traj = gradient(q_traj, time_vec(2) - time_vec(1));
    plot(time_vec, qd_traj);
    xlabel('Time [s]');
    ylabel('Joint Velocity [rad/s]');
    title('Joint Velocity Trajectories');
    legend('J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'Location', 'best');
    grid on;
    
    % Joint accelerations
    subplot(3,1,3);
    qdd_traj = gradient(qd_traj, time_vec(2) - time_vec(1));
    plot(time_vec, qdd_traj);
    xlabel('Time [s]');
    ylabel('Joint Acceleration [rad/s²]');
    title('Joint Acceleration Trajectories');
    legend('J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'Location', 'best');
    grid on;
    
    % End-effector trajectory
    figure('Name', 'End-Effector Analysis', 'Position', [200, 200, 800, 600]);
    
    subplot(2,2,1);
    plot3(plate_traj.position(:,1), plate_traj.position(:,2), plate_traj.position(:,3), 'b-', 'LineWidth', 2);
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    title('3D End-Effector Trajectory');
    grid on; axis equal;
    
    subplot(2,2,2);
    plot(time_vec, vecnorm(plate_traj.acceleration(:,1:2), 2, 2));
    xlabel('Time [s]');
    ylabel('Horizontal Acceleration [m/s²]');
    title('Horizontal Acceleration Magnitude');
    grid on;
    
    subplot(2,2,3);
    plot(time_vec, sloshing_analysis.sloshing_amplitude * 1000, 'r-', 'LineWidth', 2);
    hold on;
    yline(25, 'k--', 'Critical Height (25mm)', 'LineWidth', 1.5);
    yline(sloshing_analysis.safety_margin * 1000, 'g--', 'Safety Margin', 'LineWidth', 1);
    xlabel('Time [s]');
    ylabel('Sloshing Amplitude [mm]');
    title('Predicted Sloshing Amplitude');
    grid on;
    legend('Sloshing', 'Critical', 'Safety Margin', 'Location', 'best');
    
    subplot(2,2,4);
    bar([sloshing_analysis.max_amplitude * 1000, 25, sloshing_analysis.safety_margin * 1000]);
    set(gca, 'XTickLabel', {'Max Sloshing', 'Critical Limit', 'Safety Margin'});
    ylabel('Distance [mm]');
    title('Safety Analysis');
    grid on;
    
    % Color code the bars
    if sloshing_analysis.spill_risk
        bar_colors = [1, 0, 0; 0.8, 0.8, 0; 0, 1, 0]; % Red, yellow, green
    else
        bar_colors = [0, 1, 0; 0.8, 0.8, 0; 0, 1, 0]; % Green, yellow, green
    end
    
    h = bar([sloshing_analysis.max_amplitude * 1000, 25, sloshing_analysis.safety_margin * 1000]);
    h.FaceColor = 'flat';
    h.CData = bar_colors;
end