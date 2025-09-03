%% Main script to run sloshing control example
% Make sure SloshingController.m is in the MATLAB path

clc; clear; close all;

% Create controller instance
controller = SloshingController();

% Define start and end configurations (joint angles in radians)
q_start = [1.56, 0.57, -0.20, 0, -0.36, 0]';    % From paper
q_end   = [-1.56, 0.57, -0.20, 0, -0.36, 0]';   % From paper

% Set optimization options
options.max_time = 2.0;            % Maximum allowed time
options.include_sloshing = true;   % Include sloshing constraints

% Run optimization
[q_opt, t_opt, success] = controller.optimizeTrajectory(q_start, q_end, options);

% Post-processing
if success
    % Plot trajectory
    controller.plotTrajectory(q_opt, t_opt);

    % Validate trajectory
    controller.validateTrajectory(q_opt, t_opt);

    fprintf('Optimization completed successfully!\n');
    fprintf('Optimal execution time: %.3f seconds\n', t_opt);
else
    fprintf('Optimization failed!\n');
end
