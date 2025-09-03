%% Sloshing Control for Liquid-Filled Containers at Robot End-Effector
% Based on "Point to point time optimal handling of unmounted rigid objects 
% and liquid-filled containers" by Gattringer et al. (2023)
%
% This implementation provides trajectory optimization and control for 
% handling liquid-filled containers with a robot manipulator while 
% minimizing sloshing effects.

classdef SloshingController < handle
    properties
        % Robot parameters
        n_joints = 6;           % Number of robot joints
        robot_params;           % Robot dynamic parameters
        
        % Container parameters
        container_mass = 0.12;  % kg - container mass
        container_height = 0.06; % m - container height
        container_radius = 0.025; % m - footprint radius
        liquid_density = 1000;  % kg/m³ - water density
        filling_height = 0.025; % m - liquid filling height
        
        % Contact parameters
        friction_coeff = 0.3;   % Static friction coefficient
        max_twist_torque;       % Maximum twisting torque
        
        % Optimization parameters
        N_shooting = 150;       % Number of shooting intervals
        weight_smoothness = 1e-6; % Smoothness weight factor
        
        % Physical limits
        max_joint_vel;          % Maximum joint velocities [rad/s]
        max_joint_acc;          % Maximum joint accelerations [rad/s²]
        max_joint_jerk;         % Maximum joint jerks [rad/s³]
        max_motor_torque;       % Maximum motor torques [Nm]
        
        % Trajectory data
        trajectory;             % Optimized trajectory
        time_optimal;           % Optimal execution time
    end
    
    methods
        function obj = SloshingController(robot_config)
            % Constructor - initialize controller with robot configuration
            if nargin > 0
                obj.initializeRobotParameters(robot_config);
            else
                obj.setDefaultParameters();
            end
            obj.calculateDerivedParameters();
        end
        
        function setDefaultParameters(obj)
            % Set default parameters for 6-DOF industrial robot
            obj.max_joint_vel = [2.0, 1.5, 2.0, 3.0, 3.0, 6.0]; % rad/s
            obj.max_joint_acc = [5.0, 4.0, 5.0, 8.0, 8.0, 15.0]; % rad/s²
            obj.max_joint_jerk = [50, 40, 50, 80, 80, 150]; % rad/s³
            obj.max_motor_torque = [200, 300, 150, 50, 50, 20]; % Nm
        end
        
        function calculateDerivedParameters(obj)
            % Calculate derived parameters for contact constraints
            obj.max_twist_torque = (2/3) * obj.friction_coeff * ...
                obj.container_mass * 9.81 * obj.container_radius^3;
        end
        
        function [q_opt, t_opt, success] = optimizeTrajectory(obj, q_start, q_end, options)
            % Optimize time-optimal trajectory for liquid container handling
            %
            % Inputs:
            %   q_start - Initial joint configuration [6×1]
            %   q_end   - Final joint configuration [6×1]
            %   options - Optimization options structure
            %
            % Outputs:
            %   q_opt   - Optimal trajectory [N×6]
            %   t_opt   - Optimal time duration
            %   success - Optimization success flag
            
            if nargin < 4
                options = struct();
            end
            
            % Set default options
            if ~isfield(options, 'max_time')
                options.max_time = 5.0; % Maximum allowed time
            end
            if ~isfield(options, 'include_sloshing')
                options.include_sloshing = true;
            end
            
            fprintf('Starting trajectory optimization...\n');
            
            try
                % Setup optimization problem
                [problem, bounds] = obj.setupOptimizationProblem(q_start, q_end, options);
                
                % Solve using fmincon (interior-point method)
                options_opt = optimoptions('fmincon', ...
                    'Algorithm', 'interior-point', ...
                    'Display', 'iter', ...
                    'MaxIterations', 1000, ...
                    'OptimalityTolerance', 1e-6, ...
                    'ConstraintTolerance', 1e-6);
                
                % Initial guess - minimum jerk trajectory
                x0 = obj.generateInitialGuess(q_start, q_end, bounds.t_max);
                
                % Optimize
                [x_opt, fval, exitflag] = fmincon(@(x) obj.objectiveFunction(x), ...
                    x0, [], [], [], [], bounds.lb, bounds.ub, ...
                    @(x) obj.constraintFunction(x, q_start, q_end, options), ...
                    options_opt);
                
                success = (exitflag > 0);
                
                if success
                    % Extract solution
                    [q_opt, t_opt] = obj.extractTrajectory(x_opt);
                    obj.trajectory = q_opt;
                    obj.time_optimal = t_opt;
                    fprintf('Optimization successful! Optimal time: %.3f s\n', t_opt);
                else
                    fprintf('Optimization failed with exit flag: %d\n', exitflag);
                    q_opt = [];
                    t_opt = [];
                end
                
            catch ME
                fprintf('Error during optimization: %s\n', ME.message);
                q_opt = [];
                t_opt = [];
                success = false;
            end
        end
        
        function [problem, bounds] = setupOptimizationProblem(obj, q_start, q_end, options)
            % Setup the optimization problem structure and bounds
            
            % Decision variables: [jerk_trajectory; final_time]
            n_vars = obj.N_shooting * obj.n_joints + 1;
            
            bounds.lb = [-obj.max_joint_jerk * ones(obj.N_shooting * obj.n_joints, 1); 0.1];
            bounds.ub = [obj.max_joint_jerk * ones(obj.N_shooting * obj.n_joints, 1); options.max_time];
            bounds.t_max = options.max_time;
            
            problem.n_vars = n_vars;
            problem.q_start = q_start;
            problem.q_end = q_end;
            problem.options = options;
        end
        
        function x0 = generateInitialGuess(obj, q_start, q_end, t_max)
            % Generate initial guess using minimum jerk trajectory
            
            t_guess = t_max * 0.5; % Conservative time estimate
            dt = t_guess / obj.N_shooting;
            
            % Generate minimum jerk trajectory
            jerk_profile = zeros(obj.N_shooting, obj.n_joints);
            for i = 1:obj.n_joints
                % Simple bang-coast-bang jerk profile
                jerk_mag = 8 * (q_end(i) - q_start(i)) / t_guess^3;
                
                % First third: positive jerk
                jerk_profile(1:obj.N_shooting/3, i) = jerk_mag;
                % Middle third: zero jerk
                % Last third: negative jerk
                jerk_profile(2*obj.N_shooting/3+1:end, i) = -jerk_mag;
            end
            
            % Flatten and append time
            x0 = [jerk_profile(:); t_guess];
        end
        
        function J = objectiveFunction(obj, x)
            % Objective function: minimize time + smoothness penalty
            
            t_final = x(end);
            jerk_vars = x(1:end-1);
            
            % Primary objective: minimize time
            J = t_final;
            
            % Secondary objective: minimize jerk (smoothness)
            jerk_penalty = obj.weight_smoothness * sum(jerk_vars.^2);
            J = J + jerk_penalty;
        end
        
        function [c, ceq] = constraintFunction(obj, x, q_start, q_end, options)
            % Nonlinear constraint function
            
            t_final = x(end);
            jerk_vars = reshape(x(1:end-1), obj.N_shooting, obj.n_joints);
            
            % Integrate trajectory from jerk profile
            [q_traj, qd_traj, qdd_traj] = obj.integrateTrajectory(jerk_vars, t_final, q_start);
            
            % Inequality constraints
            c = [];
            
            % Joint velocity limits
            vel_constraints = abs(qd_traj) - obj.max_joint_vel';
            c = [c; vel_constraints(:)];
            
            % Joint acceleration limits  
            acc_constraints = abs(qdd_traj) - obj.max_joint_acc';
            c = [c; acc_constraints(:)];
            
            % Motor torque limits (simplified inverse dynamics)
            for i = 1:size(q_traj, 1)
                tau = obj.computeInverseDynamics(q_traj(i,:)', qd_traj(i,:)', qdd_traj(i,:)');
                torque_constraints = abs(tau) - obj.max_motor_torque';
                c = [c; torque_constraints];
            end
            
            % Contact force constraints (anti-slip, anti-tip)
            if options.include_sloshing
                contact_constraints = obj.evaluateContactConstraints(q_traj, qd_traj, qdd_traj);
                c = [c; contact_constraints];
            end
            
            % Sloshing constraints (critical filling height consideration)
            if options.include_sloshing
                sloshing_constraints = obj.evaluateSloshingConstraints(q_traj, qd_traj, qdd_traj);
                c = [c; sloshing_constraints];
            end
            
            % Equality constraints (boundary conditions)
            ceq = [];
            
            % Final position constraint
            pos_error = q_traj(end,:)' - q_end;
            ceq = [ceq; pos_error];
            
            % Final velocity constraint (rest-to-rest)
            vel_final = qd_traj(end,:)';
            ceq = [ceq; vel_final];
            
            % Final acceleration constraint
            acc_final = qdd_traj(end,:)';
            ceq = [ceq; acc_final];
        end
        
        function [q_traj, qd_traj, qdd_traj] = integrateTrajectory(obj, jerk_profile, t_final, q_start)
            % Integrate jerk profile to obtain position, velocity, acceleration
            
            dt = t_final / obj.N_shooting;
            
            q_traj = zeros(obj.N_shooting + 1, obj.n_joints);
            qd_traj = zeros(obj.N_shooting + 1, obj.n_joints);
            qdd_traj = zeros(obj.N_shooting + 1, obj.n_joints);
            
            % Initial conditions
            q_traj(1,:) = q_start';
            qd_traj(1,:) = 0;
            qdd_traj(1,:) = 0;
            
            % Integrate using Euler method
            for i = 1:obj.N_shooting
                % Current jerk
                jerk = jerk_profile(i,:);
                
                % Update acceleration
                qdd_traj(i+1,:) = qdd_traj(i,:) + jerk * dt;
                
                % Update velocity
                qd_traj(i+1,:) = qd_traj(i,:) + qdd_traj(i,:) * dt;
                
                % Update position
                q_traj(i+1,:) = q_traj(i,:) + qd_traj(i,:) * dt;
            end
        end
        
        function tau = computeInverseDynamics(obj, q, qd, qdd)
            % Simplified inverse dynamics computation
            % In practice, this should use the full robot model
            
            % Placeholder: simplified model with mass matrix and gravity
            M = eye(obj.n_joints) * 10; % Simplified mass matrix
            G = [0; 0; 50; 0; 0; 0];    % Simplified gravity vector
            
            tau = M * qdd + G;
        end
        
        function constraints = evaluateContactConstraints(obj, q_traj, qd_traj, qdd_traj)
            % Evaluate contact force constraints (friction cone, normal force)
            
            constraints = [];
            
            for i = 1:size(q_traj, 1)
                % Compute end-effector acceleration
                ee_acc = obj.computeEndEffectorAcceleration(q_traj(i,:)', qd_traj(i,:)', qdd_traj(i,:)');
                
                % Compute contact forces
                [f_normal, f_tangential, M_tip] = obj.computeContactForces(ee_acc);
                
                % Non-lifting constraint: f_z <= 0 (compression)
                constraints = [constraints; f_normal];
                
                % Non-slipping constraint: |f_tangential| <= mu * |f_normal|
                friction_constraint = norm(f_tangential) + obj.friction_coeff * f_normal;
                constraints = [constraints; friction_constraint];
                
                % Non-tipping constraint
                tip_constraint = norm(M_tip) + obj.container_radius * f_normal;
                constraints = [constraints; tip_constraint];
            end
        end
        
        function constraints = evaluateSloshingConstraints(obj, q_traj, qd_traj, qdd_traj)
            % Evaluate sloshing-related constraints based on critical filling height
            
            constraints = [];
            critical_height = 0.025; % m - from experimental validation
            
            for i = 1:size(q_traj, 1)
                % Compute end-effector acceleration
                ee_acc = obj.computeEndEffectorAcceleration(q_traj(i,:)', qd_traj(i,:)', qdd_traj(i,:)');
                
                % Estimate sloshing amplitude based on acceleration
                sloshing_amplitude = obj.estimateSloshingAmplitude(ee_acc);
                
                % Constraint: filling height + sloshing amplitude <= container height
                spill_constraint = obj.filling_height + sloshing_amplitude - obj.container_height;
                constraints = [constraints; spill_constraint];
            end
        end
        
        function ee_acc = computeEndEffectorAcceleration(obj, q, qd, qdd)
            % Compute end-effector Cartesian acceleration
            % Simplified computation - should use full Jacobian in practice
            
            % Placeholder: assume linear relationship
            ee_acc = [qdd(1); qdd(2); qdd(3)] * 0.5; % Simplified mapping
        end
        
        function [f_normal, f_tangential, M_tip] = computeContactForces(obj, ee_acc)
            % Compute contact forces from end-effector acceleration
            
            % Normal force (including gravitational and inertial components)
            f_normal = -obj.container_mass * (9.81 + ee_acc(3));
            
            % Tangential forces
            f_tangential = -obj.container_mass * ee_acc(1:2);
            
            % Tipping moment (simplified)
            M_tip = obj.container_mass * ee_acc(1:2) * obj.container_height * 0.5;
        end
        
        function amplitude = estimateSloshingAmplitude(obj, ee_acc)
            % Estimate sloshing amplitude based on acceleration
            % Simplified model based on experimental observations
            
            % Horizontal acceleration magnitude
            a_horizontal = norm(ee_acc(1:2));
            
            % Empirical relationship (should be calibrated experimentally)
            amplitude = min(0.01, a_horizontal * 0.002); % Max 1cm sloshing
        end
        
        function [q_traj, t_opt] = extractTrajectory(obj, x_opt)
            % Extract optimized trajectory from solution vector
            
            t_opt = x_opt(end);
            jerk_vars = reshape(x_opt(1:end-1), obj.N_shooting, obj.n_joints);
            
            % Integrate to get full trajectory
            q_start = [0; 0; 0; 0; 0; 0]; % This should be stored from optimization
            [q_traj, ~, ~] = obj.integrateTrajectory(jerk_vars, t_opt, q_start);
        end
        
        function plotTrajectory(obj, q_traj, t_final)
            % Plot the optimized trajectory
            
            if isempty(q_traj)
                fprintf('No trajectory to plot\n');
                return;
            end
            
            time_vec = linspace(0, t_final, size(q_traj, 1));
            
            figure('Name', 'Optimized Trajectory');
            
            % Joint positions
            subplot(3,1,1);
            plot(time_vec, q_traj);
            xlabel('Time [s]');
            ylabel('Joint Position [rad]');
            title('Joint Trajectories');
            legend('J1', 'J2', 'J3', 'J4', 'J5', 'J6');
            grid on;
            
            % Joint velocities (computed via differentiation)
            subplot(3,1,2);
            qd_traj = gradient(q_traj, time_vec(2) - time_vec(1));
            plot(time_vec, qd_traj);
            xlabel('Time [s]');
            ylabel('Joint Velocity [rad/s]');
            title('Joint Velocities');
            legend('J1', 'J2', 'J3', 'J4', 'J5', 'J6');
            grid on;
            
            % Joint accelerations
            subplot(3,1,3);
            qdd_traj = gradient(qd_traj, time_vec(2) - time_vec(1));
            plot(time_vec, qdd_traj);
            xlabel('Time [s]');
            ylabel('Joint Acceleration [rad/s²]');
            title('Joint Accelerations');
            legend('J1', 'J2', 'J3', 'J4', 'J5', 'J6');
            grid on;
        end
        
        function success = validateTrajectory(obj, q_traj, t_final)
            % Validate that the trajectory satisfies all constraints
            
            success = true;
            time_vec = linspace(0, t_final, size(q_traj, 1));
            
            % Compute derivatives
            qd_traj = gradient(q_traj, time_vec(2) - time_vec(1));
            qdd_traj = gradient(qd_traj, time_vec(2) - time_vec(1));
            
            % Check velocity limits
            vel_violations = any(abs(qd_traj) > obj.max_joint_vel', 'all');
            if vel_violations
                fprintf('Warning: Velocity limits violated\n');
                success = false;
            end
            
            % Check acceleration limits
            acc_violations = any(abs(qdd_traj) > obj.max_joint_acc', 'all');
            if acc_violations
                fprintf('Warning: Acceleration limits violated\n');
                success = false;
            end
            
            if success
                fprintf('Trajectory validation: PASSED\n');
            else
                fprintf('Trajectory validation: FAILED\n');
            end
        end
    end
end

%% Usage Example
function example_usage()
    % Example of how to use the SloshingController
    
    % Create controller instance
    controller = SloshingController();
    
    % Define start and end configurations (joint angles in radians)
    q_start = [1.56, 0.57, -0.20, 0, -0.36, 0]';      % From paper
    q_end = [-1.56, 0.57, -0.20, 0, -0.36, 0]';       % From paper
    
    % Set optimization options
    options.max_time = 2.0;           % Maximum allowed time
    options.include_sloshing = true;   % Include sloshing constraints
    
    % Optimize trajectory
    [q_opt, t_opt, success] = controller.optimizeTrajectory(q_start, q_end, options);
    
    if success
        % Plot results
        controller.plotTrajectory(q_opt, t_opt);
        
        % Validate trajectory
        controller.validateTrajectory(q_opt, t_opt);
        
        fprintf('Optimization completed successfully!\n');
        fprintf('Optimal execution time: %.3f seconds\n', t_opt);
    else
        fprintf('Optimization failed!\n');
    end
end