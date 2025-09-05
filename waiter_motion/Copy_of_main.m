    clear
    close all 
    clc
   
    % Create controller instance with perfect kinematics assumption
    controller = Copy_of_SloshingController();
    
    % Define task parameters
    task_params = defineTaskParameters();
    
    % Generate multiple test scenarios
    scenarios = generateTestScenarios();
    
    % Run optimization for each scenario
    results = cell(length(scenarios), 1);
    
    for i = 1:length(scenarios)
        fprintf('--- Scenario %d: %s ---\n', i, scenarios{i}.name);
        
        % Set current scenario parameters
        controller.filling_height = scenarios{i}.filling_height;
        controller.container_mass = scenarios{i}.container_mass;
        
        % Optimize trajectory
        [q_opt, t_opt, success] = controller.optimizeTrajectory(...
            scenarios{i}.q_start, scenarios{i}.q_end, scenarios{i}.options);
        
        % Store results
        results{i} = struct('scenario', scenarios{i}, 'q_opt', q_opt, ...
                           't_opt', t_opt, 'success', success);
        
        if success
            % Generate plate trajectory for sloshing compensation
            plate_traj = generatePlateTrajectory(controller, q_opt, t_opt);
            results{i}.plate_trajectory = plate_traj;
            
            % Analyze sloshing behavior
            sloshing_analysis = analyzeSloshingBehavior(controller, q_opt, t_opt);
            results{i}.sloshing_analysis = sloshing_analysis;
   
            fprintf('  ✓ Success! Optimal time: %.3f s\n', t_opt);
            fprintf('  ✓ Max sloshing amplitude: %.1f mm\n', sloshing_analysis.max_amplitude * 1000);
            fprintf('  ✓ Safety margin: %.1f mm\n', sloshing_analysis.safety_margin * 1000);
        else
            fprintf('  ✗ Optimization failed\n');
        end
        fprintf('\n');
    end
    
    % Visualize results
    visualizeResults(results);
    
    % Generate comparison report
    generateComparisonReport(results)

    controller.plotTrajectory(q_opt, t_opt);
    controller.plotNormalizedVelocities(q_opt, t_opt);

function task_params = defineTaskParameters()
    % Define the task parameters for liquid container handling
    
    task_params = struct();
    
    % Workspace limits (assuming perfect kinematics)
    task_params.workspace_limits = struct(...
        'x', [-2.0, 2.0], ...  % meters
        'y', [-2.0, 2.0], ...
        'z', [0.2, 2.0], ...
        'roll', [-pi, pi], ... % radians
        'pitch', [-pi/2, pi/2], ...
        'yaw', [-pi, pi]);
    
    % Container specifications
    task_params.container = struct(...
        'type', 'coffee_cup', ...
        'lower_diameter', 0.046, ...  % meters
        'upper_diameter', 0.072, ...
        'height', 0.088, ...
        'wall_thickness', 0.002);
    
    % Liquid properties
    task_params.liquid = struct(...
        'density', 1000, ...      % kg/m³ (water)
        'viscosity', 1e-6, ...    % m²/s (kinematic)
        'surface_tension', 0.072); % N/m
end

function scenarios = generateTestScenarios()
    % Generate different test scenarios for validation
    
    scenarios = {};
    
    % Scenario 1: Conservative filling (safe case)
    scenarios{1} = struct(...
        'name', 'Conservative Filling (20mm)', ...
        'filling_height', 0.020, ...  % 20mm - below critical
        'container_mass', 0.12, ...   % kg
        'q_start', [0, 0, 0, 0, 0, 0]', ...
        'q_end', [pi/2, 0, 0, 0, 0, 0]', ...
        'options', struct('max_time', 3.0, 'include_sloshing', true));

    % % Scenario 2: Critical filling (challenging case)
    % scenarios{2} = struct(...
    %     'name', 'Critical Filling (25mm)', ...
    %     'filling_height', 0.025, ...  % 25mm - at critical limit
    %     'container_mass', 0.14, ...   % kg (heavier with more liquid)
    %     'q_start', [pi/6, 0, 0, 0, 0, 0]', ...
    %     'q_end', [-pi/6, 0, 0, 0, 0, pi/4]', ...
    %     'options', struct('max_time', 4.0, 'include_sloshing', true));
    % 
    % % Scenario 3: High filling (requires careful trajectory)
    % scenarios{3} = struct(...
    %     'name', 'High Filling (30mm)', ...
    %     'filling_height', 0.030, ...  % 30mm - above critical
    %     'container_mass', 0.16, ...   % kg
    %     'q_start', [0, pi/4, 0, 0, 0, 0]', ...
    %     'q_end', [pi/3, -pi/6, 0, 0, 0, pi/2]', ...
    %     'options', struct('max_time', 5.0, 'include_sloshing', true));
    % 
    % % Scenario 4: Fast motion requirement (time-critical)
    % scenarios{4} = struct(...
    %     'name', 'Fast Motion (22mm)', ...
    %     'filling_height', 0.022, ...  % 22mm
    %     'container_mass', 0.13, ...   % kg
    %     'q_start', [-pi/4, 0, 0, 0, 0, 0]', ...
    %     'q_end', [pi/4, 0, 0, 0, 0, -pi/3]', ...
    %     'options', struct('max_time', 1.5, 'include_sloshing', true)); % Aggressive time limit
end

function visualizeResults(results)
    % Visualize optimization results and sloshing analysis
    
    % Create comprehensive visualization
    figure('Name', 'Sloshing Compensation Results', 'Position', [100, 100, 1200, 800]);
    
    n_scenarios = length(results);
    
    for i = 1:n_scenarios
        if ~results{i}.success
            continue;
        end
        
        % Plot trajectory and sloshing analysis
        subplot(3, n_scenarios, i);
        if ~isempty(results{i}.plate_trajectory)
            plot3(results{i}.plate_trajectory.position(:,1), ...
                  results{i}.plate_trajectory.position(:,2), ...
                  results{i}.plate_trajectory.position(:,3), 'b-', 'LineWidth', 2);
            xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
            title(sprintf('Scenario %d: Plate Trajectory', i));
            grid on; axis equal;
        end
        
        % Plot sloshing amplitude over time
        subplot(3, n_scenarios, i + n_scenarios);
        if ~isempty(results{i}.sloshing_analysis)
            plot(results{i}.sloshing_analysis.time, ...
                 results{i}.sloshing_analysis.sloshing_amplitude * 1000, 'r-', 'LineWidth', 2);
            hold on;
            yline(25, 'k--', 'Critical Height', 'LineWidth', 1.5);
            xlabel('Time [s]');
            ylabel('Sloshing Amplitude [mm]');
            title(sprintf('Sloshing Analysis - Scenario %d', i));
            grid on;
        end

                % Plot sloshing amplitude over time
        subplot(3, n_scenarios, i + 2*n_scenarios);
        if ~isempty(results{i}.sloshing_analysis)
            plot(results{i}.plate_trajectory.time, ...
                 results{i}.plate_trajectory.velocity, 'r-', 'LineWidth', 2);
            hold on;
            yline(25, 'k--', 'Critical Height', 'LineWidth', 1.5);
            xlabel('Time [s]');
            ylabel('Sloshing Amplitude [mm]');
            title(sprintf('Sloshing Analysis - Scenario %d', i));
            legend();
            grid on;
        end
    end
    
    % Create summary comparison plot
    figure('Name', 'Performance Comparison', 'Position', [200, 200, 800, 600]);
    
    % Extract performance metrics
    exec_times = [];
    max_sloshing = [];
    scenario_names = {};
    
    for i = 1:n_scenarios
        if results{i}.success
            exec_times(end+1) = results{i}.t_opt;
            max_sloshing(end+1) = results{i}.sloshing_analysis.max_amplitude * 1000;
            scenario_names{end+1} = results{i}.scenario.name;
        end
    end
    
    % Performance comparison plots
    subplot(2,1,1);
    bar(exec_times);
    set(gca, 'XTickLabel', scenario_names);
    xlabel('Scenarios');
    ylabel('Execution Time [s]');
    title('Trajectory Execution Time Comparison');
    xtickangle(45);
    
    subplot(2,1,2);
    bar(max_sloshing);
    hold on;
    yline(25, 'r--', 'Critical Threshold', 'LineWidth', 2);
    set(gca, 'XTickLabel', scenario_names);
    xlabel('Scenarios');
    ylabel('Max Sloshing [mm]');
    title('Maximum Sloshing Amplitude Comparison');
    xtickangle(45);
end

function generateComparisonReport(results)
    % Generate a detailed comparison report
    
    fprintf('\n=== SLOSHING COMPENSATION PERFORMANCE REPORT ===\n\n');
    
    successful_scenarios = 0;
    total_scenarios = length(results);
    
    for i = 1:total_scenarios
        fprintf('Scenario %d: %s\n', i, results{i}.scenario.name);
        fprintf('  Filling Height: %.1f mm\n', results{i}.scenario.filling_height * 1000);
        
        if results{i}.success
            successful_scenarios = successful_scenarios + 1;
            fprintf('  Status: SUCCESS ✓\n');
            fprintf('  Execution Time: %.3f s\n', results{i}.t_opt);
            fprintf('  Max Sloshing: %.1f mm\n', results{i}.sloshing_analysis.max_amplitude * 1000);
            fprintf('  Safety Margin: %.1f mm\n', results{i}.sloshing_analysis.safety_margin * 1000);
            
            if results{i}.sloshing_analysis.spill_risk
                fprintf('  ⚠️  WARNING: High spill risk detected!\n');
            else
                fprintf('  ✓ Safe operation confirmed\n');
            end
        else
            fprintf('  Status: FAILED ✗\n');
            fprintf('  Reason: Optimization constraints could not be satisfied\n');
        end
        fprintf('\n');
    end
    
    fprintf('=== SUMMARY ===\n');
    fprintf('Successful scenarios: %d/%d (%.1f%%)\n', ...
        successful_scenarios, total_scenarios, ...
        100 * successful_scenarios / total_scenarios);
    
    if successful_scenarios > 0
        % Find best performing scenario
        best_time = inf;
        best_scenario_idx = 0;
        for i = 1:total_scenarios
            if results{i}.success && results{i}.t_opt < best_time
                best_time = results{i}.t_opt;
                best_scenario_idx = i;
            end
        end
        
        fprintf('Fastest execution: %.3f s (Scenario %d: %s)\n', ...
            best_time, best_scenario_idx, results{best_scenario_idx}.scenario.name);
    end
    
    fprintf('\n=== RECOMMENDATIONS ===\n');
    fprintf('1. For filling heights >25mm, increase execution time budget\n');
    fprintf('2. Consider adaptive tilt compensation for high accelerations\n');
    fprintf('3. Implement real-time sloshing monitoring for safety\n');
    fprintf('4. Validate results with physical experiments\n');
end
function plate_traj = generatePlateTrajectory(controller, q_traj, t_final)
    % Generate the plate trajectory with sloshing compensation
    
    if isempty(q_traj)
        plate_traj = [];
        return;
    end
    
    time_vec = linspace(0, t_final, size(q_traj, 1));
    
    % Initialize plate trajectory structure
    plate_traj = struct();
    plate_traj.time = time_vec';
    plate_traj.position = zeros(length(time_vec), 3);    % [x, y, z]
    plate_traj.orientation = zeros(length(time_vec), 3); % [roll, pitch, yaw]
    plate_traj.velocity = zeros(length(time_vec), 3);
    plate_traj.acceleration = zeros(length(time_vec), 3);
    
    % Compute end-effector pose for each time step (simplified forward kinematics)
    for i = 1:length(time_vec)
        % Forward kinematics (simplified - assumes perfect robot)
        [pos, orient] = forwardKinematics(q_traj(i, :));
        plate_traj.position(i, :) = pos;
        plate_traj.orientation(i, :) = orient;
    end
    
    % Compute velocities and accelerations via numerical differentiation
    dt = time_vec(2) - time_vec(1);
    plate_traj.velocity = gradient(plate_traj.position, dt);
    plate_traj.acceleration = gradient(plate_traj.velocity, dt);
    
    % Apply sloshing compensation to orientation
    plate_traj = applySloshingCompensation(plate_traj, controller);
end
function [position, orientation] = forwardKinematics(q)
    % Simplified forward kinematics for perfect robot
    % Returns end-effector position and orientation
    
    % Link lengths (example values)
    L1 = 0.5; L2 = 0.4; L3 = 0.3;
    
    % Simplified 6-DOF forward kinematics
    x = L1*cos(q(1)) + L2*cos(q(1)+q(2)) + L3*cos(q(1)+q(2)+q(3));
    y = L1*sin(q(1)) + L2*sin(q(1)+q(2)) + L3*sin(q(1)+q(2)+q(3));
    z = 0.5 + 0.1*q(3); % Simplified vertical component
    
    position = [x, y, z];
    orientation = [q(4), q(5), q(6)]; % Direct mapping for simplicity
end
function plate_traj = applySloshingCompensation(plate_traj, controller)
    % Apply orientation compensation to reduce sloshing
    
    % Compensation parameters
    compensation_gain = 0.1;  % Tuning parameter
    max_tilt_angle = 0.1;     % Maximum allowed tilt (rad)
    
    for i = 1:length(plate_traj.time)
        % Get horizontal acceleration
        a_horizontal = plate_traj.acceleration(i, 1:2);
        
        % Calculate compensation tilt angles
        % Tilt opposite to acceleration direction to counteract liquid motion
        tilt_x = -compensation_gain * a_horizontal(1);
        tilt_y = -compensation_gain * a_horizontal(2);
        
        % Limit tilt angles
        tilt_x = max(-max_tilt_angle, min(max_tilt_angle, tilt_x));
        tilt_y = max(-max_tilt_angle, min(max_tilt_angle, tilt_y));
        
        % Apply compensation to plate orientation
        plate_traj.orientation(i, 1) = plate_traj.orientation(i, 1) + tilt_x; % roll
        plate_traj.orientation(i, 2) = plate_traj.orientation(i, 2) + tilt_y; % pitch
    end
    
    % Store compensation info
    plate_traj.compensation = struct(...
        'gain', compensation_gain, ...
        'max_tilt', max_tilt_angle, ...
        'applied', true);
end

function analysis = analyzeSloshingBehavior(controller, q_traj, t_final)
    % Analyze expected sloshing behavior along the trajectory
    
    if isempty(q_traj)
        analysis = [];
        return;
    end
    
    time_vec = linspace(0, t_final, size(q_traj, 1));
    sloshing_amplitude = zeros(size(time_vec));
    
    % Compute derivatives for acceleration analysis
    qd_traj = gradient(q_traj, time_vec(2) - time_vec(1));
    qdd_traj = gradient(qd_traj, time_vec(2) - time_vec(1));
    
    for i = 1:length(time_vec)
        % Compute end-effector acceleration
        ee_acc = controller.computeEndEffectorAcceleration(...
            q_traj(i,:)', qd_traj(i,:)', qdd_traj(i,:)');
        
        % Estimate sloshing amplitude
        sloshing_amplitude(i) = controller.estimateSloshingAmplitude(ee_acc);
    end

    function plotNormalizedVelocities(obj, q_traj, t_final)
            % Plot normalized joint velocities (as fraction of max limits)
            
            if isempty(q_traj)
                fprintf('No trajectory to plot\n');
                return;
            end
            
            time_vec = linspace(0, t_final, size(q_traj, 1));
            
            % Compute velocities via differentiation
            qd_traj = gradient(q_traj, time_vec(2) - time_vec(1));
            
            % Normalize velocities by their respective limits
            qd_normalized = zeros(size(qd_traj));
            for i = 1:obj.n_joints
                qd_normalized(:,i) = qd_traj(:,i) / obj.max_joint_vel(i);
            end
            
            figure('Name', 'Normalized Joint Velocities');
            
            % Plot normalized velocities
            plot(time_vec, qd_normalized, 'LineWidth', 1.5);
            xlabel('Time [s]');
            ylabel('Normalized Velocity (fraction of limit)');
            title('Normalized Joint Velocities');
            legend('J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'Location', 'best');
            grid on;
            
            % Add reference lines
            hold on;
            plot([time_vec(1), time_vec(end)], [1, 1], 'r--', 'LineWidth', 2, 'DisplayName', '+100% limit');
            plot([time_vec(1), time_vec(end)], [-1, -1], 'r--', 'LineWidth', 2, 'DisplayName', '-100% limit');
            plot([time_vec(1), time_vec(end)], [1.5, 1.5], 'k:', 'LineWidth', 1, 'DisplayName', 'Constraint limit (150%)');
            plot([time_vec(1), time_vec(end)], [-1.5, -1.5], 'k:', 'LineWidth', 1, 'DisplayName', 'Constraint limit (-150%)');
            hold off;
            
            % Set y-axis limits for better visualization
            ylim([-2, 2]);
            
            % Print maximum normalized velocities
            fprintf('\nMaximum normalized velocities:\n');
            for i = 1:obj.n_joints
                max_norm_vel = max(abs(qd_normalized(:,i)));
                fprintf('  Joint %d: %.2f%% of limit (%.3f/%.3f rad/s)\n', ...
                    i, max_norm_vel*100, max(abs(qd_traj(:,i))), obj.max_joint_vel(i));
            end
            
            % Check for violations
            violations = any(abs(qd_normalized) > 1.5, 1);
            if any(violations)
                fprintf('\nVelocity constraint violations (>150%% limit):\n');
                for i = find(violations)
                    fprintf('  Joint %d: VIOLATION\n', i);
                end
            else
                fprintf('\nAll joints within velocity constraints (150%% limit)\n');
            end
        end
    
    % Analysis results
    analysis = struct();
    analysis.time = time_vec';
    analysis.sloshing_amplitude = sloshing_amplitude';
    analysis.max_amplitude = max(sloshing_amplitude);
    analysis.mean_amplitude = mean(sloshing_amplitude);
    analysis.safety_margin = controller.container_height - controller.filling_height - analysis.max_amplitude;
    analysis.spill_risk = analysis.safety_margin < 0.005; % 5mm safety threshold
    
    % Find critical time points
    [~, max_idx] = max(sloshing_amplitude);
    analysis.critical_time = time_vec(max_idx);
    analysis.critical_acceleration = controller.computeEndEffectorAcceleration(...
        q_traj(max_idx,:)', qd_traj(max_idx,:)', qdd_traj(max_idx,:)');
end