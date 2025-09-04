
fprintf('=== Sloshing Compensation Controller for Robot Arm ===\n\n');

controller = SloshingController();
task_params = defineTaskParameters();
scenarios = generateTestScenarios();

results = cell(length(scenarios), 1);

for i = 1:length(scenarios)
    fprintf('--- Scenario %d: %s ---\n', i, scenarios{i}.name);
    
    controller.filling_height = scenarios{i}.filling_height;
    controller.container_mass = scenarios{i}.container_mass;
    [q_opt, t_opt, success] = controller.optimizeTrajectory(...
        scenarios{i}.q_start, scenarios{i}.q_end, scenarios{i}.options);
    results{i} = struct('scenario', scenarios{i}, 'q_opt', q_opt, ...
                       't_opt', t_opt, 'success', success);
    
    if success
        plate_traj = generatePlateTrajectory(controller, q_opt, t_opt);
        results{i}.plate_trajectory = plate_traj;
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

function task_params = defineTaskParameters()
    task_params = struct();
    task_params.workspace_limits = struct(...
        'x', [-2.0, 2.0], ...  % meters
        'y', [-2.0, 2.0], ...
        'z', [0.2, 2.0], ...
        'roll', [-pi, pi], ... % radians
        'pitch', [-pi/2, pi/2], ...
        'yaw', [-pi, pi]);

    task_params.container = struct(...
        'type', 'coffee_cup', ...
        'lower_diameter', 0.046, ...  % meters
        'upper_diameter', 0.072, ...
        'height', 0.088, ...
        'wall_thickness', 0.002);

    task_params.liquid = struct(...
        'density', 1000, ...      % kg/m³ (water)
        'viscosity', 1e-6, ...    % m²/s (kinematic)
        'surface_tension', 0.072); % N/m
end

function scenarios = generateTestScenarios()
    scenarios = {};
    scenarios{1} = struct(...
        'name', 'Conservative Filling (20mm)', ...
        'filling_height', 0.020, ...  % 20mm - below critical
        'container_mass', 0.12, ...   % kg
        'q_start', [0, 0, 0, 0, 0, 0]', ...
        'q_end', [pi/2, 0, 0, 0, 0, 0]', ...
        'options', struct('debug', true, 'max_time', 3.0, 'include_sloshing', true));
end