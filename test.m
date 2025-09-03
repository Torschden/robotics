function test()
    % Nonlinear Sloshing Dynamics Simulation
    % Based on spring-damper approximation with coupled equations of motion
    
    clear; clc; close all;
    
    %% Parameters
    % Physical parameters
    omega_n = 2*pi*0.5;  % Natural frequency (rad/s)
    zeta_n = 0.05;       % Damping ratio
    g = 9.81;            % Gravitational acceleration (m/s^2)
    R = 1.0;             % Tank radius (m)
    alpha_n = 1.0;       % Nonlinear coefficient
    
    % Simulation parameters
    t_end = 20;          % Simulation time (s)
    dt = 0.001;          % Time step (s)
    t = 0:dt:t_end;      % Time vector
    
    %% Excitation Definition
    % Define arbitrary excitation (you can modify this)
    excitation_type = 'random'; % Options: 'harmonic', 'chirp', 'step', 'random'
    
    switch excitation_type
        case 'harmonic'
            freq = 0.8;  % Excitation frequency (Hz)
            amp = 0.1;   % Amplitude (m)
            x0_ddot = amp * (2*pi*freq)^2 * sin(2*pi*freq*t);
            y0_ddot = 0.5 * amp * (2*pi*freq)^2 * cos(2*pi*freq*t);
            
        case 'chirp'
            f0 = 0.1; f1 = 2.0; % Frequency sweep from f0 to f1 Hz
            amp = 0.05;
            x0_ddot = amp * chirp(t, f0, t_end, f1, 'linear');
            y0_ddot = 0.5 * amp * chirp(t, f0, t_end, f1, 'linear');
            
        case 'step'
            amp = 0.2;
            x0_ddot = amp * (t > 2 & t < 4);
            y0_ddot = 0;
            
        case 'random'
            amp = 0.05;
            x0_ddot = amp * randn(size(t));
            y0_ddot = amp * randn(size(t));

        case 'point_to_point'
            % Point-to-point motion from A to B
            point_A = [0, 0];           % Starting point (m)
            point_B = [2, 0];           % Target point (m)
            max_accel = 0.1;            % Maximum acceleration (m/s^2)
            
            % Calculate motion profile and adjust simulation time
            [t, x0_ddot, y0_ddot, motion_data] = generate_point_to_point_motion(point_A, point_B, max_accel, dt);
            t_end = t(end);
            
            fprintf('Point-to-point motion:\n');
            fprintf('  From: (%.1f, %.1f) to (%.1f, %.1f) m\n', point_A(1), point_A(2), point_B(1), point_B(2));
            fprintf('  Total distance: %.2f m\n', motion_data.total_distance);
            fprintf('  Motion time: %.2f s\n', motion_data.motion_time);
            fprintf('  Settling time: %.2f s\n', motion_data.settling_time);
            fprintf('  Total simulation time adjusted to: %.2f s\n', t_end);
    end
    
    %% Initial Conditions
    % [x_n, x_n_dot, y_n, y_n_dot]
    initial_conditions = [0.01, 0, 0.01, 0];
    
    %% Solve ODE System
    fprintf('Solving nonlinear sloshing equations...\n');
    
    % Create parameter structure
    params.omega_n = omega_n;
    params.zeta_n = zeta_n;
    params.g = g;
    params.R = R;
    params.alpha_n = alpha_n;
    params.t = t;
    params.x0_ddot = x0_ddot;
    params.y0_ddot = y0_ddot;
    
    % Solve using ode45
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t_sol, y_sol] = ode45(@(t, y) sloshing_ode(t, y, params), t, initial_conditions, options);
    
    % Extract solutions
    x_n = y_sol(:, 1);
    x_n_dot = y_sol(:, 2);
    y_n = y_sol(:, 3);
    y_n_dot = y_sol(:, 4);
    
    %% Visualization
    create_visualization(t_sol, x_n, y_n,x0_ddot, x_n_dot, y_n_dot, params, excitation_type, motion_data);
    
    %% Animation
    create_animation(t_sol, x_n, y_n, params);
    
    fprintf('Simulation completed!\n');
end

function dydt = sloshing_ode(t, y, params)
    % ODE system for nonlinear sloshing dynamics
    
    % Extract state variables
    x_n = y(1);      % Position in x
    x_n_dot = y(2);  % Velocity in x
    y_n = y(3);      % Position in y
    y_n_dot = y(4);  % Velocity in y
    
    % Extract parameters
    omega_n = params.omega_n;
    zeta_n = params.zeta_n;
    g = params.g;
    R = params.R;
    alpha_n = params.alpha_n;
    
    % Interpolate excitation at current time
    x0_ddot = interp1(params.t, params.x0_ddot, t, 'linear', 0);
    y0_ddot = interp1(params.t, params.y0_ddot, t, 'linear', 0);
    
    % Common terms
    omega_g_ratio = omega_n^4 / g^2;
    nonlinear_term = alpha_n / R^2 * (x_n^2 + y_n^2);
    
    % Mass matrix coefficients
    M11 = 1 + omega_g_ratio;
    M12 = omega_g_ratio * x_n * y_n;
    M21 = omega_g_ratio * x_n * y_n;
    M22 = 1 + omega_g_ratio * y_n^2;
    
    % Damping and nonlinear terms for x equation
    C1 = 2 * omega_n * zeta_n * (1 + omega_g_ratio * x_n^2) + omega_g_ratio * x_n * x_n_dot;
    C2 = 2 * omega_n * zeta_n * omega_g_ratio * x_n * y_n + omega_g_ratio * x_n * y_n_dot;
    
    % Damping and nonlinear terms for y equation
    C3 = 2 * omega_n * zeta_n * omega_g_ratio * x_n * y_n + omega_g_ratio * y_n * x_n_dot;
    C4 = 2 * omega_n * zeta_n * (1 + omega_g_ratio * y_n^2) + omega_g_ratio * y_n * y_n_dot;
    
    % Spring terms
    K1 = omega_n^2 * x_n * (1 + nonlinear_term);
    K2 = omega_n^2 * y_n * (1 + nonlinear_term);
    
    % Right-hand side
    F1 = -x0_ddot - C1 * x_n_dot - C2 * y_n_dot - K1;
    F2 = -y0_ddot - C3 * x_n_dot - C4 * y_n_dot - K2;
    
    % Solve for accelerations using mass matrix
    det_M = M11 * M22 - M12 * M21;
    
    if abs(det_M) < 1e-10
        warning('Mass matrix is nearly singular');
        x_n_ddot = 0;
        y_n_ddot = 0;
    else
        x_n_ddot = (M22 * F1 - M12 * F2) / det_M;
        y_n_ddot = (M11 * F2 - M21 * F1) / det_M;
    end
    
    % Return derivatives
    dydt = [x_n_dot; x_n_ddot; y_n_dot; y_n_ddot];
end

function create_visualization(t, x_n, y_n,x0_ddot, x_n_dot, y_n_dot, params, excitation_type,motion_data)
    % Create comprehensive visualization
    
    figure('Position', [100, 100, 1200, 800]);
    
    % Time series plots
    subplot(2, 3, 1);
    plot(t, x_n, 'b-', t, motion_data.x0_pos, 'g-');
    xlabel('Time (s)');
    ylabel('x_n (m)');
    title('Sloshing Displacement - X Direction');
    grid on;
    
    subplot(2, 3, 2);
    plot(t, y_n, 'r-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('y_n (m)');
    title('Sloshing Displacement - Y Direction');
    grid on;
    
    % Phase portraits
    subplot(2, 3, 3);
    plot(x_n, x_n_dot, 'b-', 'LineWidth', 1);
    xlabel('x_n (m)');
    ylabel('dx_n/dt (m/s)');
    title('Phase Portrait - X Direction');
    grid on;
    
    subplot(2, 3, 4);
    plot(y_n, y_n_dot, 'r-', 'LineWidth', 1);
    xlabel('y_n (m)');
    ylabel('dy_n/dt (m/s)');
    title('Phase Portrait - Y Direction');
    grid on;
    
    % 2D trajectory
    subplot(2, 3, 5);
    plot(x_n, y_n, 'g-', 'LineWidth', 1.5);
    hold on;
    plot(x_n(1), y_n(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    plot(x_n(end), y_n(end), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    xlabel('x_n (m)');
    ylabel('y_n (m)');
    title('2D Sloshing Trajectory');
    legend('Trajectory', 'Start', 'End', 'Location', 'best');
    grid on;
    axis equal;
    
    % Energy analysis
    subplot(2, 3, 6);
    kinetic_energy = 0.5 * (x_n_dot.^2 + y_n_dot.^2);
    potential_energy = 0.5 * params.omega_n^2 * (x_n.^2 + y_n.^2);
    total_energy = kinetic_energy + potential_energy;
    
    plot(t, kinetic_energy, 'b-', t, potential_energy, 'r-', t, total_energy, 'k--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Energy');
    title('Energy Analysis');
    legend('Kinetic', 'Potential', 'Total', 'Location', 'best');
    grid on;
    
    sgtitle(sprintf('Nonlinear Sloshing Dynamics - %s Excitation', excitation_type));
end

function [t, x0_ddot, y0_ddot, motion_data] = generate_point_to_point_motion(point_A, point_B, max_accel, dt)
    % Generate smooth point-to-point motion profile with trapezoidal velocity
    % and S-curve acceleration to minimize jerk
    
    % Calculate displacement vector
    displacement = point_B - point_A;
    total_distance = norm(displacement);
    direction = displacement / total_distance;
    
    if total_distance == 0
        error('Points A and B are identical');
    end
    
    % Design motion profile parameters
    % Use trapezoidal velocity profile with S-curve acceleration
    jerk_limit = 2.0;  % Jerk limit (m/s^3) for smooth acceleration
    
    % Calculate acceleration ramp time based on jerk limit
    t_ramp = max_accel / jerk_limit;
    
    % Calculate time to reach maximum velocity
    % Assume we reach max_accel, then maintain it, then decelerate
    v_max = sqrt(max_accel * total_distance / 2);  % Triangular profile max velocity
    
    % Check if we can reach max acceleration
    if v_max > max_accel * sqrt(total_distance / max_accel)
        % Can't reach full acceleration - use triangular velocity profile
        v_max = sqrt(max_accel * total_distance);
        t_accel = v_max / max_accel;
        t_const = 0;
        t_decel = t_accel;
    else
        % Can reach full acceleration - use trapezoidal velocity profile
        t_accel = v_max / max_accel;
        t_const = (total_distance - v_max^2 / max_accel) / v_max;
        t_decel = t_accel;
    end
    
    % Total motion time
    motion_time = 2 * t_ramp + t_accel + t_const + t_decel + 2 * t_ramp;  % Include S-curve ramps
    
    % Add settling time for sloshing to settle
    settling_time = max(20, 5 / (2*pi*0.5));  % At least 5 periods of natural frequency
    total_time = motion_time + settling_time;
    
    % Create time vector
    t = 0:dt:total_time;
    
    % Initialize motion vectors
    s = zeros(size(t));      % Position along path
    s_dot = zeros(size(t));  % Velocity along path  
    s_ddot = zeros(size(t)); % Acceleration along path
    
    % Generate S-curve motion profile
    for i = 1:length(t)
        time = t(i);
        
        if time <= t_ramp
            % First jerk phase (acceleration ramp up)
            j = jerk_limit;
            s_ddot(i) = 0.5 * j * time^2;
            s_dot(i) = (1/6) * j * time^3;
            s(i) = (1/24) * j * time^4;
            
        elseif time <= t_ramp + t_accel
            % Constant acceleration phase
            t_phase = time - t_ramp;
            a0 = 0.5 * jerk_limit * t_ramp^2;
            v0 = (1/6) * jerk_limit * t_ramp^3;
            s0 = (1/24) * jerk_limit * t_ramp^4;
            
            s_ddot(i) = max_accel;
            s_dot(i) = v0 + max_accel * t_phase;
            s(i) = s0 + v0 * t_phase + 0.5 * max_accel * t_phase^2;
            
        elseif time <= 2*t_ramp + t_accel
            % Second jerk phase (acceleration ramp down)
            t_phase = time - t_ramp - t_accel;
            t_prev = t_ramp + t_accel;
            
            % Get values at start of this phase
            a0 = max_accel;
            v0 = (1/6) * jerk_limit * t_ramp^3 + max_accel * t_accel;
            s0 = (1/24) * jerk_limit * t_ramp^4 + (1/6) * jerk_limit * t_ramp^3 * t_accel + 0.5 * max_accel * t_accel^2;
            
            j = -jerk_limit;
            s_ddot(i) = a0 + j * t_phase;
            s_dot(i) = v0 + a0 * t_phase + 0.5 * j * t_phase^2;
            s(i) = s0 + v0 * t_phase + 0.5 * a0 * t_phase^2 + (1/6) * j * t_phase^3;
            
        elseif time <= 2*t_ramp + t_accel + t_const
            % Constant velocity phase
            t_phase = time - 2*t_ramp - t_accel;
            
            % Get values at start of constant velocity phase
            v_const = (1/6) * jerk_limit * t_ramp^3 + max_accel * t_accel;
            s0 = (1/24) * jerk_limit * t_ramp^4 + (1/6) * jerk_limit * t_ramp^3 * t_accel + 0.5 * max_accel * t_accel^2 + ...
                 v_const * t_ramp - 0.5 * max_accel * t_ramp + (1/6) * jerk_limit * t_ramp^3;
            
            s_ddot(i) = 0;
            s_dot(i) = v_const;
            s(i) = s0 + v_const * t_phase;
            
        elseif time <= motion_time
            % Deceleration phase (mirror of acceleration)
            t_remaining = motion_time - time;
            
            if t_remaining <= t_ramp
                % Final jerk phase
                j = -jerk_limit;
                s_ddot(i) = 0.5 * j * t_remaining^2;
                s_dot(i) = -(1/6) * j * t_remaining^3;
                
            elseif t_remaining <= t_ramp + t_decel  
                % Constant deceleration
                t_phase = t_remaining - t_ramp;
                s_ddot(i) = -max_accel;
                s_dot(i) = (1/6) * jerk_limit * t_ramp^3 + max_accel * t_phase;
                
            else
                % Initial deceleration jerk
                t_phase = t_remaining - t_ramp - t_decel;
                s_ddot(i) = -0.5 * jerk_limit * t_phase^2;
                s_dot(i) = (1/6) * jerk_limit * t_phase^3 + (1/6) * jerk_limit * t_ramp^3 + max_accel * t_decel;
            end
            
            % Calculate position by integrating from the end
            s(i) = total_distance - (s(end) - s(i));
            
        else
            % Settling phase - container at rest
            s_ddot(i) = 0;
            s_dot(i) = 0;
            s(i) = total_distance;
        end
    end
    
    % Ensure we end exactly at the target distance
    s = s * total_distance / max(s);
    
    % Convert to x and y coordinates
    x0_pos = point_A(1) + s * direction(1);
    y0_pos = point_A(2) + s * direction(2);
    x0_vel = s_dot * direction(1);
    y0_vel = s_dot * direction(2);
    x0_ddot = s_ddot * direction(1);
    y0_ddot = s_ddot * direction(2);
    
    % Store motion data for visualization
    motion_data.x0_pos = x0_pos;
    motion_data.y0_pos = y0_pos;
    motion_data.x0_vel = x0_vel;
    motion_data.y0_vel = y0_vel;
    motion_data.total_distance = total_distance;
    motion_data.motion_time = motion_time;
    motion_data.settling_time = settling_time;
    motion_data.max_velocity = max(abs(s_dot));
    motion_data.max_acceleration = max(abs(s_ddot));
end

function create_animation(t, x_n, y_n, params)
    % Create animation of liquid sloshing
    
    figure('Position', [200, 200, 800, 600]);
    
    % Tank parameters for visualization
    tank_radius = params.R;
    
    % Create tank outline
    theta = linspace(0, 2*pi, 100);
    tank_x = tank_radius * cos(theta);
    tank_y = tank_radius * sin(theta);
    
    % Animation parameters
    skip_frames = 50; % Show every 50th frame for smoother animation
    
    for i = 1:skip_frames:length(t)
        clf;
        
        % Plot tank
        plot(tank_x, tank_y, 'k-', 'LineWidth', 3);
        hold on;
        
        % Plot liquid surface (simplified representation)
        surface_amplitude = max(abs([x_n(i), y_n(i)]));
        surface_x = linspace(-tank_radius, tank_radius, 50);
        surface_y = x_n(i) * sin(pi * surface_x / tank_radius) .* exp(-abs(surface_x)/tank_radius);
        
        % Fill liquid area
        liquid_bottom = -tank_radius * 0.8;
        fill([surface_x, fliplr(surface_x)], ...
             [surface_y, liquid_bottom * ones(size(surface_x))], ...
             [0.3, 0.6, 1.0], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
        
        % Plot liquid surface
        plot(surface_x, surface_y, 'b-', 'LineWidth', 2);
        
        % Plot sloshing point
        plot(x_n(i), y_n(i), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        
        % Plot trajectory up to current time
        if i > 1
            plot(x_n(1:i), y_n(1:i), 'r--', 'LineWidth', 1, 'Color', [1, 0, 0, 0.5]);
        end
        
        % Formatting
        axis equal;
        xlim([-tank_radius*1.2, tank_radius*1.2]);
        ylim([-tank_radius*1.2, tank_radius*1.2]);
        xlabel('X Position (m)');
        ylabel('Y Position (m)');
        title(sprintf('Liquid Sloshing Animation - Time: %.2f s', t(i)));
        grid on;
        
        % Add text information
        text(-tank_radius, tank_radius*1.1, sprintf('x_n = %.3f m', x_n(i)), 'FontSize', 10);
        text(-tank_radius, tank_radius*1.0, sprintf('y_n = %.3f m', y_n(i)), 'FontSize', 10);
        
        drawnow;
        pause(0.05);
    end
end