function u = trajectory_input(t)
    % Define waypoints: [time, x_position, y_position]
    waypoints = [
        10,   1,  0.0;    % Move to (0.5, 0) by t=5s
        15,  5,  0;    % Move to (0.5, 0.3) by t=10s
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