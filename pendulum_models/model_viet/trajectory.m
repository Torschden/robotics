function u = trajectory(t)
    % Generate smooth point-to-point motion profile with trapezoidal velocity
    % and S-curve acceleration to minimize jerk
    
    

    point_A = [0, 0];
    point_B = [1.0, 0];
    max_accel = 0.15;

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