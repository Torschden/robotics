function [theta, theta_dot, x_c, x_pendulum, y_pendulum] = simulation(t,p)
    %[x_c, x_c_dot, x_c_ddot] = cart_motion(t, params);
    %[theta, theta_dot] =  simulate_pendulum_angular(t, x_c_ddot, p);
    
    
    [x_c, x_c_dot, theta, theta_dot] = simulate_closed_loop(t,p);

    [x_pendulum, y_pendulum] = calculate_positions(theta, x_c, p);
end

function x_c_ddot = pendulum_feedback_controller(theta, theta_dot, p)
    Kp = 2;
    Kd = - 0.4;
    x_c_ddot = -Kp * theta - Kd * theta_dot;
end

function [x_c, x_c_dot, x_c_ddot, theta, theta_dot] = simulate_closed_loop(t,p)
    n_steps = length(t);
    x_c = zeros(n_steps, 1);
    x_c_dot = zeros(n_steps, 1);
    x_c_ddot = zeros(n_steps, 1);
    theta = zeros(n_steps, 1);
    theta_dot = zeros(n_steps, 1);

    theta(1) = p.theta_0;
    theta_dot(1) = p.theta_dot_0;
    x_c(1) = 0;
    x_c_dot(1) = 0;

    for i = 1:n_steps-1
        % Controller
        u = pendulum_feedback_controller(theta(i), theta_dot(i), p);
        x_c_ddot(i) = u;

        % Cart dynamics (Euler integration)
        x_c_dot(i+1) = x_c_dot(i) + x_c_ddot(i) * p.dt;
        x_c(i+1) = x_c(i) + x_c_dot(i+1) * p.dt;

        % Pendulum dynamics
        theta_ddot = - (p.b / (p.m_p * p.L^2)) * theta_dot(i) ...
                     - (p.g / p.L) * sin(theta(i)) ...
                     - (x_c_ddot(i) / p.L) * cos(theta(i));
        theta_dot(i+1) = theta_dot(i) + theta_ddot * p.dt;
        theta(i+1) = theta(i) + theta_dot(i+1) * p.dt;
    end
    x_c_ddot(end) = pendulum_feedback_controller(theta(end), theta_dot(end), p);
end

function [x_c, x_c_dot, x_c_ddot] = cart_motion(t, p)
    % Define cart motion pattern

    % % oscillating
    % freq = 0.5;         % frequency (Hz)
    % amplitude = 0.3;    % amplitude (m)
    % x_c = amplitude * sin(2*pi*freq*t);
    % x_c_dot = amplitude * 2*pi*freq * cos(2*pi*freq*t);
    % x_c_ddot = -amplitude * (2*pi*freq)^2 * sin(2*pi*freq*t);

    % % constant vel
    % x_c = 0.1 * t;
    % x_c_dot = 0.1 * ones(size(t));
    % x_c_ddot = zeros(size(t));

    % speeding up / slowing down from 0->B
    x_goal = 2;
    t_acc = p.v_max / p.a_max;
    x_acc = 0.5 * p.a_max * t_acc^2;

    if 2 * x_acc >= x_goal
        x_acc = 0.5 * x_goal;
        t_acc = sqrt(x_acc * 2 / p.a_max);
        t_flat = 0;
        t_total = 2 * t_acc;
        v_peak = p.a_max * t_acc; 
    else
        x_flat = x_goal - 2 * x_acc;
        t_flat = x_flat / p.v_max;
        t_total = 2 * t_acc + t_flat;
        v_peak = p.v_max;
    end

    t = t(:); % Ensure column vector

    x_c = zeros(size(t));
    x_c_dot = zeros(size(t));
    x_c_ddot = zeros(size(t));

    for i = 1:length(t)
        ti = t(i);

        if ti < t_acc
            x_c(i) = 0.5 * p.a_max * ti^2;
            x_c_dot(i) = p.a_max * ti;
            x_c_ddot(i) = p.a_max;
        elseif ti < (t_acc + t_flat)
            dt = ti - t_acc;
            x_c(i) = x_acc + p.v_max * dt;
            x_c_dot(i) = p.v_max;
            x_c_ddot(i) = 0;
        elseif ti < t_total
            dt = ti - (t_acc + t_flat);
            x_c(i) = x_acc + p.v_max * t_flat + v_peak * dt - 0.5 * p.a_max * dt^2;
            x_c_dot(i) = p.v_max - p.a_max * dt;
            x_c_ddot(i) = -p.a_max;
        else
            x_c(i) = x_goal;
            x_c_dot(i) = 0;
            x_c_ddot(i) = 0;
        end
    end
end

function [theta, theta_dot] = simulate_pendulum_angular(t, x_c_ddot, p)
    % Simulate pendulum using angular equation of motion (cart drives pendulum)
    n_steps = length(t);
    theta = zeros(n_steps, 1);
    theta_dot = zeros(n_steps, 1);
    theta(1) = p.theta_0;
    theta_dot(1) = p.theta_dot_0;
    dt = p.dt;
    
    for i = 1:n_steps-1
        % Equation: θ'' + (b/mL^2) θ' + (g/L) sinθ = -x_c_ddot/L * cosθ
        theta_ddot = - (p.b/(p.m_p*p.L^2)) * theta_dot(i) ...
                     - (p.g/p.L) * sin(theta(i)) ...
                     - (x_c_ddot(i)/p.L) * cos(theta(i));
        % Euler integration (can be replaced with RK4 for more accuracy)
        theta_dot(i+1) = theta_dot(i) + theta_ddot * dt;
        theta(i+1) = theta(i) + theta_dot(i+1) * dt;
    end
end

function [x_pendulum, y_pendulum] = calculate_positions(theta, x_c, p)
    % Calculate pendulum bob positions
    x_pendulum = x_c(:) + p.L * sin(theta(:));
    y_pendulum = -p.L * cos(theta(:));
end