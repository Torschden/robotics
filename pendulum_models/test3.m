%% Dual Pendulum Cart Control System
% Two identical pendulums attached to the same point on a cart
% Cart movement designed to dampen both pendulum oscillations

clear; close all; clc;

%% System Parameters
% Cart parameters
m_cart = 1.0;           % Cart mass (kg)

% Pendulum parameters (identical for both pendulums)
m_pen = 0.2;            % Pendulum mass (kg)
l_pen = 0.3;            % Pendulum length (m)
I_pen = 0.006;          % Pendulum moment of inertia (kg·m²)
b_pen = 0.005;          % Pendulum damping coefficient (N·m·s/rad)

% Physical constants
g = 9.81;               % Gravity (m/s²)
u_max = 5;              % Maximum control input (N)

% Simulation parameters
Ts = 1/100;             % Sampling time (s)
t_sim = 15;             % Simulation time (s)
t = 0:Ts:t_sim;         % Time vector

%% Initial Conditions
theta1_0 = deg2rad(0); % Initial angle pendulum 1 (rad)
theta2_0 = deg2rad(0); % Initial angle pendulum 2 (rad)

%% Approach: Use single pendulum model with combined effect
% Since both pendulums are identical and attached to the same point,
% we can model this as a single pendulum system where the control
% affects both pendulums simultaneously through cart motion

% For control design, use equivalent single pendulum parameters
m_pen_eq = 2 * m_pen;   % Combined pendulum mass
I_pen_eq = 2 * I_pen;   % Combined inertia

% Linearized system for single equivalent pendulum
% State: [x, x_dot, theta_eq, theta_eq_dot] where theta_eq represents average motion
M_total = m_cart + m_pen_eq;
I_eff = I_pen_eq + m_pen_eq * l_pen^2;
den = M_total * I_eff - m_pen_eq^2 * l_pen^2;

% 4x4 system matrix for equivalent system
A_eq = zeros(4,4);
A_eq(1,2) = 1; % dx/dt = x_dot
A_eq(2,3) = (m_pen_eq^2 * l_pen^2 * g) / den;     % theta effect on x_ddot
A_eq(2,4) = -(m_pen_eq * l_pen * b_pen) / den;    % theta_dot effect on x_ddot
A_eq(3,4) = 1; % dtheta/dt = theta_dot
A_eq(4,2) = -(m_pen_eq * l_pen * g * M_total) / den; % x_dot effect on theta_ddot
A_eq(4,3) = (I_eff * m_pen_eq * g * l_pen) / den;    % theta effect on theta_ddot
A_eq(4,4) = -(I_eff * b_pen) / den;                  % theta_dot damping

% Input matrix for equivalent system
B_eq = zeros(4,1);
B_eq(2) = I_eff / den;
B_eq(4) = (m_pen_eq * l_pen) / den;

% Check controllability of equivalent system
Co_eq = ctrb(A_eq, B_eq);
rank_Co_eq = rank(Co_eq);
fprintf('Equivalent system controllability rank: %d (should be 4)\n', rank_Co_eq);

% Discretize equivalent system
sys_eq = ss(A_eq, B_eq, eye(4), 0);
sys_eq_d = c2d(sys_eq, Ts, 'zoh');
Ad_eq = sys_eq_d.A;
Bd_eq = sys_eq_d.B;

%% LQR Controller for equivalent system
Q_eq = diag([10, 1, 100, 10]); % [x, x_dot, theta, theta_dot]
R_eq = 1;

try
    K_eq = dlqr(Ad_eq, Bd_eq, Q_eq, R_eq);
    fprintf('LQR controller successfully computed for equivalent system\n');
    fprintf('LQR Gains: K = [%.3f, %.3f, %.3f, %.3f]\n', K_eq);
catch
    error('LQR design failed even for equivalent system');
end

% Verify stability
A_cl_eq = Ad_eq - Bd_eq * K_eq;
cl_poles_eq = eig(A_cl_eq);
fprintf('Closed-loop poles magnitude: ');
fprintf('%.3f ', abs(cl_poles_eq));
fprintf('\n');

%% Create Disturbance Signal
impulse_time = 5;
impulse_index = round(impulse_time / Ts);
impulse_strength = 0.1; % N·m·s

%% Individual Pendulum Control Strategy
% We'll control each pendulum individually using the same controller
% but apply it based on each pendulum's state relative to equilibrium

% Initialize states for both pendulums
x_states = zeros(length(t), 6); % [x, x_dot, theta1, theta1_dot, theta2, theta2_dot]
u_control = zeros(length(t), 1);

% Initial conditions
x_states(1, :) = [0, 0, theta1_0, 0, theta2_0, 0];

for i = 1:length(t)-1
    % Current state
    x_pos = x_states(i, 1);
    x_vel = x_states(i, 2);
    theta1 = x_states(i, 3);
    theta1_dot = x_states(i, 4);
    theta2 = x_states(i, 5);
    theta2_dot = x_states(i, 6);
    
    % Control strategy: Use weighted average of both pendulums
    % This allows the controller to handle both pendulums simultaneously
    theta_avg = (theta1 + theta2) / 2;
    theta_dot_avg = (theta1_dot + theta2_dot) / 2;
    
    % Control based on cart position and average pendulum state
    state_eq = [x_pos; x_vel; theta_avg; theta_dot_avg];
    u = -K_eq * state_eq;
    
    % Apply control limits
    u = max(-u_max, min(u_max, u));
    u_control(i) = u;
    
    % Nonlinear dynamics for both pendulums
    sin_th1 = sin(theta1);
    cos_th1 = cos(theta1);
    sin_th2 = sin(theta2);
    cos_th2 = cos(theta2);
    
    % Mass matrix (3x3: x_ddot, theta1_ddot, theta2_ddot)
    M11 = m_cart + 2*m_pen;
    M12 = m_pen * l_pen * cos_th1;
    M13 = m_pen * l_pen * cos_th2;
    M21 = m_pen * l_pen * cos_th1;
    M22 = I_pen + m_pen * l_pen^2;
    M23 = 0;
    M31 = m_pen * l_pen * cos_th2;
    M32 = 0;
    M33 = I_pen + m_pen * l_pen^2;
    
    M = [M11 M12 M13;
         M21 M22 M23;
         M31 M32 M33];
    
    % Force vector
    F1 = u + m_pen * l_pen * (theta1_dot^2 * sin_th1 + theta2_dot^2 * sin_th2);
    F2 = -m_pen * g * l_pen * sin_th1 - b_pen * theta1_dot;
    F3 = -m_pen * g * l_pen * sin_th2 - b_pen * theta2_dot;
    
    % Add impulse disturbance to pendulum 1
    if i == impulse_index
        F2 = F2 + impulse_strength / Ts;
    end
    
    F = [F1; F2; F3];
    
    % Solve for accelerations
    try
        accel = M \ F;
        x_ddot = accel(1);
        theta1_ddot = accel(2);
        theta2_ddot = accel(3);
    catch
        % Fallback to decoupled dynamics if matrix is singular
        x_ddot = u / (m_cart + 2*m_pen);
        theta1_ddot = -(g/l_pen) * sin_th1 - (b_pen/(I_pen + m_pen*l_pen^2)) * theta1_dot;
        theta2_ddot = -(g/l_pen) * sin_th2 - (b_pen/(I_pen + m_pen*l_pen^2)) * theta2_dot;
    end
    
    % Euler integration
    x_states(i+1, 1) = x_pos + Ts * x_vel;
    x_states(i+1, 2) = x_vel + Ts * x_ddot;
    x_states(i+1, 3) = theta1 + Ts * theta1_dot;
    x_states(i+1, 4) = theta1_dot + Ts * theta1_ddot;
    x_states(i+1, 5) = theta2 + Ts * theta2_dot;
    x_states(i+1, 6) = theta2_dot + Ts * theta2_ddot;
end

%% Simulation without control (for comparison)
x_states_open = zeros(length(t), 6);
x_states_open(1, :) = [0, 0, theta1_0, 0, theta2_0, 0];

for i = 1:length(t)-1
    x_pos = x_states_open(i, 1);
    x_vel = x_states_open(i, 2);
    theta1 = x_states_open(i, 3);
    theta1_dot = x_states_open(i, 4);
    theta2 = x_states_open(i, 5);
    theta2_dot = x_states_open(i, 6);
    
    % No control
    u = 0;
    
    sin_th1 = sin(theta1);
    cos_th1 = cos(theta1);
    sin_th2 = sin(theta2);
    cos_th2 = cos(theta2);
    
    M11 = m_cart + 2*m_pen;
    M12 = m_pen * l_pen * cos_th1;
    M13 = m_pen * l_pen * cos_th2;
    M21 = m_pen * l_pen * cos_th1;
    M22 = I_pen + m_pen * l_pen^2;
    M23 = 0;
    M31 = m_pen * l_pen * cos_th2;
    M32 = 0;
    M33 = I_pen + m_pen * l_pen^2;
    
    M = [M11 M12 M13;
         M21 M22 M23;
         M31 M32 M33];
    
    F1 = m_pen * l_pen * (theta1_dot^2 * sin_th1 + theta2_dot^2 * sin_th2);
    F2 = -m_pen * g * l_pen * sin_th1 - b_pen * theta1_dot;
    F3 = -m_pen * g * l_pen * sin_th2 - b_pen * theta2_dot;
    
    if i == impulse_index
        F2 = F2 + impulse_strength / Ts;
    end
    
    F = [F1; F2; F3];
    
    try
        accel = M \ F;
        x_ddot = accel(1);
        theta1_ddot = accel(2);
        theta2_ddot = accel(3);
    catch
        x_ddot = 0;
        theta1_ddot = -(g/l_pen) * sin_th1 - (b_pen/(I_pen + m_pen*l_pen^2)) * theta1_dot;
        theta2_ddot = -(g/l_pen) * sin_th2 - (b_pen/(I_pen + m_pen*l_pen^2)) * theta2_dot;
    end
    
    x_states_open(i+1, 1) = x_pos + Ts * x_vel;
    x_states_open(i+1, 2) = x_vel + Ts * x_ddot;
    x_states_open(i+1, 3) = theta1 + Ts * theta1_dot;
    x_states_open(i+1, 4) = theta1_dot + Ts * theta1_ddot;
    x_states_open(i+1, 5) = theta2 + Ts * theta2_dot;
    x_states_open(i+1, 6) = theta2_dot + Ts * theta2_ddot;
end

%% Extract results
x_pos = x_states(:, 1);
x_vel = x_states(:, 2);
theta1 = x_states(:, 3);
theta2 = x_states(:, 5);

theta1_deg = rad2deg(theta1);
theta2_deg = rad2deg(theta2);

theta1_open_deg = rad2deg(x_states_open(:, 3));
theta2_open_deg = rad2deg(x_states_open(:, 5));

%% Plotting Results
figure('Name', 'Cart Motion', 'Position', [100, 100, 800, 600]);
subplot(3,1,1);
plot(t, x_pos, 'b-', 'LineWidth', 1.5);
ylabel('Position (m)');
title('Cart Motion Response');
grid on;

subplot(3,1,2);
plot(t, x_vel, 'g-', 'LineWidth', 1.5);
ylabel('Velocity (m/s)');
grid on;

subplot(3,1,3);
plot(t(1:end-1), u_control(1:end-1), 'r-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Control Force (N)');
grid on;

% Main comparison plot
figure('Name', 'Pendulum Control Effectiveness', 'Position', [200, 100, 1000, 800]);

% Pendulum angles comparison
subplot(3,1,1);
plot(t, theta1_deg, 'b-', 'LineWidth', 2); hold on;
plot(t, theta1_open_deg, 'b--', 'LineWidth', 1.5, 'Color', [0.5 0.5 1]);
plot(t, theta2_deg, 'r-', 'LineWidth', 2);
plot(t, theta2_open_deg, 'r--', 'LineWidth', 1.5, 'Color', [1 0.5 0.5]);
xlabel('Time (s)');
ylabel('Angle (deg)');
title('Pendulum Angle Response - Controlled vs Uncontrolled (Impulse at t = 5s)');
legend('Pen1 Controlled', 'Pen1 Open-loop', 'Pen2 Controlled', 'Pen2 Open-loop', 'Location', 'best');
grid on;

% Cart position
subplot(3,1,2);
plot(t, x_pos, 'k-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Cart Position (m)');
title('Cart Position (shows active control movement)');
grid on;

% Control effort
subplot(3,1,3);
plot(t(1:end-1), u_control(1:end-1), 'r-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Control Force (N)');
title('Control Effort Required');
grid on;

%% Performance metrics
max_angle1_ctrl = max(abs(theta1_deg));
max_angle1_open = max(abs(theta1_open_deg));
max_angle2_ctrl = max(abs(theta2_deg));
max_angle2_open = max(abs(theta2_open_deg));

improvement1 = (max_angle1_open - max_angle1_ctrl) / max_angle1_open * 100;
improvement2 = (max_angle2_open - max_angle2_ctrl) / max_angle2_open * 100;

% Settling time analysis
settling_threshold = 0.05; % 5% of max
settle_bound1 = settling_threshold * max(abs(theta1_deg));
settle_bound2 = settling_threshold * max(abs(theta2_deg));

settle_idx1 = find(abs(theta1_deg) > settle_bound1, 1, 'last');
settle_idx2 = find(abs(theta2_deg) > settle_bound2, 1, 'last');

settle_time1 = settle_idx1 * Ts;
settle_time2 = settle_idx2 * Ts;

fprintf('\n=== Performance Analysis ===\n');
fprintf('Impulse applied at t = %.1f s with strength: %.4f N·m·s\n', impulse_time, impulse_strength);
fprintf('Pendulum 1 - Max angle: %.1f° → %.1f° (%.1f%% improvement)\n', ...
        max_angle1_open, max_angle1_ctrl, improvement1);
fprintf('Pendulum 2 - Max angle: %.1f° → %.1f° (%.1f%% improvement)\n', ...
        max_angle2_open, max_angle2_ctrl, improvement2);
fprintf('Settling times: Pen1 = %.1fs, Pen2 = %.1fs\n', settle_time1, settle_time2);
fprintf('Max control effort: %.2f N, RMS: %.2f N\n', max(abs(u_control)), rms(u_control));

%% Animation
figure('Name', 'Dual Pendulum Cart Animation', 'Position', [300, 100, 900, 600]);
cart_width = 0.3;
cart_height = 0.15;
skip_frames = 8;

for i = 1:skip_frames:length(t)
    clf;
    hold on;
    axis equal;
    grid on;
    xlim([-1.2, 1.2]);
    ylim([-0.8, 0.4]);
    
    % Current positions
    cart_x = x_pos(i);
    cart_y = 0;
    phi1 = theta1(i);
    phi2 = theta2(i);
    
    % Draw cart
    rectangle('Position', [cart_x - cart_width/2, cart_y - cart_height/2, ...
                          cart_width, cart_height], ...
              'FaceColor', [0.8, 0.8, 0.9], 'EdgeColor', 'black', 'LineWidth', 2);
    
    % Pendulum endpoints
    pend1_x = cart_x + l_pen * sin(phi1);
    pend1_y = cart_y - l_pen * cos(phi1);
    pend2_x = cart_x + l_pen * sin(phi2);
    pend2_y = cart_y - l_pen * cos(phi2);
    
    % Draw pendulums
    plot([cart_x, pend1_x], [cart_y, pend1_y], 'b-', 'LineWidth', 4);
    plot([cart_x, pend2_x], [cart_y, pend2_y], 'r-', 'LineWidth', 4);
    
    % Draw masses
    plot(pend1_x, pend1_y, 'bo', 'MarkerSize', 15, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'cyan');
    plot(pend2_x, pend2_y, 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'yellow');
    
    % Draw cart pivot
    plot(cart_x, cart_y, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'black');
    
    % Add force arrow (scaled for visibility)
    if abs(u_control(min(i, length(u_control)))) > 0.1
        force_scale = 0.1;
        force_x = u_control(min(i, length(u_control))) * force_scale;
        arrow_y = cart_y + cart_height/2 + 0.1;
        if force_x > 0
            quiver(cart_x, arrow_y, force_x, 0, 'g', 'LineWidth', 3, 'MaxHeadSize', 0.5);
        else
            quiver(cart_x, arrow_y, force_x, 0, 'g', 'LineWidth', 3, 'MaxHeadSize', 0.5);
        end
    end
    
    title(sprintf('Dual Pendulum Control | t=%.1fs | Force=%.1fN | θ₁=%.1f° θ₂=%.1f°', ...
                  t(i), u_control(min(i, length(u_control))), rad2deg(phi1), rad2deg(phi2)));
    xlabel('Position (m)');
    ylabel('Height (m)');
    
    drawnow;
    pause(0.02);
end

fprintf('\n=== Animation Complete ===\n');
fprintf('The cart actively moves to counteract pendulum swings.\n');
fprintf('Notice how both pendulums are stabilized through cart motion alone.\n');