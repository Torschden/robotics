% Slosh Simulation with Position, Velocity, and Slosh Angle
clear; clc;

% === System Parameters ===
omega_n = 2*pi*1.68;     % Natural frequency (rad/s)
zeta = 0.0083;           % Damping ratio
rm = 0.06;               % Virtual pendulum length (m)
T = 0.001;               % Sampling time (s)
K = 2;                   % Time compression factor (for filter design)

% === Time Vector ===
t_end = 5;               % Simulation time (s)
t = 0:T:t_end;
N = length(t);

% === Input: Reference Acceleration (e.g., square pulse) ===
u = zeros(1, N);
u(t >= 0.5 & t <= 1.5) = 1;   % Unit acceleration pulse for 1 second
u(t > 3.5 & t <= 4.5) = -1;
% === Generate Filter: Critically Damped 3rd Order / Underdamped 2nd Order ===
G_underdamped = tf(1, [1 2*zeta*omega_n omega_n^2]);
omega_c = omega_n;
G_crit_damped = tf(1, [1 3*omega_c 3*omega_c^2 omega_c^3]);

% Input shaping filter: G_shaper = G_target / G_original
G_shaper = G_crit_damped / G_underdamped;
G_shaper_d = c2d(G_shaper, T, 'tustin');

% === Apply Filter ===
acc_filtered = lsim(G_shaper_d, u, t);  % Filtered acceleration

% === Integrate to Get Velocity and Position ===
velocity = cumtrapz(t, acc_filtered);      % First integral
position = cumtrapz(t, velocity);          % Second integral

% === Simulate Slosh Angle Dynamics ===
theta = zeros(1, N);
theta_dot = 0;

for k = 2:N
    acc = acc_filtered(k);
    theta_ddot = -2*zeta*omega_n*theta_dot - omega_n^2*theta(k-1) - acc/rm;
    theta_dot = theta_dot + theta_ddot*T;
    theta(k) = theta(k-1) + theta_dot*T;
end

% === Plot Results ===

figure;
subplot(3,1,1);
plot(t, position, 'LineWidth', 1.5);
ylabel('Position (m)');
title('Filtered Motion Profile & Slosh Angle');
grid on;

subplot(3,1,2);
plot(t, velocity, 'LineWidth', 1.5);
ylabel('Velocity (m/s)');
grid on;

subplot(3,1,3);
plot(t, rad2deg(theta), 'LineWidth', 1.5);
ylabel('\theta (deg)');
xlabel('Time (s)');
grid on;
