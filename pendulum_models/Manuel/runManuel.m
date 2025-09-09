container_params.V = 0.5 / 1000;  % 0.5L volume
container_params.rho_l = 1000;     % Water density
container_params.l = 0.021;        % Pendulum length
container_params.d = 1.51;         % Damping constant

%% traj gen
% Points and times
p0 = [0; 0; 0];
p1 = [0.5; 0.8; 0];
p2 = [1; 1; 0];
t0 = 0;  t2 = 2;

% Choose endpoint velocities (try nonzero to shape curvature)
v0 = [0; 0; 0];       % start at rest
v2 = [0; 0; 0];       % end at rest

[u_pp, traj] = accelSpline3PointsClamped(p0,p1,p2,t0,t2,v0,v2);

% Plot per-axis acceleration
tt = linspace(t0,t2,300);
U  = u_pp(tt);
figure;
plot(tt, U(1,:), 'LineWidth', 2); hold on;
plot(tt, U(2,:), 'LineWidth', 2);
plot(tt, U(3,:), 'LineWidth', 2); grid on;
xlabel('Time [s]'); ylabel('Acceleration [units/s^2]');
legend('a_x','a_y','a_z'); title('Per-axis acceleration (clamped spline)');

%% step functions
u_step = @(t) [0 ; 1.0 * (t >= 0.5 & t <= 1.0); 0];  % 1 m/s² step from 0.5s to 1.0s
u_step2 = @(t) [1.0 * (t >= 1.5 & t <= 2.0) - 1.0  * (t >= 2 & t <= 2.5); ...
               1.0 * (t >= 1.5 & t <= 2.0) - 1.0  * (t >= 2 & t <= 2.5); ...
               0];
u = u_step2;

% Create model and run demo
model = SloshingModel3D(container_params);

% model.runImpulsDemo(); % theoretical impuls response 1 dim
model.runDemo(u_pp);