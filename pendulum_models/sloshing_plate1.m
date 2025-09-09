clear
close all

% Parameters
g = 9.81; % gravity [m/s^2]
mass1 = 0.2; % vessel 1 liquid mass [kg] ~500ml water
mass2 = 0.3; % vessel 2 liquid mass [kg] ~300ml water
h_glass = 0.06; % glass height [m]
d_glass = 0.08; % glass diameter [m]
plate_mass = 0.50; % plate mass [kg]
plate_size = 0; % plate size [m]
plate_damp = 0; % plate damping [N*s/m]

% IMPROVED: More realistic pendulum lengths (liquid sloshing usually ~1/3 to 1/2 of height)
l1 = h_glass; % pendulum length for vessel 1 [m] - increased from h/4
l2 = h_glass/2; % pendulum length for vessel 2 [m] - increased from h/4

I1 = mass1 * l1^2; % moment of inertia vessel 1
I2 = mass2 * l2^2; % moment of inertia vessel 2

% IMPROVED: More realistic sloshing parameters
% Natural frequency should be around √(g/l) for small oscillations
k1 = mass1 * g / l1; % CHANGED: removed /2 factor - this was making it too stiff
k2 = mass2 * g / l2; % CHANGED: removed /2 factor

% IMPROVED: Much lower damping for realistic liquid behavior
c1 = 0.00255; % REDUCED: was 0.05, now 0.005 (10x less damping)
c2 = 0; % REDUCED: was 0.01, now 0.003

% Initial conditions
theta10 = deg2rad(10); % initial deflection [rad]
dtheta10 = 0; % initial angular velocity [rad/s]
theta20 = 0; % vessel 2 starts at rest
dtheta20 = 0;
x_plate0 = 0; % plate position [m]
dx_plate0 = 0; % plate velocity [m/s]

% IMPROVED: Reduced coupling coefficients for more realistic interaction
alpha_impact = 0; % REDUCED: was 1.0, now 0.1 (less aggressive wall impact)
beta_plate = 0.5; % REDUCED: was 1.0, now 0.5

Tmax = 6; % simulation time [s]

% Pack parameters into a structure
params = struct('g', g, 'mass1', mass1, 'mass2', mass2, 'l1', l1, 'l2', l2, ...
    'I1', I1, 'I2', I2, 'k1', k1, 'k2', k2, 'c1', c1, 'c2', c2, ...
    'plate_mass', plate_mass, 'plate_damp', plate_damp, ...
    'alpha_impact', alpha_impact, 'beta_plate', beta_plate);

% Initial state vector
x0 = [theta10; dtheta10; theta20; dtheta20; x_plate0; dx_plate0];

% Display expected natural frequencies
f1_expected = sqrt(g/l1)/(2*pi);
f2_expected = sqrt(g/l2)/(2*pi);
fprintf('Expected natural frequencies:\n');
fprintf('Vessel 1: %.2f Hz (period = %.2f s)\n', f1_expected, 1/f1_expected);
fprintf('Vessel 2: %.2f Hz (period = %.2f s)\n', f2_expected, 1/f2_expected);

% Integrate system
[t, sol] = ode45(@(t,x) coupled_vessel_slosh_improved(t, x, params), [0 Tmax], x0);

% Visualization
figure('Position', [100, 100, 800, 600]);
subplot(2,1,1);
plot(t, rad2deg(sol(:,1)), 'b-', 'LineWidth', 2);
hold on;
plot(t, rad2deg(sol(:,3)), 'r-', 'LineWidth', 2);
legend('Vessel 1 Liquid Angle [deg]','Vessel 2 Liquid Angle [deg]', 'Location', 'best');
xlabel('Time [s]');
ylabel('Liquid Angle [deg]');
title('Improved Pendulum Liquid Oscillation');
grid on;

subplot(2,1,2);
plot(t, sol(:,5), 'k-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Plate Position [m]');
title('Plate Movement due to Sloshing');
grid on;

% ODE system function with improvements
function dx = coupled_vessel_slosh_improved(t, x, params)
    % Extract parameters
    g = params.g;
    mass1 = params.mass1;
    mass2 = params.mass2;
    l1 = params.l1;
    l2 = params.l2;
    I1 = params.I1;
    I2 = params.I2;
    k1 = params.k1;
    k2 = params.k2;
    c1 = params.c1;
    c2 = params.c2;
    plate_mass = params.plate_mass;
    plate_damp = params.plate_damp;
    alpha_impact = params.alpha_impact;
    beta_plate = params.beta_plate;

    % State variables
    theta1 = x(1);
    dtheta1 = x(2);
    theta2 = x(3);
    dtheta2 = x(4);
    x_plate = x(5);
    dx_plate = x(6);

    % IMPROVED: More realistic wall impact model
    % Only activate when angle is large (near walls) and moving toward wall
    wall_threshold = deg2rad(45); % activate near walls
    if abs(theta1) > wall_threshold
        wall_impact1 = alpha_impact * mass1 * l1 * dtheta1^2 * sign(theta1);
    else
        wall_impact1 = 0;
    end

    % Plate movement gets force from wall impact
    ddx_plate = (wall_impact1 - plate_damp * dx_plate) / plate_mass;

    % Vessel 2 base follows plate
    base_acc2 = beta_plate * ddx_plate;

    % IMPROVED: Vessel 1 angular acceleration (small angle approximation for realistic behavior)
    if abs(theta1) < deg2rad(30) % small angle approximation
        ddtheta1 = (-c1*dtheta1 - k1*theta1) / I1;
    else % large angle - nonlinear
        ddtheta1 = (-c1*dtheta1 - k1*theta1 + mass1*g*l1*sin(theta1)) / I1;
    end

    % Vessel 2 angular acceleration (excitation by plate movement)
    force_plate2 = base_acc2 * l2 * mass2;
    if abs(theta2) < deg2rad(30)
        ddtheta2 = (-c2*dtheta2 - k2*theta2 + force_plate2) / I2;
    else
        ddtheta2 = (-c2*dtheta2 - k2*theta2 + mass2*g*l2*sin(theta2) + force_plate2) / I2;
    end

    % Return derivatives
    dx = zeros(6,1);
    dx(1) = dtheta1;
    dx(2) = ddtheta1;
    dx(3) = dtheta2;
    dx(4) = ddtheta2;
    dx(5) = dx_plate;
    dx(6) = ddx_plate;
end