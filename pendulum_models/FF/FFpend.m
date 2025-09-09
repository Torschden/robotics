% Parameters
M = 1;      % cart mass (kg)
m = 1;      % pendulum mass (kg)
g = 9.81;   % gravity (m/s^2)
l = 0.5;    % pendulum length (m)
Kd = 10;    % cart damping coefficient

% Force profile struct
profile.times  = [3,3.5,4, 10, 10.5,11];
profile.values = [-10,17,7.5, 17, -10,0];
profile.k      = 25;      % logistic steepness
profile.t0     = 2;      % force zero before t0

F = make_logistic_force(profile);  % function from previous messages

% Initial state: [cart pos; cart vel; pendulum angle; pendulum ang vel]
x0 = [0; 0; 0; 0];

tspan = [0 30];

% Pack parameters in struct for ODE function
params.M = M;
params.m = m;
params.g = g;
params.l = l;
params.Kd = Kd;

% ODE options with tighter tolerances for stability
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

% Simulate
[t, x] = ode45(@(t,x) pendulum_cart_dynamics(t,x,F,params), tspan, x0, opts);

% Plot
figure;
subplot(4,1,1);
plot(t, x(:,1), 'LineWidth', 1.5);
ylabel('Cart Position z (m)');
grid on;

subplot(4,1,2);
plot(t, x(:,2), 'LineWidth', 1.5);
ylabel('Cart Velocity \dot{z} (m/s)');
grid on;

subplot(4,1,3);
plot(t, x(:,3)*180, 'LineWidth', 1.5);
ylabel('Pendulum Angle \theta (°)');
grid on;

subplot(4,1,4);
fplot(F, [0 30], 'r', 'LineWidth', 1.5);
ylabel('Force F(t) (N)');
xlabel('Time (s)');
grid on;

% --- Dynamics function ---
function dxdt = pendulum_cart_dynamics(t, x, F_func, p)
    z_dot = x(2);
    theta = x(3);
    theta_dot = x(4);
    F = F_func(t);

    % Shorthands
    M = p.M;
    m = p.m;
    g = p.g;
    l = p.l;
    Kd = p.Kd;

    b_theta = 0.1; % small pendulum angular damping coefficient

    denom1 = M + m * sin(theta)^2;
    denom2 = l - (m * l * cos(theta)^2)/(M + m);

    % Compute derivatives
    dzdt = z_dot;

    dz_dotdt = ( F ...
               - Kd * z_dot ...
               - m * l * theta_dot^2 * sin(theta) ...
               - m * g * sin(theta) * cos(theta) ) / denom1;

    dthetadt = theta_dot;

    dtheta_dotdt = ( -g * sin(theta) ...
                   + ( F - Kd * z_dot - m * l * theta_dot^2 * sin(theta) ) ...
                     * cos(theta) / (M + m) ) / denom2 ;

    dxdt = [dzdt; dz_dotdt; dthetadt; dtheta_dotdt];
end

% --- Logistic force generator ---
function F_handle = make_logistic_force(profile)
    times = profile.times;
    values = profile.values;
    k = profile.k;
    t0 = profile.t0;

    F_handle = @(t) logistic_sum(t, times, values, k, t0);
end

function F = logistic_sum(t, times, values, k, t0)
    t = t(:);
    F = zeros(size(t));
    F_current = 0;

    for i = 1:length(times)
        F_next = values(i);
        delta = F_next - F_current;
        F = F + delta ./ (1 + exp(-k * (t - times(i))));
        F_current = F_next;
    end

    F(t < t0) = 0;
end

% xp = x(:,1) + l * sin(x(:,3)); % pendulum tip x-position
% yp = -l * cos(x(:,3));         % pendulum tip y-position (down is negative)
% figure;
% plot(xp, yp);
% axis equal;
% grid on;
% xlabel('Pendulum tip x (m)');
% ylabel('Pendulum tip y (m)');
% title('Pendulum tip trajectory');