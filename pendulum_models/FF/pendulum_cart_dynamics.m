function dxdt = pendulum_cart_dynamics(t, x, F_func, p)
    % Unpack parameters
    m = p.m;
    M = p.M;
    L = p.L;
    g = p.g;
    Kd = p.Kd;

    % States
    x_pos = x(1);
    x_dot = x(2);
    theta = x(3);
    theta_dot = x(4);

    % Input force
    F = F_func(t);

    % Useful terms
    S = sin(theta);
    C = cos(theta);

    denom = M + m - m * C^2;

    % Accelerations
    x_ddot = (F + m * L * theta_dot^2 * S - m * g * C * S) / denom;
    theta_ddot = (-F * C - m * L * theta_dot^2 * S * C + (M + m) * g * S) / (L * denom);

    % State derivatives
    dxdt = [x_dot;
            x_ddot;
            theta_dot;
            theta_ddot];
end
