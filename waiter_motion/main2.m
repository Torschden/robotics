function waiterMotionOptimalControl()
    % Waiter Motion Optimal Control Problem - MATLAB Implementation
    % Based on "Point to point time optimal handling of unmounted rigid objects"
    % by Gattringer et al. (2023)
    
    % Problem parameters
    params = setupParameters();
    
    % Initial and final configurations (from paper)
    q0 = [1.56; 0.57; -0.20; 0; -0.36; 0];  % rad
    qe = [-1.56; 0.57; -0.20; 0; -0.36; 0]; % rad
    
    % Solve optimization problem
    [solution, fval, exitflag] = solveOptimalControl(params, q0, qe);
    
    % Plot results
    if exitflag > 0
        plotResults(solution, params);
        fprintf('Optimization completed successfully!\n');
        fprintf('Optimal execution time: %.3f seconds\n', solution.T_opt);
    else
        fprintf('Optimization failed with exit flag: %d\n', exitflag);
    end
end

function params = setupParameters()
    % Setup optimization parameters
    
    % Robot parameters
    params.n_dof = 6;           % Degrees of freedom
    params.n_obj = 2;           % Number of objects
    params.N = 150;             % Number of shooting intervals
    
    % Physical limits
    params.q_dot_max = [2.0; 2.0; 2.0; 3.0; 3.0; 3.0];     % rad/s
    params.u_max = [10.0; 10.0; 10.0; 15.0; 15.0; 15.0];   % rad/s^3
    params.tau_max = [100.0; 100.0; 80.0; 40.0; 40.0; 20.0]; % Nm
    
    % Object parameters
    params.m_obj = 0.12;        % kg (mass of each cup)
    params.r_o = 0.025;         % m (footprint radius)
    params.mu_0 = 0.7;          % static friction coefficient
    params.M_z_max = (2/3) * params.mu_0 * params.r_o^3; % max twisting torque
    
    % Optimization parameters
    params.k = 1e-6;            % smoothness weight (small for time-optimal)
    
    % State and control dimensions
    params.n_states = 3 * params.n_dof;    % [q; q_dot; q_ddot]
    params.n_controls = params.n_dof;       % jerk
end

function [solution, fval, exitflag] = solveOptimalControl(params, q0, qe)
    % Solve the optimal control problem using fmincon
    
    % Decision variables: [X(:); U(:); T]
    % X: (n_states x N+1) states at each node
    % U: (n_controls x N) controls at each interval  
    % T: final time
    
    n_x_vars = params.n_states * (params.N + 1);
    n_u_vars = params.n_controls * params.N;
    n_vars = n_x_vars + n_u_vars + 1;
    
    % Initial guess
    x0 = initialGuess(params, q0, qe);
    
    % Bounds
    [lb, ub] = setupBounds(params);
    
    % Constraints
    [A, b, Aeq, beq] = setupLinearConstraints(params, q0, qe);
    
    % Nonlinear constraints
    nonlcon = @(x) nonlinearConstraints(x, params);
    
    % Objective function
    objective = @(x) objectiveFunction(x, params);
    
    % Optimization options
    options = optimoptions('fmincon', ...
        'Algorithm', 'interior-point', ...
        'MaxIterations', 3000, ...
        'MaxFunctionEvaluations', 50000, ...
        'OptimalityTolerance', 1e-6, ...
        'ConstraintTolerance', 1e-6, ...
        'Display', 'iter');
    
    % Solve optimization
    fprintf('Solving waiter motion optimization problem...\n');
    [x_opt, fval, exitflag] = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    
    % Parse solution
    solution = parseSolution(x_opt, params);
end

function x0 = initialGuess(params, q0, qe)
    % Generate initial guess for optimization variables
    
    n_x_vars = params.n_states * (params.N + 1);
    n_u_vars = params.n_controls * params.N;
    
    x0 = zeros(n_x_vars + n_u_vars + 1, 1);
    
    % State initial guess (linear interpolation for positions, zero for vel/acc)
    for k = 0:params.N
        idx_start = k * params.n_states + 1;
        alpha = k / params.N;
        
        % Position (linear interpolation)
        q_guess = (1 - alpha) * q0 + alpha * qe;
        x0(idx_start:idx_start + params.n_dof - 1) = q_guess;
        
        % Velocity and acceleration (zero)
        % Already initialized to zero
    end
    
    % Control initial guess (zero)
    % Already initialized to zero
    
    % Time initial guess
    x0(end) = 2.0;
end

function [lb, ub] = setupBounds(params)
    % Setup variable bounds
    
    n_x_vars = params.n_states * (params.N + 1);
    n_u_vars = params.n_controls * params.N;
    n_vars = n_x_vars + n_u_vars + 1;
    
    lb = -inf(n_vars, 1);
    ub = inf(n_vars, 1);
    
    % Control bounds
    u_start_idx = n_x_vars + 1;
    for k = 1:params.N
        for i = 1:params.n_controls
            idx = u_start_idx + (k-1) * params.n_controls + i - 1;
            lb(idx) = -params.u_max(i);
            ub(idx) = params.u_max(i);
        end
    end
    
    % Time bounds
    lb(end) = 0.1;  % minimum time
    ub(end) = 10.0; % maximum time
end

function [A, b, Aeq, beq] = setupLinearConstraints(params, q0, qe)
    % Setup linear equality constraints for boundary conditions
    
    n_x_vars = params.n_states * (params.N + 1);
    n_u_vars = params.n_controls * params.N;
    n_vars = n_x_vars + n_u_vars + 1;
    
    % No linear inequality constraints
    A = [];
    b = [];
    
    % Boundary conditions (equality constraints)
    n_eq = 6 * params.n_dof; % 3 conditions at start, 3 at end
    Aeq = zeros(n_eq, n_vars);
    beq = zeros(n_eq, 1);
    
    eq_idx = 1;
    
    % Initial conditions: q(0) = q0, q_dot(0) = 0, q_ddot(0) = 0
    for i = 1:params.n_dof
        % q(0) = q0(i)
        Aeq(eq_idx, i) = 1;
        beq(eq_idx) = q0(i);
        eq_idx = eq_idx + 1;
    end
    
    for i = 1:params.n_dof
        % q_dot(0) = 0
        Aeq(eq_idx, params.n_dof + i) = 1;
        beq(eq_idx) = 0;
        eq_idx = eq_idx + 1;
    end
    
    for i = 1:params.n_dof
        % q_ddot(0) = 0
        Aeq(eq_idx, 2*params.n_dof + i) = 1;
        beq(eq_idx) = 0;
        eq_idx = eq_idx + 1;
    end
    
    % Final conditions: q(T) = qe, q_dot(T) = 0, q_ddot(T) = 0
    final_idx = params.N * params.n_states;
    
    for i = 1:params.n_dof
        % q(T) = qe(i)
        Aeq(eq_idx, final_idx + i) = 1;
        beq(eq_idx) = qe(i);
        eq_idx = eq_idx + 1;
    end
    
    for i = 1:params.n_dof
        % q_dot(T) = 0
        Aeq(eq_idx, final_idx + params.n_dof + i) = 1;
        beq(eq_idx) = 0;
        eq_idx = eq_idx + 1;
    end
    
    for i = 1:params.n_dof
        % q_ddot(T) = 0
        Aeq(eq_idx, final_idx + 2*params.n_dof + i) = 1;
        beq(eq_idx) = 0;
        eq_idx = eq_idx + 1;
    end
end

function J = objectiveFunction(x, params)
    % Objective function: minimize time + smoothness penalty
    
    T = x(end);
    
    % Extract controls
    n_x_vars = params.n_states * (params.N + 1);
    U = reshape(x(n_x_vars+1:n_x_vars + params.n_controls*params.N), params.n_controls, params.N);
    
    % Time-optimal with smoothness penalty
    smoothness_penalty = params.k * sum(U(:).^2) * (T / params.N);
    J = T + smoothness_penalty;
end

function [c, ceq] = nonlinearConstraints(x, params)
    % Nonlinear constraints
    
    % Parse decision variables
    [X, U, T] = parseDecisionVars(x, params);
    
    c = [];   % Inequality constraints
    ceq = []; % Equality constraints
    
    dt = T / params.N;
    
    % Dynamics constraints (multiple shooting)
    for k = 1:params.N
        x_k = X(:, k);
        u_k = U(:, k);
        x_k1 = X(:, k+1);
        
        % RK4 integration
        x_next = rk4_step(x_k, u_k, dt, params);
        
        % Continuity constraint
        ceq = [ceq; x_k1 - x_next];
        
        % State constraints
        q_k = x_k(1:params.n_dof);
        qd_k = x_k(params.n_dof+1:2*params.n_dof);
        qdd_k = x_k(2*params.n_dof+1:3*params.n_dof);
        
        % Joint velocity limits: |q_dot| <= q_dot_max
        for i = 1:params.n_dof
            c = [c; qd_k(i) - params.q_dot_max(i)];
            c = [c; -qd_k(i) - params.q_dot_max(i)];
        end
        
        % Motor torque limits
        tau = inverseDynamicsSimplified(q_k, qd_k, qdd_k, params);
        for i = 1:params.n_dof
            c = [c; tau(i) - params.tau_max(i)];
            c = [c; -tau(i) - params.tau_max(i)];
        end
        
        % Contact constraints for each object
        for obj_i = 1:params.n_obj
            [f_x, f_y, f_z, M_x, M_y, M_z] = contactForces(q_k, qd_k, qdd_k, obj_i, params);
            
            % Non-lifting: f_z <= 0
            c = [c; f_z];
            
            % Non-slipping: sqrt(f_x^2 + f_y^2) <= -mu*f_z
            c = [c; sqrt(f_x^2 + f_y^2) + params.mu_0 * f_z];
            
            % Non-tip-over: sqrt(M_x^2 + M_y^2) <= -r_o*f_z
            c = [c; sqrt(M_x^2 + M_y^2) + params.r_o * f_z];
            
            % Non-twisting: |M_z| <= M_z_max*(-f_z)
            c = [c; abs(M_z) + params.M_z_max * f_z];
        end
    end
end

function [X, U, T] = parseDecisionVars(x, params)
    % Parse decision variables from optimization vector
    
    n_x_vars = params.n_states * (params.N + 1);
    n_u_vars = params.n_controls * params.N;
    
    X = reshape(x(1:n_x_vars), params.n_states, params.N + 1);
    U = reshape(x(n_x_vars+1:n_x_vars+n_u_vars), params.n_controls, params.N);
    T = x(end);
end

function x_next = rk4_step(x, u, dt, params)
    % RK4 integration step for dynamics
    
    k1 = dynamics(x, u, params);
    k2 = dynamics(x + dt/2*k1, u, params);
    k3 = dynamics(x + dt/2*k2, u, params);
    k4 = dynamics(x + dt*k3, u, params);
    
    x_next = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

function x_dot = dynamics(x, u, params)
    % System dynamics: integrator chain
    
    q = x(1:params.n_dof);
    q_dot = x(params.n_dof+1:2*params.n_dof);
    q_ddot = x(2*params.n_dof+1:3*params.n_dof);
    q_dddot = u;
    
    x_dot = [q_dot; q_ddot; q_dddot];
end

function tau = inverseDynamicsSimplified(q, q_dot, q_ddot, params)
    % Simplified inverse dynamics model
    % In practice, this would use identified robot parameters
    
    % Simplified diagonal mass matrix and friction
    M_diag = [10.0; 8.0; 6.0; 2.0; 1.5; 1.0];
    tau = M_diag .* q_ddot + 0.1 * abs(q_dot);
end

function [f_x, f_y, f_z, M_x, M_y, M_z] = contactForces(q, q_dot, q_ddot, obj_idx, params)
    % Compute contact forces for object
    % Simplified model - in practice would use full kinematics
    
    % Object positions relative to EE
    if obj_idx == 1
        r_obj = [0.185; 0.075; 0.025];
    else
        r_obj = [0.185; -0.075; 0.025];
    end
    
    % Simplified EE acceleration (would use forward kinematics in practice)
    a_ee = q_ddot(1:3);
    
    % Contact forces
    f_x = -params.m_obj * a_ee(1);
    f_y = -params.m_obj * a_ee(2);
    f_z = -params.m_obj * (9.81 + a_ee(3));
    
    % Contact moments
    M_x = -params.m_obj * r_obj(3) * a_ee(2);
    M_y = params.m_obj * r_obj(3) * a_ee(1);
    M_z = 0; % Simplified
end

function solution = parseSolution(x_opt, params)
    % Parse optimization solution
    
    [X_opt, U_opt, T_opt] = parseDecisionVars(x_opt, params);
    
    solution.X = X_opt;
    solution.U = U_opt;
    solution.T_opt = T_opt;
    solution.time = linspace(0, T_opt, params.N + 1);
    solution.time_u = linspace(0, T_opt, params.N);
end

function plotResults(solution, params)
    % Plot optimization results
    
    X = solution.X;
    U = solution.U;
    time = solution.time;
    time_u = solution.time_u;
    
    figure('Position', [100, 100, 1200, 800]);
    
    % Joint positions
    subplot(3, 2, 1);
    hold on;
    for i = 1:params.n_dof
        plot(time, X(i, :), 'LineWidth', 1.5, 'DisplayName', sprintf('q_%d', i));
    end
    xlabel('Time [s]');
    ylabel('Position [rad]');
    title('Joint Positions');
    legend('Location', 'best');
    grid on;
    
    % Joint velocities
    subplot(3, 2, 2);
    hold on;
    for i = 1:params.n_dof
        plot(time, X(params.n_dof + i, :), 'LineWidth', 1.5, 'DisplayName', sprintf('\\dot{q}_%d', i));
    end
    xlabel('Time [s]');
    ylabel('Velocity [rad/s]');
    title('Joint Velocities');
    legend('Location', 'best');
    grid on;
    
    % Joint accelerations
    subplot(3, 2, 3);
    hold on;
    for i = 1:params.n_dof
        plot(time, X(2*params.n_dof + i, :), 'LineWidth', 1.5, 'DisplayName', sprintf('\\ddot{q}_%d', i));
    end
    xlabel('Time [s]');
    ylabel('Acceleration [rad/s²]');
    title('Joint Accelerations');
    legend('Location', 'best');
    grid on;
    
    % Joint jerks (controls)
    subplot(3, 2, 4);
    hold on;
    for i = 1:params.n_controls
        plot(time_u, U(i, :), 'LineWidth', 1.5, 'DisplayName', sprintf('u_%d', i));
    end
    xlabel('Time [s]');
    ylabel('Jerk [rad/s³]');
    title('Joint Jerks (Control)');
    legend('Location', 'best');
    grid on;
    
    % Velocity limits check
    subplot(3, 2, 5);
    hold on;
    for i = 1:params.n_dof
        vel_norm = abs(X(params.n_dof + i, :)) / params.q_dot_max(i);
        plot(time, vel_norm, 'LineWidth', 1.5, 'DisplayName', sprintf('|\\dot{q}_%d|/max', i));
    end
    plot(time, ones(size(time)), 'r--', 'LineWidth', 2, 'DisplayName', 'Limit');
    xlabel('Time [s]');
    ylabel('Normalized Velocity');
    title('Velocity Limits Check');
    legend('Location', 'best');
    grid on;
    
    % Control limits check
    subplot(3, 2, 6);
    hold on;
    for i = 1:params.n_controls
        u_norm = abs(U(i, :)) / params.u_max(i);
        plot(time_u, u_norm, 'LineWidth', 1.5, 'DisplayName', sprintf('|u_%d|/max', i));
    end
    plot(time_u, ones(size(time_u)), 'r--', 'LineWidth', 2, 'DisplayName', 'Limit');
    xlabel('Time [s]');
    ylabel('Normalized Control');
    title('Control Limits Check');
    legend('Location', 'best');
    grid on;
    
    sgtitle(sprintf('Waiter Motion Optimization Results (T_{opt} = %.3f s)', solution.T_opt));
end

% Additional utility functions for extended analysis
function analyzeContactForces(solution, params)
    % Analyze contact forces throughout the trajectory
    
    X = solution.X;
    time = solution.time;
    
    % Initialize force arrays
    forces = zeros(6, params.n_obj, length(time)); % [f_x, f_y, f_z, M_x, M_y, M_z]
    
    for k = 1:length(time)
        q_k = X(1:params.n_dof, k);
        qd_k = X(params.n_dof+1:2*params.n_dof, k);
        qdd_k = X(2*params.n_dof+1:3*params.n_dof, k);
        
        for obj_i = 1:params.n_obj
            [f_x, f_y, f_z, M_x, M_y, M_z] = contactForces(q_k, qd_k, qdd_k, obj_i, params);
            forces(:, obj_i, k) = [f_x; f_y; f_z; M_x; M_y; M_z];
        end
    end
    
    % Plot contact forces
    figure('Position', [150, 150, 1200, 600]);
    
    for obj_i = 1:params.n_obj
        % Friction forces
        subplot(2, params.n_obj, obj_i);
        hold on;
        f_tangent = squeeze(sqrt(forces(1, obj_i, :).^2 + forces(2, obj_i, :).^2));
        f_friction_limit = -params.mu_0 * squeeze(forces(3, obj_i, :));
        plot(time, f_tangent, 'b-', 'LineWidth', 2, 'DisplayName', '||f_T||');
        plot(time, f_friction_limit, 'r--', 'LineWidth', 2, 'DisplayName', '-\mu f_z');
        xlabel('Time [s]');
        ylabel('Force [N]');
        title(sprintf('Friction Constraint - Object %d', obj_i));
        legend('Location', 'best');
        grid on;
        
        % Normal force
        subplot(2, params.n_obj, params.n_obj + obj_i);
        plot(time, squeeze(forces(3, obj_i, :)), 'g-', 'LineWidth', 2);
        xlabel('Time [s]');
        ylabel('Force [N]');
        title(sprintf('Normal Force - Object %d', obj_i));
        grid on;
    end
    
    sgtitle('Contact Force Analysis');
end