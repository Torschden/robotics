classdef WaiterMotionOptimizer < handle
    properties
        n_dof
        n_obj
        n_states
        n_controls
        
        % Robot parameters
        q_dot_max
        u_max
        tau_max
        
        % Object parameters
        m_obj
        r_o
        mu_0
        M_z_max
        
        % Optimization parameters
        N
        k
        
        % CasADi symbolic variables
        x
        q
        q_dot
        q_ddot
        u
        dt
        
        % Robot kinematic parameters (simplified 6-DOF robot)
        dh_params  % DH parameters for forward kinematics
    end
    
    methods
        function self = WaiterMotionOptimizer(n_dof, n_obj)
            if nargin < 1, n_dof = 6; end
            if nargin < 2, n_obj = 2; end
            
            self.n_dof = n_dof;
            self.n_obj = n_obj;
            self.n_states = 3 * n_dof; % q, q_dot, q_ddot
            self.n_controls = n_dof;   % jerk
            
            % Robot limits
            self.q_dot_max = [2.0, 2.0, 2.0, 3.0, 3.0, 3.0];
            self.u_max = [10.0, 10.0, 10.0, 15.0, 15.0, 15.0];
            self.tau_max = [100.0, 100.0, 80.0, 40.0, 40.0, 20.0];
            
            % Object parameters
            self.m_obj = 0.12;
            self.r_o = 0.025;
            self.mu_0 = 0.7;
            self.M_z_max = (2/3) * self.mu_0 * self.r_o^3;
            
            % Optimization setup
            self.N = 150;
            self.k = 1e-3;
            
            % Initialize simplified DH parameters (example for 6-DOF robot)
            % [a, alpha, d, theta_offset] for each joint
            self.dh_params = [
                0,     pi/2,   0.5,   0;      % Joint 1
                0.5,   0,      0,     0;      % Joint 2
                0.3,   pi/2,   0,     0;      % Joint 3
                0,     -pi/2,  0.4,   0;      % Joint 4
                0,     pi/2,   0,     0;      % Joint 5
                0,     0,      0.1,   0;      % Joint 6
            ];
            
            % Setup symbolic variables
            self.setup_symbolic_variables();
        end
        
        function setup_symbolic_variables(self)
            import casadi.*
            self.x = SX.sym('x', self.n_states, 1);
            self.q = self.x(1:self.n_dof);
            self.q_dot = self.x(self.n_dof+1:2*self.n_dof);
            self.q_ddot = self.x(2*self.n_dof+1:end);
            
            self.u = SX.sym('u', self.n_controls, 1);
            self.dt = SX.sym('dt');
        end
        
        function T = dh_transform(self, a, alpha, d, theta)
            % Create DH transformation matrix
            import casadi.*
            T = [cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta);
                 sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
                 0,           sin(alpha),             cos(alpha),            d;
                 0,           0,                      0,                     1];
        end
        
        function T = forward_kinematics(self, q)
            % Forward kinematics using DH parameters
            import casadi.*
            T = DM.eye(4);
            
            for i = 1:self.n_dof
                a = self.dh_params(i, 1);
                alpha = self.dh_params(i, 2);
                d = self.dh_params(i, 3);
                theta_offset = self.dh_params(i, 4);
                theta = q(i) + theta_offset;
                
                T_i = self.dh_transform(a, alpha, d, theta);
                T = T * T_i;
            end
        end
        
        function [pos, euler] = extract_pose(self, T)
            % Extract position and Euler angles from transformation matrix
            import casadi.*
            pos = T(1:3, 4);
            
            % Extract ZYX Euler angles (roll-pitch-yaw)
            sy = sqrt(T(1,1)^2 + T(2,1)^2);
            singular = sy < 1e-6;
            
            if_else = @(cond, if_true, if_false) if_true*cond + if_false*(1-cond);
            
            x = if_else(singular, 0, atan2(T(3,2), T(3,3)));
            y = if_else(singular, atan2(-T(3,1), sy), atan2(-T(3,1), sy));
            z = if_else(singular, atan2(T(2,1), T(1,1)), atan2(T(2,1), T(1,1)));
            
            euler = [x; y; z];  % [roll, pitch, yaw]
        end
        
        function q_solution = inverse_kinematics(self, pos_des, euler_des, q_init)
            % Solve inverse kinematics using optimization
            import casadi.*
            
            if nargin < 4
                q_init = zeros(self.n_dof, 1);  % Default initial guess
            end
            
            % Setup optimization for IK
            q = SX.sym('q', self.n_dof, 1);
            T = self.forward_kinematics(q);
            [pos, euler] = self.extract_pose(T);
            
            % Cost function: minimize pose error
            pos_error = pos - pos_des(:);
            euler_error = euler - euler_des(:);
            
            % Wrap angle differences to [-pi, pi]
            for i = 1:3
                euler_error(i) = atan2(sin(euler_error(i)), cos(euler_error(i)));
            end
            
            cost = sum1(pos_error.^2) + sum1(euler_error.^2);
            
            % Joint limits
            q_min = [-pi; -pi/2; -pi; -pi; -pi/2; -pi];
            q_max = [pi; pi/2; pi; pi; pi/2; pi];
            
            % Create NLP
            nlp = struct('x', q, 'f', cost);
            opts = struct;
            opts.ipopt.print_level = 0;
            opts.print_time = false;
            
            solver = nlpsol('ik_solver', 'ipopt', nlp, opts);
            
            % Solve
            sol = solver('x0', q_init, 'lbx', q_min, 'ubx', q_max);
            q_solution = full(sol.x);
            
            % Check if solution is good
            if solver.stats.return_status ~= "Solve_Succeeded"
                warning('Inverse kinematics may not have converged properly');
            end
        end
        
        function x_dot = dynamics(self, x, u)
            q_dot = x(self.n_dof+1:2*self.n_dof);
            q_ddot = x(2*self.n_dof+1:end);
            q_dddot = u;
            x_dot = vertcat(q_dot, q_ddot, q_dddot);
        end
        
        function tau = inverse_dynamics_simplified(self, q, q_dot, q_ddot)
            import casadi.*
            M_diag = DM([10.0; 8.0; 6.0; 2.0; 1.5; 1.0]);  
            tau = M_diag .* q_ddot + 0.1 .* abs(q_dot);
        end
        
        function [f_x, f_y, f_z, M_x, M_y, M_z] = contact_forces(self, q, q_dot, q_ddot, obj_idx)
            import casadi.*
            if obj_idx == 0
                r_obj = DM([0.185, 0.075, 0.025]);
            else
                r_obj = DM([0.185, -0.075, 0.025]);
            end
            
            % Use forward kinematics to get proper end-effector acceleration
            T = self.forward_kinematics(q);
            
            % Simplified acceleration calculation
            % In practice, this should use proper Jacobian computation
            a_ee = q_ddot(1:3);  % Simplified
            
            f_x = -self.m_obj * a_ee(1);
            f_y = -self.m_obj * a_ee(2);
            f_z = -self.m_obj * (9.81 + a_ee(3));
            
            M_x = -self.m_obj * r_obj(3) * a_ee(2);
            M_y =  self.m_obj * r_obj(3) * a_ee(1);
            M_z = 0.0;
        end
        
        function [nlp, lbx, ubx, lbg, ubg, X, U, T] = setup_optimization_problem(self, q0, qe)
            import casadi.*
            
            X = SX.sym('X', self.n_states, self.N+1);
            U = SX.sym('U', self.n_controls, self.N);
            T = SX.sym('T');
            
            % Objective
            J = T + self.k * sum(sum(U.^2)) * (T/self.N);
            
            g = {};
            lbg = [];
            ubg = [];
            
            % Initial constraints
            g{end+1} = X(1:self.n_dof,1) - q0(:);
            lbg = [lbg; zeros(self.n_dof,1)];
            ubg = [ubg; zeros(self.n_dof,1)];
            
            g{end+1} = X(self.n_dof+1:2*self.n_dof,1);
            lbg = [lbg; zeros(self.n_dof,1)];
            ubg = [ubg; zeros(self.n_dof,1)];
            
            g{end+1} = X(2*self.n_dof+1:end,1);
            lbg = [lbg; zeros(self.n_dof,1)];
            ubg = [ubg; zeros(self.n_dof,1)];
            
            % Final constraints
            g{end+1} = X(1:self.n_dof,end) - qe(:);
            lbg = [lbg; zeros(self.n_dof,1)];
            ubg = [ubg; zeros(self.n_dof,1)];
            
            g{end+1} = X(self.n_dof+1:2*self.n_dof,end);
            lbg = [lbg; zeros(self.n_dof,1)];
            ubg = [ubg; zeros(self.n_dof,1)];
            
            g{end+1} = X(2*self.n_dof+1:end,end);
            lbg = [lbg; zeros(self.n_dof,1)];
            ubg = [ubg; zeros(self.n_dof,1)];
            
            % Dynamics constraints (multiple shooting)
            dt = T / self.N;
            for k = 1:self.N
                x_k = X(:,k);
                u_k = U(:,k);
                x_k1 = X(:,k+1);
                
                k1 = self.dynamics(x_k, u_k);
                k2 = self.dynamics(x_k + dt/2*k1, u_k);
                k3 = self.dynamics(x_k + dt/2*k2, u_k);
                k4 = self.dynamics(x_k + dt*k3, u_k);
                
                x_next = x_k + dt/6*(k1 + 2*k2 + 2*k3 + k4);
                
                g{end+1} = x_k1 - x_next;
                lbg = [lbg; zeros(self.n_states,1)];
                ubg = [ubg; zeros(self.n_states,1)];
                
                % Joint velocity limits
                q_dot_k = x_k(self.n_dof+1:2*self.n_dof);
                for i = 1:self.n_dof
                    g{end+1} = q_dot_k(i);
                    lbg = [lbg; -self.q_dot_max(i)];
                    ubg = [ubg;  self.q_dot_max(i)];
                end
                
                % Control limits
                for i = 1:self.n_controls
                    g{end+1} = u_k(i);
                    lbg = [lbg; -self.u_max(i)];
                    ubg = [ubg;  self.u_max(i)];
                end
                
                % Torque limits
                q_k = x_k(1:self.n_dof);
                q_ddot_k = x_k(2*self.n_dof+1:end);
                tau = self.inverse_dynamics_simplified(q_k, q_dot_k, q_ddot_k);
                for i = 1:self.n_dof
                    g{end+1} = tau(i);
                    lbg = [lbg; -self.tau_max(i)];
                    ubg = [ubg;  self.tau_max(i)];
                end
                
                % Contact constraints
                for obj_i = 0:self.n_obj-1
                    safe_sqrt = @(z) sqrt(z + 1e-12);
                    
                    [f_x, f_y, f_z, M_x, M_y, M_z] = self.contact_forces(q_k, q_dot_k, q_ddot_k, obj_i);
                    
                    % Non-lifting
                    g{end+1} = f_z; lbg = [lbg; -inf]; ubg = [ubg; 0];
                    
                    % Non-slipping
                    g{end+1} = safe_sqrt(f_x^2 + f_y^2) + self.mu_0 * f_z;
                    lbg = [lbg; -inf]; ubg = [ubg; 0];
                    
                    % Non-tip-over
                    g{end+1} = safe_sqrt(M_x^2 + M_y^2) + self.r_o * f_z;
                    lbg = [lbg; -inf]; ubg = [ubg; 0];
                    
                    % Non-twisting
                    g{end+1} = abs(M_z) + self.M_z_max * f_z;
                    lbg = [lbg; -inf]; ubg = [ubg; 0];
                end
            end
            
            g = vertcat(g{:});
            
            % Bounds
            x_min = -inf(self.n_states,1);
            x_max =  inf(self.n_states,1);
            u_min = -self.u_max(:);
            u_max =  self.u_max(:);
            
            lbx = [];
            ubx = [];
            for k = 1:self.N+1
                lbx = [lbx; x_min];
                ubx = [ubx; x_max];
            end
            for k = 1:self.N
                lbx = [lbx; u_min];
                ubx = [ubx; u_max];
            end
            lbx = [lbx; 0.1]; % T min
            ubx = [ubx; 10.0]; % T max
            
            % NLP
            decision_vars = vertcat( ...
                reshape(X, self.n_states*(self.N+1), 1), ...
                reshape(U, self.n_controls*self.N, 1), ...
                T);
            nlp = struct('x', decision_vars, 'f', J, 'g', g);
        end
        
        function [X_opt, U_opt, T_opt, stats] = solve_optimization_cartesian(self, r0, rE, phi0, phiE)
            % Main solve function using Cartesian coordinates
            import casadi.*
            
            fprintf('Converting Cartesian poses to joint configurations...\n');
            
            % Convert Cartesian poses to joint configurations using IK
            q0 = self.inverse_kinematics(r0, phi0);
            qe = self.inverse_kinematics(rE, phiE, q0); % Use q0 as initial guess for qe
            
            % Verify the IK solutions
            fprintf('Initial pose IK solution:\n');
            T0 = self.forward_kinematics(q0);
            [pos0_check, euler0_check] = self.extract_pose(T0);
            fprintf('  Desired pos: [%.3f, %.3f, %.3f], Got: [%.3f, %.3f, %.3f]\n', ...
                    r0(1), r0(2), r0(3), full(pos0_check(1)), full(pos0_check(2)), full(pos0_check(3)));
            fprintf('  Desired ori: [%.3f, %.3f, %.3f], Got: [%.3f, %.3f, %.3f]\n', ...
                    phi0(1), phi0(2), phi0(3), full(euler0_check(1)), full(euler0_check(2)), full(euler0_check(3)));
            
            fprintf('Final pose IK solution:\n');
            Te = self.forward_kinematics(qe);
            [pose_check, eulere_check] = self.extract_pose(Te);
            fprintf('  Desired pos: [%.3f, %.3f, %.3f], Got: [%.3f, %.3f, %.3f]\n', ...
                    rE(1), rE(2), rE(3), full(pose_check(1)), full(pose_check(2)), full(pose_check(3)));
            fprintf('  Desired ori: [%.3f, %.3f, %.3f], Got: [%.3f, %.3f, %.3f]\n', ...
                    phiE(1), phiE(2), phiE(3), full(eulere_check(1)), full(eulere_check(2)), full(eulere_check(3)));
            
            % Solve optimization with joint configurations
            [X_opt, U_opt, T_opt, stats] = self.solve_optimization(q0, qe);
        end
        
        function [X_opt, U_opt, T_opt, stats] = solve_optimization(self, q0, qe)
            import casadi.*
            
            [nlp, lbx, ubx, lbg, ubg, X, U, T] = self.setup_optimization_problem(q0, qe);
            
            opts = struct;
            opts.ipopt.max_iter = 3000;
            opts.ipopt.print_level = 0;
            opts.print_time = false;
            
            solver = nlpsol('solver', 'ipopt', nlp, opts);
            
            % Initial guess
            x0_guess = [];
            for k = 0:self.N
                alpha = k/self.N;
                q_guess = (1-alpha)*q0(:) + alpha*qe(:);
                qd_guess = zeros(self.n_dof,1);
                qdd_guess = zeros(self.n_dof,1);
                x0_guess = [x0_guess; q_guess; qd_guess; qdd_guess];
            end
            for k = 1:self.N
                x0_guess = [x0_guess; zeros(self.n_controls,1)];
            end
            x0_guess = [x0_guess; 2.0]; % time guess
            
            fprintf('Solving optimization problem...\n');
            sol = solver('x0', x0_guess, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);
            
            x_opt = full(sol.x);
            n_x_vars = (self.N+1) * self.n_states;
            n_u_vars = self.N * self.n_controls;
            
            X_opt = reshape(x_opt(1:n_x_vars), self.n_states, self.N+1);
            U_opt = reshape(x_opt(n_x_vars+1:n_x_vars+n_u_vars), self.n_controls, self.N);
            T_opt = x_opt(end);
            
            stats = solver.stats();
        end
        
        function plot_results(self, X_opt, U_opt, T_opt)
            time = linspace(0, double(T_opt), self.N+1);
            time_u = linspace(0, double(T_opt), self.N);
        
            X_np = full(X_opt);
            U_np = full(U_opt);
        
            figure;
            tiledlayout(3,2,'TileSpacing','compact');
        
            % Joint positions
            nexttile;
            hold on;
            title('Joint Positions');
            for i = 1:self.n_dof
                plot(time, X_np(i,:), 'DisplayName', sprintf('q%d', i));
            end
            xlabel('Time [s]');
            ylabel('Joint angle [rad]');
            legend; grid on;
        
            % Joint velocities
            nexttile;
            hold on;
            title('Joint Velocities');
            for i = 1:self.n_dof
                plot(time, X_np(self.n_dof+i,:), 'DisplayName', sprintf('qd%d', i));
            end
            xlabel('Time [s]');
            ylabel('Joint velocity [rad/s]');
            legend; grid on;
        
            % Joint accelerations
            nexttile;
            hold on;
            title('Joint Accelerations');
            for i = 1:self.n_dof
                plot(time, X_np(2*self.n_dof+i,:), 'DisplayName', sprintf('qdd%d', i));
            end
            xlabel('Time [s]');
            ylabel('Joint acceleration [rad/s^2]');
            legend; grid on;
        
            % Joint jerks
            nexttile;
            hold on;
            title('Joint Jerks (Control inputs)');
            for i = 1:self.n_controls
                plot(time_u, U_np(i,:), 'DisplayName', sprintf('u%d', i));
            end
            xlabel('Time [s]');
            ylabel('Joint jerk [rad/s^3]');
            legend; grid on;
        
            % Velocity limits check
            nexttile;
            hold on;
            title('Normalized Joint Velocity Limits');
            for i = 1:self.n_dof
                plot(time, abs(X_np(self.n_dof+i,:))/self.q_dot_max(i), ...
                     'DisplayName', sprintf('|qd%d|/max', i));
            end
            yline(1,'r--','Limit');
            xlabel('Time [s]');
            ylabel('Normalized velocity [-]');
            legend; grid on;
        
            % Control limits check
            nexttile;
            hold on;
            title('Normalized Jerk (Control) Limits');
            for i = 1:self.n_controls
                plot(time_u, abs(U_np(i,:))/self.u_max(i), ...
                     'DisplayName', sprintf('|u%d|/max', i));
            end
            yline(1,'r--','Limit');
            xlabel('Time [s]');
            ylabel('Normalized jerk [-]');
            legend; grid on;
        
            sgtitle('Waiter Motion Optimization Results');
            fprintf('Optimal execution time: %.3f seconds\n', double(T_opt));
        end
        
        function plot_end_effector_trajectory(self, X_opt)
            % Plot the end-effector trajectory in Cartesian space
            time = linspace(0, size(X_opt,2)-1, size(X_opt,2));
            
            % Extract joint positions
            Q = X_opt(1:self.n_dof, :);
            
            % Compute end-effector positions
            positions = zeros(3, size(Q,2));
            orientations = zeros(3, size(Q,2));
            
            for k = 1:size(Q,2)
                T = self.forward_kinematics(Q(:,k));
                [pos, euler] = self.extract_pose(T);
                positions(:,k) = full(pos);
                orientations(:,k) = full(euler);
            end
            
            figure;
            subplot(2,1,1);
            hold on;
            plot(time, positions(1,:), 'r-', 'LineWidth', 2, 'DisplayName', 'X');
            plot(time, positions(2,:), 'g-', 'LineWidth', 2, 'DisplayName', 'Y');
            plot(time, positions(3,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Z');
            xlabel('Time step');
            ylabel('Position [m]');
            title('End-Effector Position');
            legend; grid on;
            
            subplot(2,1,2);
            hold on;
            plot(time, orientations(1,:), 'r-', 'LineWidth', 2, 'DisplayName', 'Roll');
            plot(time, orientations(2,:), 'g-', 'LineWidth', 2, 'DisplayName', 'Pitch');
            plot(time, orientations(3,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Yaw');
            xlabel('Time step');
            ylabel('Orientation [rad]');
            title('End-Effector Orientation');
            legend; grid on;
            
            % 3D trajectory
            figure;
            plot3(positions(1,:), positions(2,:), positions(3,:), 'b-', 'LineWidth', 2);
            hold on;
            plot3(positions(1,1), positions(2,1), positions(3,1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
            plot3(positions(1,end), positions(2,end), positions(3,end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
            xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
            title('3D End-Effector Trajectory');
            legend('Trajectory', 'Start', 'End');
            grid on; axis equal;
        end
    end
end