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
            
            a_ee = q_ddot(1:3);
            
            f_x = -self.m_obj * a_ee(1);
            f_y = -self.m_obj * a_ee(2);
            f_z = -self.m_obj * (9.81 + a_ee(3));
            
            M_x = -self.m_obj * r_obj(3) * a_ee(2);
            M_y =  self.m_obj * r_obj(3) * a_ee(1);
            M_z = 0.0;
        end
        
        %% setup
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
        
        %% solve
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
            
            disp('Solving optimization problem...');
            sol = solver('x0', x0_guess, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);
            
            x_opt = full(sol.x);
            n_x_vars = (self.N+1) * self.n_states;
            n_u_vars = self.N * self.n_controls;
            
            X_opt = reshape(x_opt(1:n_x_vars), self.n_states, self.N+1);
            U_opt = reshape(x_opt(n_x_vars+1:n_x_vars+n_u_vars), self.n_controls, self.N);
            T_opt = x_opt(end);
            
            stats = solver.stats();
        end

        %% plotting
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
                plot(time, X_np(self.n_dof+i,:), 'DisplayName', sprintf('q̇%d', i));
            end
            xlabel('Time [s]');
            ylabel('Joint velocity [rad/s]');
            legend; grid on;
        
            % Joint accelerations
            nexttile;
            hold on;
            title('Joint Accelerations');
            for i = 1:self.n_dof
                plot(time, X_np(2*self.n_dof+i,:), 'DisplayName', sprintf('q̈%d', i));
            end
            xlabel('Time [s]');
            ylabel('Joint acceleration [rad/s^2]');
            legend; grid on;
        
            % Joint jerks
            nexttile;
            hold on;
            title('Joint Jerks (Control inputs)');
            for i = 1:self.n_controls
                plot(time_u, U_np(i,:), 'DisplayName', sprintf('q⃛%d', i));
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
                     'DisplayName', sprintf('|q̇%d|/max', i));
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
    end
end