%% 3D Liquid Sloshing Model Implementation
% Based on "A Dynamic Optimization Approach for Sloshing Free Transport
% of Liquid Filled Containers using an Industrial Robot"
% Reinhold et al., 2019

classdef SloshingModel3D < handle
    properties
        % Physical parameters
        m           % Mass of liquid [kg]
        l_phi       % Pendulum length for phi rotation [m]
        l_theta     % Pendulum length for theta rotation [m]  
        d_phi       % Damping constant for phi [N*s*m^-1]
        d_theta     % Damping constant for theta [N*s*m^-1]
        g           % Gravitational constant [m/s^2]
        
        % Container properties
        V           % Container volume [m^3]
        rho_l       % Liquid density [kg/m^3]
        
        % Model parameters for cylindrical container
        omega       % Natural frequency [rad/s]
        zeta        % Damping ratio [-]
        f_natural   % Natural frequency [Hz]
    end
    
    methods
        function obj = SloshingModel3D(container_params)
            % Constructor with container parameters
            % container_params: struct with fields V, rho_l, (optional: m, l, d)
            
            obj.g = 9.81; % [m/s^2]
            obj.V = container_params.V;
            obj.rho_l = container_params.rho_l;
            obj.m = obj.rho_l * obj.V;

            % Use provided parameters or defaults from paper
            if isfield(container_params, 'l')
                obj.l_phi = container_params.l;
                obj.l_theta = container_params.l;
            else
                obj.l_phi = 0.021;   % [m] from Table I
                obj.l_theta = 0.021; % [m] from Table I
            end
            
            if isfield(container_params, 'd')
                obj.d_phi = container_params.d;
                obj.d_theta = container_params.d;
            else
                obj.d_phi = 1.51;    % [N*s*m^-1] from Table I
                obj.d_theta = 1.51;  % [N*s*m^-1] from Table I
            end
            
            % Calculate derived parameters
            obj.omega = sqrt(obj.g / obj.l_phi);
            obj.zeta = obj.d_phi / (2 * obj.m * obj.omega);
            obj.f_natural = obj.omega / (2 * pi);
            
            fprintf('Sloshing Model Parameters:\n');
            fprintf('Mass: %.3f kg\n', obj.m);
            fprintf('Pendulum length: %.3f m\n', obj.l_phi);
            fprintf('Damping constant: %.3f N*s*m^-1\n', obj.d_phi);
            fprintf('Natural frequency: %.3f Hz\n', obj.f_natural);
            fprintf('Damping ratio: %.3f\n', obj.zeta);
        end
        
        function dydt = dynamics3D(obj, t, y, u_interp)
            % 3D liquid dynamics according to equation (3)
            % y = [x_l0, y_l0, z_l0, phi_l0, theta_l0, psi_l0, 
            %      x_dot_l0, y_dot_l0, z_dot_l0, phi_dot_l0, theta_dot_l0, psi_dot_l0]
            % u = [u_x, u_y, u_z] - container accelerations
            
            % Interpolate control input at current time
            u = u_interp(t);
            u_x = u(1); u_y = u(2); u_z = u(3);
            
            % Extract states
            x_l0 = y(1); y_l0 = y(2); z_l0 = y(3);
            phi_l0 = y(4); theta_l0 = y(5); psi_l0 = y(6);
            x_dot_l0 = y(7); y_dot_l0 = y(8); z_dot_l0 = y(9);
            phi_dot_l0 = y(10); theta_dot_l0 = y(11); psi_dot_l0 = y(12);
            
            % State derivatives according to equation (3)
            dydt = zeros(12, 1);
            
            % Position derivatives
            dydt(1) = x_dot_l0;      % dx_l0/dt
            dydt(2) = y_dot_l0;      % dy_l0/dt  
            dydt(3) = z_dot_l0;      % dz_l0/dt
            dydt(4) = phi_dot_l0;    % dphi_l0/dt
            dydt(5) = theta_dot_l0;  % dtheta_l0/dt
            dydt(6) = psi_dot_l0;    % dpsi_l0/dt
            
            % Velocity derivatives
            dydt(7) = u_x;           % d²x_l0/dt²
            dydt(8) = u_y;           % d²y_l0/dt²
            dydt(9) = u_z;           % d²z_l0/dt²
            
            % Angular acceleration for phi (rotation about x-axis, excited by u_y)
            dydt(10) = (u_y / obj.l_phi) * cos(phi_l0) - ...
                      ((obj.g - u_z) / obj.l_phi) * sin(phi_l0) - ...
                      (obj.d_phi / obj.m) * phi_dot_l0 * cos(phi_l0)^2;
            
            % Angular acceleration for theta (rotation about y-axis, excited by u_x)  
            dydt(11) = (u_x / obj.l_theta) * cos(theta_l0) - ...
                      ((obj.g - u_z) / obj.l_theta) * sin(theta_l0) - ...
                      (obj.d_theta / obj.m) * theta_dot_l0 * cos(theta_l0)^2;
            
            % No rotation around z-axis
            dydt(12) = 0;            % d²psi_l0/dt²
        end
        
        function [t, y] = simulate(obj, t_span, y0, u_func)
            % Simulate the 3D liquid dynamics
            % t_span: time vector [t0, t1] or [t0, t1, ..., tn]
            % y0: initial state vector (12x1)
            % u_func: function handle for control input u(t) -> [u_x, u_y, u_z]
            
            % Create interpolation function for control input
            if isa(u_func, 'function_handle')
                u_interp = u_func;
            else
                error('u_func must be a function handle');
            end
            
            % ODE options
            options = odeset('RelTol', 1e-6, 'AbsTol', 1e-7);
            
            % Solve ODE
            [t, y] = ode45(@(t, y) obj.dynamics3D(t, y, u_interp), t_span, y0, options);
        end
        
        function H = getTransferFunction(obj)
            % Get transfer function for linearized model (equation 5)
            % H(s) = phi(s)/u(s) for small angle approximation
            
            s = tf('s');
            K = 1/obj.g * obj.l_phi;
            omega_n = obj.omega;
            zeta_n = obj.zeta;
            
            H = K / (1/omega_n^2 * s^2 + 2*zeta_n/omega_n * s + 1);
        end
        
        function response = impulseResponse(obj, t_sim)
            % Analytical impulse response for PT2 element (equation 7)
            % Assuming unit impulse amplitude
            
            if obj.zeta < 1  % Underdamped
                omega_d = obj.omega * sqrt(1 - obj.zeta^2);
                A = 1; % Unit amplitude
                response = (A / sqrt(1 - obj.zeta^2)) * exp(-obj.zeta * obj.omega * t_sim) .* ...
                          sin(omega_d * t_sim);
            else  % Overdamped or critically damped
                response = zeros(size(t_sim));
                warning('System is not underdamped, impulse response set to zero');
            end
        end
        
        function [omega_id, zeta_id] = identifyParameters(obj, t_exp, phi_exp, plot_results)
            % Identify parameters from experimental impulse response
            % t_exp: experimental time vector
            % phi_exp: experimental angle response
            % plot_results: boolean to show plots
            
            if nargin < 4
                plot_results = false;
            end
            
            % FFT to find natural frequency
            Fs = 1 / mean(diff(t_exp));  % Sampling frequency
            Y = fft(phi_exp);
            P2 = abs(Y/length(phi_exp));
            P1 = P2(1:floor(length(phi_exp)/2)+1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = Fs*(0:(length(phi_exp)/2))/length(phi_exp);
            
            [~, max_idx] = max(P1(2:end)); % Skip DC component
            f_peak = f(max_idx + 1);
            omega_id = 2 * pi * f_peak;
            
            % Optimize damping ratio using least squares
            A_exp = max(abs(phi_exp));
            
            objective = @(zeta) obj.optimizeDamping(zeta, omega_id, A_exp, t_exp, phi_exp);
            zeta_id = fminbnd(objective, 0.001, 0.999);
            
            if plot_results
                figure;
                
                % Time domain comparison
                subplot(2,1,1);
                plot(t_exp, phi_exp, 'r-', 'LineWidth', 2);
                hold on;
                phi_model = (A_exp / sqrt(1 - zeta_id^2)) * exp(-zeta_id * omega_id * t_exp) .* ...
                           sin(omega_id * sqrt(1 - zeta_id^2) * t_exp);
                plot(t_exp, phi_model, 'b--', 'LineWidth', 2);
                xlabel('Time [s]');
                ylabel('Angle [rad]');
                title('Impulse Response Comparison');
                legend('Experimental', 'Model', 'Location', 'best');
                grid on;
                
                % Frequency domain
                subplot(2,1,2);
                semilogx(f, P1);
                xlabel('Frequency [Hz]');
                ylabel('|FFT(φ)|');
                title('Power Spectral Density');
                grid on;
                hold on;
                plot(f_peak, P1(max_idx + 1), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
                legend('FFT', sprintf('Peak at %.3f Hz', f_peak), 'Location', 'best');
            end
            
            fprintf('Identified Parameters:\n');
            fprintf('Natural frequency: %.3f Hz (%.3f rad/s)\n', f_peak, omega_id);
            fprintf('Damping ratio: %.4f\n', zeta_id);
        end
        
        function runImpulsDemo(obj)
            % Impulse response simulation
            fprintf('\n=== Impulse Response Test ===\n');
            t_sim = 0:0.01:2;
            impulse_resp = obj.impulseResponse(t_sim);

            figure('Name', 'Impulse Response', 'NumberTitle', 'off');
            plot(t_sim, impulse_resp * 180/pi, 'b-', 'LineWidth', 2);
            xlabel('Time [s]');
            ylabel('Angle [degrees]');
            title('Theoretical Impulse Response');
            grid on;
        end

        function y = runDemo(obj,u)
            % 3D simulation with step input
            fprintf('\n=== 3D Dynamics Simulation ===\n');
            
            % Initial conditions (rest state)
            y0 = zeros(12, 1);
            y0(3) = 0.2;  % Initial height z_l0 = 0.2m
            
            % Time span
            t_span = [0, 7];

            % Simulate
            fprintf('Running 3D simulation...\n');
            [t, y] = obj.simulate(t_span, y0, u);
            
            % plot
            plotting(t,y,u);    

            liquid_surface_animator(t, y)
            
            fprintf('Simulation completed. Maximum sloshing angle: %.2f degrees\n',max(abs(y(:,4))) * 180/pi);
            fprintf('\nSimulation Statistics:\n');
            fprintf('Final liquid position (y): %.4f m\n', y(end, 2));
            fprintf('Final sloshing angle (φ): %.4f degrees\n', y(end, 4) * 180/pi);
            fprintf('Maximum liquid displacement: %.4f m\n', max(abs(y(:, 2))));
            fprintf('Time to settle (within 1°): %.2f s\n', find_settle_time(t, y(:,4), 1*pi/180));
        end
    end
    
    methods (Access = private)
        function error = optimizeDamping(obj, zeta, omega, A, t_exp, phi_exp)
            % Objective function for damping parameter optimization
            if zeta >= 1 || zeta <= 0
                error = inf;
                return;
            end
            
            phi_model = (A / sqrt(1 - zeta^2)) * exp(-zeta * omega * t_exp) .* ...
                       sin(omega * sqrt(1 - zeta^2) * t_exp);
            error = sum((phi_exp - phi_model').^2);
        end
    end
end

%% Helper Functions
function t_settle = find_settle_time(t, phi, threshold)
    % Find time when oscillation settles within threshold
    phi_abs = abs(phi);
    settled_idx = find(phi_abs <= threshold, 1, 'first');
    if ~isempty(settled_idx)
        % Check if it stays settled for at least 0.5 seconds
        end_check = min(length(t), settled_idx + round(0.5/(t(2)-t(1))));
        if all(phi_abs(settled_idx:end_check) <= threshold)
            t_settle = t(settled_idx);
        else
            t_settle = t(end);  % Not settled by end of simulation
        end
    else
        t_settle = t(end);  % Not settled by end of simulation
    end
end

%% Example Usage - Uncomment the line below to run demo automatically
% SloshingModel3D.demo();