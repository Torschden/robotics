function sloshing_simulation()
    % Sloshing simulation based on Godderidge et al. (2012)
    % "A rapid method for the simulation of sloshing using a mathematical model 
    % based on the pendulum equation"
    
    clear; close all; clc;
    
    %% Tank and Fluid Parameters
    L = 1;          % Tank length (m)
    B = 1;           % Tank width (m)
    h = 1;           % Fluid depth (m)
    rho = 1025;      % Fluid density (kg/m³)
    g = 9.81;        % Gravitational acceleration (m/s²)
    mu = 1e-3;       % Dynamic viscosity (Pa·s)
    
    % Filling ratio
    h_L_ratio = h/L;
    
    %% Pendulum Model Parameters
    % Natural frequency calculation
    omega_n = sqrt(g * pi * tanh(pi * h/L) / L);
    T1 = 2*pi/omega_n;  % Natural period
    
    % Effective mass (first mode)
    m_total = rho * L * B * h;
    m_eff = m_total * (8 * L * tanh(pi * h/L)) / (pi^3 * h);
    
    % Pendulum length
    l = g / omega_n^2;
    
    % Restoring force coefficient (from paper)
    alpha = 1.0;  % Adjusted based on CFD data in practice
    
    %% Damping Parameters
    % Calculate damping coefficient using Keulegan's approach
    k = pi/L;  % Wave number
    omega = omega_n;
    
    % Boundary layer dissipation components
    DE1 = (g^2 / omega^2 * pi^2 / 4) * sqrt(mu * rho / (2 * omega)) * sinh(2*k*h) / cosh(k*h)^2;
    DE2 = (g^2 / omega^2 * pi/2) * sqrt(mu * rho / (2 * omega)) * (B * k / cosh(k*h)^2) * (sinh(2*k*h)/2 - k*h);
    DE3 = (g^2 / omega^2 * pi^2 / 2) * sqrt(mu * rho / (2 * omega)) * (B * k / cosh(k*h)^2);
    
    % Total energy
    E = pi * rho * g * B / (4 * k);   %% wave amplitude?  (5) 
    
    % Linear damping coefficient
    d_linear = (DE1 + DE2 + DE3) / (2 * E);
    
    % Third-order damping coefficient
    d3 = d_linear * 0.1;  % Empirical scaling from paper
    
    %% Impact Parameters
    theta_crit = 15 * pi/180;  % Critical angle for impact (radians)
    theta_1 = 20 * pi/180;     % Upper bound for impact force
    a_impact = 50;             % Impact force coefficient
    b_impact = 23;             % Impact force exponent (2N-1)
    c_impact = 10;             % Impact damping coefficient
    d_impact = 24;             % Impact damping exponent
    m_impact = 1;              % Impact damping power
    

    T_exc = 1.1 * T1;   % Excitation period
    omega_exc = 2*pi / T_exc;
    
    %% Simulation Parameters
    t_end = 200;         % Simulation time (s)
    dt = 0.01;          % Time step (s)
    t_vec = 0:dt:t_end;
    N = length(t_vec);
    
    %% Excitation Parameters
    % Translational excitation
    x0_vec = zeros(size(t_vec));
    x0_vec(t_vec <= 5) = 0.5;

    % for i = 0:t_end
    %     if i <= 5
    %         x0(i) = 0.5;
    %     else 
    %         x0(i) = 0;
    %     end
    % end

    %% Initial Conditions
    theta0 = 0;         % Initial angle (rad)
    theta_dot0 = 0;     % Initial angular velocity (rad/s)
    
    %% Solve ODE System
    % State vector: [theta, theta_dot]
    y0 = [theta0; theta_dot0];
    
    % Define the system of ODEs
    x0_interp = @(t) interp1(t_vec, x0_vec, t, 'previous', 0);  % zero extrapolation

    odefun = @(t, y) sloshing_ode(t, y, l, alpha, d_linear, d3, g, ...
                                  theta_crit, theta_1, a_impact, b_impact, ...
                                  c_impact, d_impact, m_impact, ...
                                  x0_interp(t), omega_exc);
    
    % Solve using ODE113 (Adams-Bashforth-Moulton method as in paper)
    options = odeset('RelTol', 3e-6, 'AbsTol', 1e-8);
    [t_sol, y_sol] = ode113(odefun, t_vec, y0, options);
    
    theta = y_sol(:, 1);
    theta_dot = y_sol(:, 2);
    
    %% Calculate Forces and Moments
    % Translational excitation
    A_t = x0_vec * omega_exc^2 * sin(omega_exc * t_sol);
    
    % Sloshing force
    F_sloshing = m_eff * A_t .* cos(theta);
    
    % Restoring force
    F_restoring = -m_eff * g * alpha * sin(theta) / l;
    
    % Impact forces
    F_impact = zeros(size(theta));
    for i = 1:length(theta)
        if abs(theta(i)) > theta_crit && abs(theta(i)) < theta_1
            F_impact(i) = -sign(theta(i)) * a_impact * (abs(theta(i))/theta_crit)^b_impact;
        end
    end
    
    %% Plotting
    figure('Position', [100, 100, 1200, 800]);
    
    % Subplot 1: Angular displacement
    subplot(2, 3, 1);
    plot(t_sol, theta * 180/pi, 'b-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Angular Displacement (deg)');
    title('Fluid Centre of Mass Angular Displacement');
    grid on;
    
    % Subplot 2: Angular velocity
    subplot(2, 3, 2);
    plot(t_sol, theta_dot, 'r-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    title('Angular Velocity');
    grid on;
    
    % Subplot 3: Phase portrait
    subplot(2, 3, 3);
    plot(theta * 180/pi, theta_dot, 'g-', 'LineWidth', 1.5);
    xlabel('Angular Displacement (deg)');
    ylabel('Angular Velocity (rad/s)');
    title('Phase Portrait');
    grid on;
    
    % Subplot 4: Sloshing force
    subplot(2, 3, 4);
    plot(t_sol, F_sloshing/1000, 'b-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Sloshing Force (kN)');
    title('Sloshing Force');
    grid on;
    
    % Subplot 5: Impact forces
    subplot(2, 3, 5);
    plot(t_sol, F_impact/1000, 'r-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Impact Force (kN)');
    title('Impact Forces');
    grid on;
    
    % Subplot 6: Energy analysis
    subplot(2, 3, 6);
    KE = 0.5 * m_eff * (l * theta_dot).^2;  % Kinetic energy
    PE = m_eff * g * l * (1 - cos(theta));   % Potential energy
    Total_E = KE + PE;
    
    plot(t_sol, KE/1000, 'b-', 'LineWidth', 1.5); hold on;
    plot(t_sol, PE/1000, 'r-', 'LineWidth', 1.5);
    plot(t_sol, Total_E/1000, 'k--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Energy (kJ)');
    title('Energy Analysis');
    legend('Kinetic', 'Potential', 'Total', 'Location', 'best');
    grid on;
    
    sgtitle('Sloshing Simulation Results - Pendulum Model');
    
    %% Display Parameters
    fprintf('=== Sloshing Simulation Parameters ===\n');
    fprintf('Tank Length (L): %.2f m\n', L);
    fprintf('Tank Width (B): %.2f m\n', B);
    fprintf('Fluid Depth (h): %.2f m\n', h);
    fprintf('Filling Ratio (h/L): %.3f\n', h_L_ratio);
    fprintf('Natural Period (T1): %.2f s\n', T1);
    fprintf('Excitation Period: %.2f s (%.2f T1)\n', T_exc, T_exc/T1);
    fprintf('Effective Mass: %.0f kg\n', m_eff);
    fprintf('Pendulum Length: %.2f m\n', l);
    fprintf('Linear Damping: %.6f\n', d_linear);
    fprintf('Critical Impact Angle: %.1f deg\n', theta_crit * 180/pi);
    
    %% Additional Analysis Plot
    figure('Position', [150, 150, 1000, 600]);
    
    % Time series comparison
    subplot(2, 2, 1);
    plot(t_sol, theta * 180/pi, 'b-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    title('Angular Displacement vs Time');
    grid on;
    
    % Excitation vs response
    subplot(2, 2, 2);
    excitation = x0_vec * sin(omega_exc * t_sol);
    plot(t_sol, excitation, 'r-', 'LineWidth', 1.5); hold on;
    plot(t_sol, theta * 180/pi / 10, 'b-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Excitation vs Response');
    legend('Excitation (m)', 'Response (deg/10)', 'Location', 'best');
    grid on;
    
    % Frequency response (simple)
    subplot(2, 2, 3);
    % Find peaks to estimate response amplitude
    [pks, locs] = findpeaks(abs(theta), 'MinPeakHeight', 0.01);
    if ~isempty(pks)
        response_amp = mean(pks) * 180/pi
        excitation_amp = 0.5
        transfer_function = response_amp / excitation_amp;
        
        bar([1, 2], [excitation_amp, response_amp]);
        set(gca, 'XTickLabel', {'Excitation (m)', 'Response (deg)'});
        ylabel('Amplitude');
        title(sprintf('Response Amplitude (TF: %.2f deg/m)', transfer_function));
    end
    grid on;
    
    % Damping analysis
    subplot(2, 2, 4);
    % Calculate envelope of oscillations
    env_upper = abs(hilbert(theta));
    plot(t_sol, env_upper * 180/pi, 'r-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Envelope (deg)');
    title('Oscillation Envelope (Damping Effect)');
    grid on;
    
    sgtitle('Detailed Sloshing Analysis');
end

function dydt = sloshing_ode(t, y, l, alpha, d_linear, d3, g, ...
                            theta_crit, theta_1, a_impact, b_impact, ...
                            c_impact, d_impact, m_impact, ...
                            x0, omega_exc)
    % ODE system for sloshing pendulum model
    % y(1) = theta (angular displacement)
    % y(2) = theta_dot (angular velocity)
    
    theta = y(1);
    theta_dot = y(2);
    
    % Translational excitation
    A_t = x0 * omega_exc^2 * sin(omega_exc * t);
    
    % Restoring force term
    restoring = -(g/l) * alpha * sin(theta);
    
    % Linear damping
    linear_damping = -d_linear * theta_dot;
    
    % Third-order damping
    nonlinear_damping = -d3 * theta_dot^3;
    
    % External forcing
    external_force = (A_t / l) * cos(theta);
    
    % Impact force
    impact_force = 0;
    if abs(theta) > theta_crit && abs(theta) < theta_1
        impact_force = -sign(theta) * (a_impact/l) * (abs(theta)/theta_crit)^b_impact;
    end
    
    % Impact damping
    impact_damping = 0;
    if abs(theta) > theta_crit
        impact_damping = -sign(theta_dot) * (c_impact/l) * ...
                        (abs(theta)/theta_crit)^(2*d_impact) * abs(theta_dot)^m_impact;
    end
    
    % System of ODEs
    dydt = zeros(2, 1);
    dydt(1) = theta_dot;
    dydt(2) = restoring + linear_damping + nonlinear_damping + external_force + ...
              impact_force + impact_damping;
end