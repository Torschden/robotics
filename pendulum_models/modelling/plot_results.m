function plot_results(t, theta, theta_dot, x_c, x_pendulum, y_pendulum, params)
    % Create analysis plots
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot 1: Pendulum angle vs time
    subplot(2,4,1);
    plot(t, theta*180/pi, 'b-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Pendulum Angle (degrees)');
    title('Pendulum Angle vs Time');
    grid on;
    
    % Plot 2: Cart position vs time
    subplot(2,4,2);
    plot(t, x_c, 'r-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Cart Position (m)');
    title('Cart Position vs Time');
    grid on;
    
    % Plot 3: Pendulum angular velocity
    subplot(2,4,3);
    plot(t, theta_dot, 'g-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    title('Pendulum Angular Velocity');
    grid on;
    
    % Plot 4: Phase portrait
    subplot(2,4,4);
    plot(theta*180/pi, theta_dot, 'b-', 'LineWidth', 1.5);
    xlabel('Angle (degrees)');
    ylabel('Angular Velocity (rad/s)');
    title('Phase Portrait');
    grid on;
    
    % Plot 5: Energy analysis
    plot_energy_analysis(t, theta, theta_dot, params);
    
    % Plot 6: Trajectory plot
    subplot(2,4,6);
    plot(x_pendulum, y_pendulum, 'b-', 'LineWidth', 1.5);
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Pendulum Trajectory');
    axis equal;
    grid on;
    
    % Plot 7: Pendulum length verification
    subplot(2,4,7);
    pendulum_length = sqrt((x_pendulum - x_c').^2 + y_pendulum.^2);
    plot(t, pendulum_length, 'r-', 'LineWidth', 2);
    hold on;
    plot(t, params.L * ones(size(t)), 'k--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Length (m)');
    title('Pendulum Length Verification');
    legend('Actual', 'Expected', 'Location', 'best');
    grid on;
    length_error = max(abs(pendulum_length - params.L));
    ylim([params.L - max(0.01, 2*length_error), params.L + max(0.01, 2*length_error)]);
    
    % Add text showing maximum error
    text(0.05, 0.95, sprintf('Max Error: %.2e m', length_error), ...
         'Units', 'normalized', 'FontSize', 8, 'BackgroundColor', 'white');
    
    % Plot 8: Position components
    subplot(2,4,8);
    plot(t, x_pendulum, 'r-', t, y_pendulum, 'b-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Position (m)');
    title('Pendulum Position Components');
    legend('X Position', 'Y Position', 'Location', 'best');
    grid on;
end

function plot_energy_analysis(t, theta, theta_dot, params)
    % Plot energy analysis in subplot 5
    KE = 0.5 * params.m_p * (params.L * theta_dot).^2;  % kinetic energy
    PE = params.m_p * params.g * params.L * (1 - cos(theta)); % potential energy
    Total_E = KE + PE;
    
    subplot(2,4,5);
    plot(t, KE, 'r-', t, PE, 'b-', t, Total_E, 'k--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Energy (J)');
    title('Energy vs Time');
    legend('Kinetic', 'Potential', 'Total', 'Location', 'best');
    grid on;
end