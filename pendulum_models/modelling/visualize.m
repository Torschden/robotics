function visualize(t, theta, theta_dot, x_c, x_pendulum, y_pendulum, p)
    plot_results(t, theta, theta_dot, x_c, x_pendulum, y_pendulum, p);
    
    % animation
    animate_system(t, theta, x_c, x_pendulum, y_pendulum, p);
    
    % Display
    display_results(theta, x_c);
end

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
    
    % % Plot 7: Pendulum length verification
    % subplot(2,4,7);
    % pendulum_length = sqrt((x_pendulum - x_c').^2 + y_pendulum.^2);
    % plot(t, pendulum_length, 'r-', 'LineWidth', 2);
    % hold on;
    % plot(t, params.L * ones(size(t)), 'k--', 'LineWidth', 1.5);
    % xlabel('Time (s)');
    % ylabel('Length (m)');
    % title('Pendulum Length Verification');
    % legend('Actual', 'Expected', 'Location', 'best');
    % grid on;
    % length_error = max(abs(pendulum_length - params.L));
    % ylim([params.L - max(0.01, 2*length_error), params.L + max(0.01, 2*length_error)]);
    % 
    % % Add text showing maximum error
    % text(0.05, 0.95, sprintf('Max Error: %.2e m', length_error), ...
    %      'Units', 'normalized', 'FontSize', 8, 'BackgroundColor', 'white');
    
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

function animate_system(t, theta, x_c, x_pendulum, y_pendulum, params)
    % Create and run animation
    figure('Position', [200, 150, 800, 600]);
    setup_animation_figure();
    
    n_steps = length(t);
    fprintf('Starting animation...\n');
    
    for i = 1:params.anim_skip:n_steps
        draw_frame(i, t, theta, x_c, x_pendulum, y_pendulum, params);
        drawnow;
        pause(0.01);
    end
end

function setup_animation_figure()
    % Set up the animation figure
    axis equal;
    xlim([-2, 2]);
    ylim([-1.5, 0.5]);
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Cart-Pendulum System Animation');
    grid on;
    hold on;
end

function draw_frame(i, t, theta, x_c, x_pendulum, y_pendulum, params)
    % Draw a single animation frame
    cla;
    
    % Draw cart
    draw_cart(x_c(i));
    
    % Draw pendulum
    draw_pendulum(x_c(i), x_pendulum(i), y_pendulum(i));
    
    % Draw trail
    draw_trail(i, x_pendulum, y_pendulum, params);
    
    % Draw ground
    plot([-2, 2], [-0.3, -0.3], 'k-', 'LineWidth', 2);
    
    % Update display
    update_display(i, t, theta, x_c);
    
    % Set axis properties
    axis equal;
    xlim([-2, 2]);
    ylim([-1.5, 0.5]);
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title(sprintf('Cart-Pendulum System - Time: %.2f s', t(i)));
    grid on;
end

function draw_cart(x_pos)
    % Draw the cart
    cart_width = 0.2;
    cart_height = 0.1;
    rectangle('Position', [x_pos-cart_width/2, -cart_height/2, cart_width, cart_height], ...
              'FaceColor', 'blue', 'EdgeColor', 'black');
end

function draw_pendulum(x_cart, x_pend, y_pend)
    plot([x_cart, x_pend], [0, y_pend], 'k-', 'LineWidth', 3);
    plot(x_pend, y_pend, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red');
    plot(x_cart, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'black');
end

function draw_trail(i, x_pendulum, y_pendulum, params)
    % Draw pendulum trail
    trail_start = max(1, i - params.trail_length);
    if i > trail_start
        plot(x_pendulum(trail_start:i), y_pendulum(trail_start:i), ...
             'b:', 'LineWidth', 1, 'Color', [0.5, 0.5, 1, 0.7]);
    end
end

function update_display(i, t, theta, x_c)
    % Update text display on animation
    text(-1.8, 0.3, sprintf('Time: %.2f s', t(i)), 'FontSize', 12, 'FontWeight', 'bold');
    text(-1.8, 0.1, sprintf('Angle: %.1f°', theta(i)*180/pi), 'FontSize', 10);
    text(-1.8, -0.1, sprintf('Cart Pos: %.2f m', x_c(i)), 'FontSize', 10);
end

function display_results(theta, x_c)
    % Display final simulation results
    fprintf('\nSimulation Results:\n');
    fprintf('Maximum pendulum angle: %.2f degrees\n', max(abs(theta))*180/pi);
    fprintf('Final pendulum angle: %.2f degrees\n', theta(end)*180/pi);
    fprintf('Cart displacement range: %.2f m\n', max(x_c) - min(x_c));
end