% Plot results
function plotting(t,y,u)
    figure('Name', '3D Liquid Sloshing Simulation', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);
    
    %% position
    subplot(5,3,1);
    plot(t, y(:,1), 'LineWidth', 2); % y position
    xlabel('Time [s]'); ylabel('x_{l0} [m]');
    title('Liquid Position X');
    grid on;

    subplot(5,3,2);
    plot(t, y(:,2), 'LineWidth', 2); % y position
    xlabel('Time [s]'); ylabel('y_{l0} [m]');
    title('Liquid Position Y');
    grid on;
    
    subplot(5,3,3);
    plot(t, y(:,3), 'LineWidth', 2); % y position
    xlabel('Time [s]'); ylabel('z_{l0} [m]');
    title('Liquid Position Z');
    grid on;

    %% velocity
    subplot(5,3,4);
    plot(t, y(:,7), 'LineWidth', 2); % x velocity
    xlabel('Time [s]'); ylabel('ẋ_{l0} [m/s]');
    title('Liquid Velocity X');
    grid on;

    subplot(5,3,5);
    plot(t, y(:,8), 'LineWidth', 2); % y velocity
    xlabel('Time [s]'); ylabel('ẏ_{l0} [m/s]');
    title('Liquid Velocity Y');
    grid on;

    subplot(5,3,6);
    plot(t, y(:,9), 'LineWidth', 2); % x velocity
    xlabel('Time [s]'); ylabel('ż_{l0} [m/s]');
    title('Liquid Velocity Z');
    grid on;
    
    %% angle
    subplot(5,3,8);
    plot(t, y(:,5) * 180/pi, 'LineWidth', 2); % phi angle
    xlabel('Time [s]'); ylabel('θ_{l0} [deg]');
    title('Liquid Angle around y-axis');
    grid on;

    subplot(5,3,7);
    plot(t, y(:,4) * 180/pi, 'LineWidth', 2); % phi angle
    xlabel('Time [s]'); ylabel('θ_{l0} [deg]');
    title('Liquid Angle around y-axis');
    grid on;

    subplot(5,3,9);
    plot(t, y(:,6) * 180/pi, 'LineWidth', 2); % phi angle
    xlabel('Time [s]'); ylabel('θ_{l0} [deg]');
    title('Liquid Angle around y-axis');
    grid on;
    

    %% angular velocities
    subplot(5,3,10);
    plot(t, y(:,10) * 180/pi, 'LineWidth', 2); % phi angular velocity
    xlabel('Time [s]'); ylabel('φ̇_{l0} [deg/s]');
    title('Liquid Angular Velocity φ');
    grid on;

    subplot(5,3,11);
    plot(t, y(:,11) * 180/pi, 'LineWidth', 2); % phi angular velocity
    xlabel('Time [s]'); ylabel('φ̇_{l0} [deg/s]');
    title('Liquid Angular Velocity φ');
    grid on;

    subplot(5,3,12);
    plot(t, y(:,12) * 180/pi, 'LineWidth', 2); % phi angular velocity
    xlabel('Time [s]'); ylabel('φ̇_{l0} [deg/s]');
    title('Liquid Angular Velocity φ');
    grid on;
    
    %% Control input
    t_ctrl = linspace(0, t(end), 1000);
    u_vals = arrayfun(@(t) u(t), t_ctrl, 'UniformOutput', false);
    u_x = cellfun(@(u) u(1), u_vals);
    u_y = cellfun(@(u) u(2), u_vals);
    u_z = cellfun(@(u) u(3), u_vals);

    subplot(5,3,13);
    plot(t_ctrl, u_x, 'LineWidth', 2);
    xlabel('Time [s]'); ylabel('u_x [m/s²]');
    ylim([min(u_x)-0.2, max(u_x)+0.2]);
    title('Control Input X');
    grid on;

    subplot(5,3,14);
    plot(t_ctrl, u_y, 'LineWidth', 2);
    xlabel('Time [s]'); ylabel('u_y [m/s²]');
    ylim([min(u_y)-0.2, max(u_y)+0.2]);
    title('Control Input Y');
    grid on;

    subplot(5,3,15);
    plot(t_ctrl, u_z, 'LineWidth', 2);
    xlabel('Time [s]'); ylabel('u_z [m/s²]');
    ylim([min(u_z)-0.2, max(u_z)+0.2]);
    title('Control Input Z');
    grid on;

    sgtitle('3D Liquid Sloshing Simulation Results');
end