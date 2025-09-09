function main()
    clear 
    close all
    clc

    fprintf('Starting Cart-Pendulum Simulation...\n');

    p = params();
    t = 0:p.dt:p.t_end;
    
    [theta, theta_dot, x_c, x_pendulum, y_pendulum] = simulation(t,p);

    visualize(t, theta, theta_dot, x_c, x_pendulum, y_pendulum, p);

    fprintf('Simulation Complete!\n');
end