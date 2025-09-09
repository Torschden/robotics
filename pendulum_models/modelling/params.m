function p = params()
    % Define all system parameters
    p.g = 9.81;           % gravity (m/s^2)
    p.L = 1.0;            % pendulum length (m)
    p.m_p = 0.1;          % pendulum mass (kg)
    p.m_c = 1.0;          % cart mass (kg) - for visualization only
    p.b = 0.01;           % damping coefficient
    
    p.v_max = 0.3;
    p.a_max = 0.05;
    
    % Simulation parameters
    p.dt = 0.01;          % time step (s)
    p.t_end = 10;         % simulation time (s)
    
    % Initial conditions
    p.theta_0 = pi/6;     % initial angle (rad)
    p.theta_dot_0 = 0;    % initial angular velocity (rad/s)
    p.x_c_0 = 0;          % initial cart position (m)
    
    % Animation parameters
    p.anim_skip = 5;      % animate every nth frame
    p.trail_length = 200; % length of pendulum trail
end