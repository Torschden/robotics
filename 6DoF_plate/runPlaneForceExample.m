function runPlaneForceExample()
    % Create simulation instance
    sim = PlaneForceSimulation();

    sim.setPlanePose([0; 0; 0], [0; 0; 0]);
    sim.addPoint([2; 1; 1]);
    
    % Apply constant force and torque
    goal_force = [5; 0; 0];
    goal_torque = [0; 0; 3];
    sim.simulateMotion(0.05, 60, [zeros(3,1), goal_force], [zeros(3,1), goal_torque]);
    
    fprintf('\nSimulation complete!\n');

end