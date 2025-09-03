% sim = PlaneForceSimulation();
% 
% % Set plane pose
% sim.setPlanePose([1; 2; 1], [0.2; 0.3; 0.1]); % position and orientation
% 
% % Add points above plane
% sim.addPoint([2; 1; 1]);      % Point in body frame
% % sim.addPoint([-1; 2; 0.5]);
% 
% % Apply force in X and torque around Z
% sim.setForceAndTorque([15; 0; 0], [0; 0; 8]);
% 
% % Calculate and visualize
% sim.updateVisualization();
% sim.printResults();

% % Run example
% runPlaneForceExample();

sim = PlaneAccelerationSimulation();

sim.setPlanePose([0; 0; 0], [0; 0; 0]);
sim.addPoint([2; 1; 1]);

% Apply constant force and torque
goal_force = [5; 0; 0];
goal_torque = [0; 0; 3];
sim.simulateMotion(0.05, 60, [zeros(3,1), goal_force], [zeros(3,1), goal_torque]);

fprintf('\nSimulation complete!\n');