classdef PlaneAccelerationSimulation < handle
    % 6-DOF Plane Acceleration Simulation
    % Simulates a plane with translational and rotational degrees of freedom
    % and calculates acceleration at points above the plane
    
    properties
        % Plane center position in Cartesian coordinates
        center_pos = [0; 0; 0]
        
        % Plane orientation (Euler angles: roll, pitch, yaw in radians)
        orientation = [0; 0; 0]
        
        % Rotation matrix from body frame to world frame
        R_body_to_world
        
        % Points above the plane (in body frame coordinates)
        points_body = []
        
        % Lever arms (from center to points in body frame)
        lever_arms = []
        
        % Applied forces and torques at center ONLY
        applied_force = [0; 0; 0]      % [Fx; Fy; Fz] in world frame
        applied_torque = [0; 0; 0]     % [Tx; Ty; Tz] in world frame
        
        % Physical properties
        mass = 10                      % Mass of the plane (kg)
        inertia = eye(3) * 2          % Inertia tensor (kg⋅m²)
        
        % Current accelerations at center
        linear_acceleration = [0; 0; 0]    % Linear acceleration at center
        angular_acceleration = [0; 0; 0]   % Angular acceleration at center
        
        % Current velocities (for simulation)
        velocity = [0; 0; 0]
        angular_velocity = [0; 0; 0]
        
        % Visualization properties
        fig_handle
        ax_handle
        plane_handle
        points_handle
        acceleration_handles
    end
    
    methods
        function obj = PlaneAccelerationSimulation()
            % Constructor
            obj.updateRotationMatrix();
            obj.setupVisualization();
        end
        
        function updateRotationMatrix(obj)
            % Update rotation matrix from Euler angles (ZYX convention)
            roll = obj.orientation(1);
            pitch = obj.orientation(2);
            yaw = obj.orientation(3);
            
            % Individual rotation matrices
            Rx = [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)];
            Ry = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)];
            Rz = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1];
            
            % Combined rotation matrix (ZYX convention)
            obj.R_body_to_world = Rz * Ry * Rx;
        end
        
        function setPlanePose(obj, position, orientation)
            % Set plane position and orientation
            % position: [x; y; z] in world coordinates
            % orientation: [roll; pitch; yaw] in radians
            obj.center_pos = position(:);
            obj.orientation = orientation(:);
            obj.updateRotationMatrix();
        end
        
        function setPlaneProperties(obj, mass, inertia)
            % Set physical properties of the plane
            % mass: scalar mass (kg)
            % inertia: 3x3 inertia tensor (kg⋅m²)
            obj.mass = mass;
            obj.inertia = inertia;
        end
        
        function addPoint(obj, point_body)
            % Add a point above the plane in body frame coordinates
            % point_body: [x; y; z] relative to plane center in body frame
            obj.points_body = [obj.points_body, point_body(:)];
            obj.lever_arms = [obj.lever_arms, point_body(:)];
        end
        
        function clearPoints(obj)
            % Clear all points
            obj.points_body = [];
            obj.lever_arms = [];
        end
        
        function setForceAndTorque(obj, force_world, torque_world)
            % Set applied force and torque at plane center ONLY
            % force_world: [Fx; Fy; Fz] in world frame
            % torque_world: [Tx; Ty; Tz] in world frame
            obj.applied_force = force_world(:);
            obj.applied_torque = torque_world(:);
            
            % Calculate accelerations at center
            obj.linear_acceleration = obj.applied_force / obj.mass;
            obj.angular_acceleration = obj.inertia \ obj.applied_torque;
        end
        
        function accelerations_world = calculatePointAccelerations(obj)
            % Calculate translational accelerations at each point due to applied force and torque at center
            % Returns:
            % accelerations_world: 3xN matrix of translational accelerations at each point in world frame
            % Points only have translational accelerations - no angular accelerations
            
            n_points = size(obj.points_body, 2);
            accelerations_world = zeros(3, n_points);
            
            if n_points == 0
                return;
            end
            
            % Calculate accelerations at center first
            obj.linear_acceleration = obj.applied_force / obj.mass;
            obj.angular_acceleration = obj.inertia \ obj.applied_torque;
            
            % For each point, calculate the translational acceleration due to:
            % 1. Direct translation of center linear acceleration
            % 2. Translational acceleration due to angular acceleration via lever arm (α × r)
            % 3. Centripetal acceleration due to current angular velocity (ω × (ω × r))
            
            for i = 1:n_points
                lever_arm_body = obj.lever_arms(:, i);
                lever_arm_world = obj.R_body_to_world * lever_arm_body;
                
                % 1. Direct linear acceleration (same as center)
                accel_translational = obj.linear_acceleration;
                
                % 2. Translational acceleration due to angular acceleration (α × r)
                angular_accel_world = obj.angular_acceleration;
                accel_angular = cross(angular_accel_world, lever_arm_world);
                
                % 3. Centripetal acceleration due to current angular velocity (ω × (ω × r))
                angular_vel_world = obj.angular_velocity;
                accel_centripetal = cross(angular_vel_world, cross(angular_vel_world, lever_arm_world));
                
                % Total translational acceleration at point (in global coordinates)
                accelerations_world(:, i) = accel_translational + accel_angular + accel_centripetal;
            end
        end
        
        function points_world = getPointsWorldFrame(obj)
            % Get point positions in world frame
            points_world = zeros(3, size(obj.points_body, 2));
            for i = 1:size(obj.points_body, 2)
                points_world(:, i) = obj.center_pos + obj.R_body_to_world * obj.points_body(:, i);
            end
        end
        
        function setupVisualization(obj)
            % Setup 3D visualization
            obj.fig_handle = figure('Name', '6-DOF Plane Acceleration Simulation', 'Position', [100 100 1000 800]);
            obj.ax_handle = axes('Parent', obj.fig_handle);
            
            hold(obj.ax_handle, 'on');
            grid(obj.ax_handle, 'on');
            axis(obj.ax_handle, 'equal');
            xlabel(obj.ax_handle, 'X (m)');
            ylabel(obj.ax_handle, 'Y (m)');
            zlabel(obj.ax_handle, 'Z (m)');
            title(obj.ax_handle, '6-DOF Plane with Acceleration Calculation');
            view(obj.ax_handle, 45, 30);
        end
     
        function updateVisualization(obj, accel_scale)
            % Update the 3D visualization
            if nargin < 2
                accel_scale = 0.1; % Scale factor for acceleration vectors
            end

            % Clear previous plots
            delete(obj.plane_handle);
            delete(obj.points_handle);
            if ~isempty(obj.acceleration_handles)
                delete(obj.acceleration_handles);
            end

            % Draw plane (as a square in the XY plane of body frame)
            plane_size = 4.5;
            plane_corners_body = [
                -plane_size/2, plane_size/2, plane_size/2, -plane_size/2;
                -plane_size/2, -plane_size/2, plane_size/2, plane_size/2;
                0, 0, 0, 0
            ];

            % Transform to world frame
            plane_corners_world = zeros(size(plane_corners_body));
            for i = 1:4
                plane_corners_world(:, i) = obj.center_pos + obj.R_body_to_world * plane_corners_body(:, i);
            end

            % Plot plane
            obj.plane_handle = patch(obj.ax_handle, ...
                plane_corners_world(1, :), plane_corners_world(2, :), plane_corners_world(3, :), ...
                'blue', 'FaceAlpha', 0.3, 'EdgeColor', 'blue', 'LineWidth', 2);

            % Draw points and their accelerations
            if ~isempty(obj.points_body)
                points_world = obj.getPointsWorldFrame();
                accelerations_world = obj.calculatePointAccelerations();

                % Plot points
                obj.points_handle = scatter3(obj.ax_handle, ...
                    points_world(1, :), points_world(2, :), points_world(3, :), ...
                    100, 'red', 'filled', 'MarkerEdgeColor', 'black');

                % Plot acceleration vectors at points
                obj.acceleration_handles = [];
                for i = 1:size(points_world, 2)
                    % Translational acceleration vector
                    accel_vec = accelerations_world(:, i) * accel_scale;
                    h = quiver3(obj.ax_handle, ...
                        points_world(1, i), points_world(2, i), points_world(3, i), ...
                        accel_vec(1), accel_vec(2), accel_vec(3), ...
                        'magenta', 'LineWidth', 2, 'MaxHeadSize', 0.2);
                    obj.acceleration_handles = [obj.acceleration_handles, h];

                    % Connect point to plane center
                    z_plane = obj.center_pos(3);  
                    plot3(obj.ax_handle, ...
                          [points_world(1, i), points_world(1, i)], ...  % same x
                          [points_world(2, i), points_world(2, i)], ...  % same y
                          [points_world(3, i), z_plane], ...            % z from point down to plane
                          'k-', 'LineWidth', 1);
                end
            end

            % Draw applied force and torque at center ONLY
            if norm(obj.applied_force) > 0
                force_vec = obj.applied_force * accel_scale * 0.5;  % Scale differently for visualization
                quiver3(obj.ax_handle, ...
                    obj.center_pos(1), obj.center_pos(2), obj.center_pos(3), ...
                    force_vec(1), force_vec(2), force_vec(3), ...
                    'cyan', 'LineWidth', 3, 'MaxHeadSize', 0.15);
            end

            if norm(obj.applied_torque) > 0
                torque_vec = obj.applied_torque * accel_scale * 0.5;  % Scale differently for visualization
                quiver3(obj.ax_handle, ...
                    obj.center_pos(1), obj.center_pos(2), obj.center_pos(3), ...
                    torque_vec(1), torque_vec(2), torque_vec(3), ...
                    'yellow', 'LineWidth', 3, 'MaxHeadSize', 0.15);
            end

            % Set axis limits
            all_points = [obj.center_pos, obj.getPointsWorldFrame()];
            if ~isempty(all_points)
                xlim(obj.ax_handle, [-4,4]);
                ylim(obj.ax_handle, [-4,4]);
                zlim(obj.ax_handle, [-4,4]);
            end

            % Add legend
            legend(obj.ax_handle, ...
                'Location', 'best', 'FontSize', 8);
        end

        
        function printResults(obj)
            % Print current state and acceleration calculations
            fprintf('\n=== Plane Acceleration Simulation Results ===\n');
            fprintf('Plane Center Position: [%.2f, %.2f, %.2f] m\n', obj.center_pos);
            fprintf('Plane Orientation (Roll, Pitch, Yaw): [%.2f, %.2f, %.2f] rad\n', obj.orientation);
            fprintf('Applied Force at Center: [%.2f, %.2f, %.2f] N\n', obj.applied_force);
            fprintf('Applied Torque at Center: [%.2f, %.2f, %.2f] N⋅m\n', obj.applied_torque);
            fprintf('Linear Acceleration at Center: [%.2f, %.2f, %.2f] m/s²\n', obj.linear_acceleration);
            fprintf('Angular Acceleration at Center: [%.2f, %.2f, %.2f] rad/s²\n', obj.angular_acceleration);
            fprintf('Current Angular Velocity: [%.2f, %.2f, %.2f] rad/s\n', obj.angular_velocity);
            
            if ~isempty(obj.points_body)
                points_world = obj.getPointsWorldFrame();
                accelerations_world = obj.calculatePointAccelerations();
                
                fprintf('\nPoint Analysis:\n');
                for i = 1:size(obj.points_body, 2)
                    fprintf('Point %d:\n', i);
                    fprintf('  Body Frame Position: [%.2f, %.2f, %.2f] m\n', obj.points_body(:, i));
                    fprintf('  World Frame Position: [%.2f, %.2f, %.2f] m\n', points_world(:, i));
                    fprintf('  Lever Arm (Body): [%.2f, %.2f, %.2f] m\n', obj.lever_arms(:, i));
                    fprintf('  Translational Acceleration: [%.2f, %.2f, %.2f] m/s²\n', accelerations_world(:, i));
                    fprintf('  Acceleration Magnitude: %.2f m/s²\n\n', norm(accelerations_world(:, i)));
                end
            end
        end
        
        function runExample(obj)
            % Run an example simulation
            fprintf('Running example simulation...\n');
            
            % Set plane pose
            obj.setPlanePose([1; 2; 1], [0.2; 0.3; 0.1]);
            
            % Add some points above the plane
            obj.addPoint([1; 0.5; 0.5]);   % Point 1
            obj.addPoint([-0.8; -0.3; 0.8]); % Point 2
            obj.addPoint([0; 1.2; 0.3]);   % Point 3
            
            % Apply forces and torques at center only
            obj.setForceAndTorque([10; 5; -3], [2; -1; 4]);
            
            % Update visualization
            obj.updateVisualization(0.2);
            obj.printResults();
        end
        
        function simulateMotion(obj, dt, n_steps, force_traj, torque_traj)
            % Simulate motion over time with time-varying force/torque

            if nargin < 2, dt = 0.1; end
            if nargin < 3, n_steps = 50; end
            if nargin < 4, force_traj = [zeros(3,1), obj.applied_force]; end
            if nargin < 5, torque_traj = [zeros(3,1), obj.applied_torque]; end

            % Interpolate if only start/goal provided
            if size(force_traj,2) == 2
                force_traj = repmat(force_traj(:,1),1,n_steps) + ...
                    (repmat((0:(n_steps-1))/(n_steps-1),3,1) .* repmat(force_traj(:,2)-force_traj(:,1),1,n_steps));
            end
            if size(torque_traj,2) == 2
                torque_traj = repmat(torque_traj(:,1),1,n_steps) + ...
                    (repmat((0:(n_steps-1))/(n_steps-1),3,1) .* repmat(torque_traj(:,2)-torque_traj(:,1),1,n_steps));
            end

            fprintf('Simulating motion for %.1f seconds...\n', dt * n_steps);
            
            % Initialize velocities
            obj.velocity = [0; 0; 0];
            obj.angular_velocity = [0; 0; 0];
            
            % Trajectory storage
            pose_traj = zeros(6, n_steps); % [x;y;z;roll;pitch;yaw]
            accel_traj = zeros(6, n_steps); % [ax;ay;az;alpha_x;alpha_y;alpha_z]
            force_log = zeros(3, n_steps);
            torque_log = zeros(3, n_steps);
            point_accel_traj = zeros(3, size(obj.points_body,2), n_steps);

            for step = 1:n_steps
                % Use interpolated force and torque for each step
                obj.setForceAndTorque(force_traj(:,step), torque_traj(:,step));

                % Update velocities using current accelerations
                obj.velocity = obj.velocity + obj.linear_acceleration * dt;
                obj.angular_velocity = obj.angular_velocity + obj.angular_acceleration * dt;
                
                % Update position and orientation   
                obj.center_pos = obj.center_pos + obj.velocity * dt;
                obj.orientation = obj.orientation + obj.angular_velocity * dt;
                obj.updateRotationMatrix();
                
                % Store trajectory
                pose_traj(:, step) = [obj.center_pos; obj.orientation];
                accel_traj(:, step) = [obj.linear_acceleration; obj.angular_acceleration];
                force_log(:, step) = obj.applied_force;
                torque_log(:, step) = obj.applied_torque;
                accelerations_world = obj.calculatePointAccelerations();
                point_accel_traj(:,:,step) = accelerations_world;

                % Update visualization every few steps
                if mod(step, 5) == 0
                    obj.updateVisualization(0.1);
                    title(obj.ax_handle, sprintf('Step %d/%d - Time: %.1fs', step, n_steps, step*dt));
                    drawnow;
                    pause(0.05);
                end
            end

            obj.printResults();

            % Print initial, middle, and final values
            fprintf('\nForce at t=0:      [%.2f, %.2f, %.2f] N\n', force_log(:,1));
            fprintf('Force at t/2:     [%.2f, %.2f, %.2f] N\n', force_log(:,round(n_steps/2)));
            fprintf('Force at t_final: [%.2f, %.2f, %.2f] N\n', force_log(:,n_steps));
            
            fprintf('\nTorque at t=0:      [%.2f, %.2f, %.2f] N⋅m\n', torque_log(:,1));
            fprintf('Torque at t/2:     [%.2f, %.2f, %.2f] N⋅m\n', torque_log(:,round(n_steps/2)));
            fprintf('Torque at t_final: [%.2f, %.2f, %.2f] N⋅m\n', torque_log(:,n_steps));

            % --- Trajectory Plots ---
            t = (0:n_steps-1)*dt;
            figure('Name','Plane Pose Trajectory');
            subplot(2,1,1);
            plot(t, pose_traj(1,:), 'r', t, pose_traj(2,:), 'g', t, pose_traj(3,:), 'b');
            xlabel('Time (s)'); ylabel('Position (m)');
            legend('X','Y','Z'); title('COG Position Trajectory');
            grid on;
            subplot(2,1,2);
            plot(t, pose_traj(4,:), 'r', t, pose_traj(5,:), 'g', t, pose_traj(6,:), 'b');
            xlabel('Time (s)'); ylabel('Orientation (rad)');
            legend('Roll','Pitch','Yaw'); title('COG Orientation Trajectory');
            grid on;

            figure('Name','Translational Acceleration Trajectories');
            subplot(2,1,1);
            plot(t, accel_traj(1,:), 'r', t, accel_traj(2,:), 'g', t, accel_traj(3,:), 'b');
            xlabel('Time (s)'); ylabel('Linear Accel at COG (m/s²)');
            legend('ax','ay','az'); title('COG Linear Acceleration Trajectory');
            grid on;
            subplot(2,1,2);
            hold on;
            colors = lines(size(obj.points_body,2));
            for i = 1:size(obj.points_body,2)
                plot(t, squeeze(point_accel_traj(1,i,:)), '-', 'Color', colors(i,:), 'DisplayName',sprintf('Point %d ax',i));
                plot(t, squeeze(point_accel_traj(2,i,:)), '--', 'Color', colors(i,:), 'DisplayName',sprintf('Point %d ay',i));
                plot(t, squeeze(point_accel_traj(3,i,:)), ':', 'Color', colors(i,:), 'DisplayName',sprintf('Point %d az',i));
            end
            xlabel('Time (s)'); ylabel('Translational Accel at Points (m/s²)');
            legend('show'); title('Translational Accelerations at Points (Global Frame)');
            grid on;
            hold off;

            % --- Angular Acceleration at Center Only ---
            figure('Name','Angular Acceleration at Center');
            plot(t, accel_traj(4,:), 'r', t, accel_traj(5,:), 'g', t, accel_traj(6,:), 'b');
            xlabel('Time (s)'); ylabel('Angular Accel at COG (rad/s²)');
            legend('αx','αy','αz'); title('COG Angular Acceleration (Converted to Translational at Points)');
            grid on;
        end
    end
end

%% Example Usage Script
function runPlaneAccelerationExample()
    % Create and run example simulation
    sim = PlaneAccelerationSimulation();
    sim.runExample();
end

%% Utility Functions
function demonstrateSphericalCoordinates()
    % Demonstrate spherical coordinate system
    fprintf('\n=== Spherical Coordinates Demonstration ===\n');
    
    % Create points in spherical coordinates
    r = 3;              % radius
    theta = pi/4;       % azimuthal angle (0 to 2π)
    phi = pi/3;         % polar angle (0 to π)
    
    % Convert to Cartesian
    x = r * sin(phi) * cos(theta);
    y = r * sin(phi) * sin(theta);
    z = r * cos(phi);
    
    fprintf('Spherical: (r=%.2f, θ=%.2f, φ=%.2f)\n', r, theta, phi);
    fprintf('Cartesian: (x=%.2f, y=%.2f, z=%.2f)\n', x, y, z);
    
    % Visualize
    figure('Name', 'Spherical Coordinates');
    [X, Y, Z] = sphere(20);
    surf(X*r, Y*r, Z*r, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.2);
    hold on;
    
    % Plot point
    plot3(x, y, z, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red');
    
    % Plot coordinate lines
    plot3([0, x], [0, 0], [0, 0], 'r--', 'LineWidth', 2); % to projection on xy
    plot3([x, x], [0, y], [0, 0], 'g--', 'LineWidth', 2);
    plot3([x, x], [y, y], [0, z], 'b--', 'LineWidth', 2);
    plot3([0, x], [0, y], [0, z], 'k-', 'LineWidth', 3);   % radius vector
    
    axis equal;
    grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Spherical Coordinate System');
    legend('Unit Sphere', 'Point', 'X projection', 'Y projection', 'Z projection', 'Radius Vector');
end

%% Main execution
% Uncomment the line below to run the example
% runPlaneAccelerationExample();

% Uncomment the line below to see spherical coordinates demo
% demonstrateSphericalCoordinates();