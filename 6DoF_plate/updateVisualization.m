function updateVisualization(obj, force_scale)
    % Update the 3D visualization
    if nargin < 2
        force_scale = 0.1; % Scale factor for force vectors
    end
    
    % Clear previous plots
    delete(obj.plane_handle);
    delete(obj.points_handle);
    if ~isempty(obj.force_handles)
        delete(obj.force_handles);
    end
    
    % Draw plane (as a square in the XY plane of body frame)
    plane_size = 2;
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
    
    % Draw coordinate system at plane center
    coord_length = 1;
    x_axis_world = obj.center_pos + obj.R_body_to_world * [coord_length; 0; 0];
    y_axis_world = obj.center_pos + obj.R_body_to_world * [0; coord_length; 0];
    z_axis_world = obj.center_pos + obj.R_body_to_world * [0; 0; coord_length];
    
    quiver3(obj.ax_handle, obj.center_pos(1), obj.center_pos(2), obj.center_pos(3), ...
        x_axis_world(1)-obj.center_pos(1), x_axis_world(2)-obj.center_pos(2), x_axis_world(3)-obj.center_pos(3), ...
        'r', 'LineWidth', 3, 'MaxHeadSize', 0.1);
    quiver3(obj.ax_handle, obj.center_pos(1), obj.center_pos(2), obj.center_pos(3), ...
        y_axis_world(1)-obj.center_pos(1), y_axis_world(2)-obj.center_pos(2), y_axis_world(3)-obj.center_pos(3), ...
        'g', 'LineWidth', 3, 'MaxHeadSize', 0.1);
    quiver3(obj.ax_handle, obj.center_pos(1), obj.center_pos(2), obj.center_pos(3), ...
        z_axis_world(1)-obj.center_pos(1), z_axis_world(2)-obj.center_pos(2), z_axis_world(3)-obj.center_pos(3), ...
        'b', 'LineWidth', 3, 'MaxHeadSize', 0.1);
    
    % Draw points and their forces
    if ~isempty(obj.points_body)
        points_world = obj.getPointsWorldFrame();
        [forces_world, ~] = obj.calculatePointForces();
        
        % Plot points
        obj.points_handle = scatter3(obj.ax_handle, ...
            points_world(1, :), points_world(2, :), points_world(3, :), ...
            100, 'red', 'filled', 'MarkerEdgeColor', 'black');
        
        % Plot force vectors at points
        obj.force_handles = [];
        for i = 1:size(points_world, 2)
            % Force vector
            force_vec = forces_world(:, i) * force_scale;
            h = quiver3(obj.ax_handle, ...
                points_world(1, i), points_world(2, i), points_world(3, i), ...
                force_vec(1), force_vec(2), force_vec(3), ...
                'magenta', 'LineWidth', 2, 'MaxHeadSize', 0.2);
            obj.force_handles = [obj.force_handles, h];
            
            % Connect point to plane center
            plot3(obj.ax_handle, ...
                [obj.center_pos(1), points_world(1, i)], ...
                [obj.center_pos(2), points_world(2, i)], ...
                [obj.center_pos(3), points_world(3, i)], ...
                'k--', 'LineWidth', 1);
        end
    end
    
    % Draw applied force and torque at center
    if norm(obj.applied_force) > 0
        force_vec = obj.applied_force * force_scale;
        quiver3(obj.ax_handle, ...
            obj.center_pos(1), obj.center_pos(2), obj.center_pos(3), ...
            force_vec(1), force_vec(2), force_vec(3), ...
            'cyan', 'LineWidth', 3, 'MaxHeadSize', 0.15);
    end
    
    if norm(obj.applied_torque) > 0
        torque_vec = obj.applied_torque * force_scale;
        quiver3(obj.ax_handle, ...
            obj.center_pos(1), obj.center_pos(2), obj.center_pos(3), ...
            torque_vec(1), torque_vec(2), torque_vec(3), ...
            'yellow', 'LineWidth', 3, 'MaxHeadSize', 0.15);
    end
    
    % Set axis limits
    all_points = [obj.center_pos, obj.getPointsWorldFrame()];
    if ~isempty(all_points)
        % margin = 3;
        % xlim(obj.ax_handle, [min(all_points(1, :)) - margin, max(all_points(1, :)) + margin]);
        % ylim(obj.ax_handle, [min(all_points(2, :)) - margin, max(all_points(2, :)) + margin]);
        % zlim(obj.ax_handle, [min(all_points(3, :)) - margin, max(all_points(3, :)) + margin]);
        % constant 
        xlim(obj.ax_handle, [-4,4]);
        ylim(obj.ax_handle, [-4,4]);
        zlim(obj.ax_handle, [-4,4]);
    end
    
    % Add legend
    legend(obj.ax_handle, {'Plane', 'X-axis', 'Y-axis', 'Z-axis', 'Points', 'Point Forces', ...
        'Lever Arms', 'Applied Force', 'Applied Torque'}, ...
        'Location', 'best', 'FontSize', 8);
end
