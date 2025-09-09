function liquid_surface_animator(t, y, varargin)
% LIQUID_SURFACE_ANIMATOR - Animates 3D liquid surface with restart button
%
% Usage: liquid_surface_animator(t, y)
%        liquid_surface_animator(t, y, 'PropertyName', PropertyValue, ...)

%% Parse input arguments
p = inputParser;
addRequired(p, 't');
addRequired(p, 'y');
addParameter(p, 'AnimationSpeed', 1.0, @isnumeric);
addParameter(p, 'DiscRadius', 0.05, @isnumeric);
addParameter(p, 'SaveVideo', false, @islogical);
addParameter(p, 'VideoName', 'liquid_animation.mp4', @ischar);
addParameter(p, 'ShowTrails', true, @islogical);
addParameter(p, 'ContainerSize', [0.05, 0.05, 0.1], @isnumeric);
parse(p, t, y, varargin{:});

% Extract parameters
anim_speed = p.Results.AnimationSpeed;
disc_radius = p.Results.DiscRadius;
save_video = p.Results.SaveVideo;
video_name = p.Results.VideoName;
show_trails = p.Results.ShowTrails;
container_size = p.Results.ContainerSize;

%% Extract data from simulation
x_pos = y(:, 1);
y_pos = y(:, 2);
z_pos = y(:, 3);
phi   = y(:, 4); % roll
theta = y(:, 5); % pitch
psi   = y(:, 6); % yaw

%% Setup figure
fig = figure('Name', '3D Liquid Analysis: Movement vs Surface Orientation', ...
             'NumberTitle', 'off', 'Position', [50, 50, 1600, 800]);

paused = false;

%% UI
% Add restart button
uicontrol('Style', 'pushbutton', 'String', 'Restart Animation', ...
          'FontSize', 12, 'BackgroundColor', [0.8, 0.9, 1], ...
          'Units', 'normalized', 'Position', [0.3, 0.02, 0.12, 0.05], ...
          'Callback', @(src, event) run_animation());

uicontrol('Style', 'pushbutton', 'String', 'Pause', ...
          'FontSize', 12, 'BackgroundColor', [1 0.7 0.7], ...
          'Units', 'normalized', 'Position', [0.45, 0.02, 0.12, 0.05], ...
          'Callback', @(src,event) set_pause(true));

% Resume button
uicontrol('Style', 'pushbutton', 'String', 'Resume', ...
          'FontSize', 12, 'BackgroundColor', [0.7 1 0.7], ...
          'Units', 'normalized', 'Position', [0.6, 0.02, 0.12, 0.05], ...
          'Callback', @(src,event) set_pause(false));

function set_pause(state)
    paused = state;
end



%% LEFT PLOT: Movement
subplot(1, 2, 1);
ax_movement = gca;
hold on; grid on; axis equal;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Liquid Center Movement in Space'); view(45, 30);

max_pos = max([max(abs(x_pos)), max(abs(y_pos)), max(abs(z_pos))]);
axis_lim = max(0.15, max_pos * 1.2);
xlim([-axis_lim, axis_lim]); ylim([-axis_lim, axis_lim]);
zlim([min(z_pos)-0.05, max(z_pos)+0.1]);

[X_cyl, Y_cyl, Z_cyl] = cylinder(container_size(1), 20);
Z_cyl = Z_cyl * container_size(3);
surf(X_cyl, Y_cyl, Z_cyl, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.3, ...
     'FaceColor', [0.7, 0.7, 0.7]);

current_pos_marker = plot3(0, 0, 0, 'ro', 'MarkerSize', 12, ...
                          'MarkerFaceColor', 'r', 'LineWidth', 2);
if show_trails
    pos_trail = plot3(NaN, NaN, NaN, 'r-', 'LineWidth', 2, 'Color', [1, 0, 0, 0.6]);
end

%% RIGHT PLOT: Orientation
subplot(1, 2, 2);
ax_orientation = gca;
hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Liquid Surface Orientation'); view(315, 30);
xlim([-0.12, 0.12]); ylim([-0.12, 0.12]); zlim([-0.12, 0.12]);
%% cylinder
[Xc, Yc, Zc] = cylinder(container_size(1), 25);
Zc = Zc * container_size(3) - 0.02;
surf(Xc, Yc, Zc - container_size(3)/2, ... % center at z=0
     'FaceColor', [0.6, 0.6, 0.6], ...
     'FaceAlpha', 0.15, ...
     'EdgeColor', [0.3, 0.3, 0.3], ...
     'EdgeAlpha', 0.05);
%%
[disc_x, disc_y] = meshgrid(linspace(-disc_radius, disc_radius, 100)); % reduce mesh for performace
disc_mask = (disc_x.^2 + disc_y.^2) <= disc_radius^2;
disc_z = zeros(size(disc_x));
disc_x(~disc_mask) = NaN; disc_y(~disc_mask) = NaN; disc_z(~disc_mask) = NaN;
liquid_surf_orient = surf(disc_x, disc_y, disc_z, 'FaceColor', [0, 0.5, 1], ...
                          'FaceAlpha', 0.7, 'EdgeColor', 'none');

vector_length = 0.08;
normal_line_orient   = plot3([0 0], [0 0], [0 vector_length], 'r-', 'LineWidth', 3);
tangent1_line_orient = plot3([0 vector_length], [0 0], [0 0], 'g-', 'LineWidth', 3);
tangent2_line_orient = plot3([0 0], [0 vector_length], [0 0], 'b-', 'LineWidth', 3);

if show_trails
    trail_normal_orient   = plot3(NaN, NaN, NaN, 'r:', 'LineWidth', 1.5);
    trail_tangent1_orient = plot3(NaN, NaN, NaN, 'g:', 'LineWidth', 1.5);
    trail_tangent2_orient = plot3(NaN, NaN, NaN, 'b:', 'LineWidth', 1.5);
end

%% Text displays
pos_text   = annotation('textbox', [0.02, 0.85, 0.2, 0.12], 'String', '', ...
                        'BackgroundColor', 'white', 'EdgeColor', 'black', ...
                        'FontSize', 10, 'FontWeight', 'bold');
angle_text = annotation('textbox', [0.52, 0.85, 0.2, 0.12], 'String', '', ...
                        'BackgroundColor', 'white', 'EdgeColor', 'black', ...
                        'FontSize', 10, 'FontWeight', 'bold');

%% Animation settings
dt = mean(diff(t));
frame_skip = max(1, round(0.05 / (dt * anim_speed)));
n_frames = length(t);

% Trail storage
if show_trails
    trail_length = min(300, n_frames);
    pos_trail_data         = NaN(trail_length, 3);
    orient_trail_normal    = NaN(trail_length, 3);
    orient_trail_tangent1  = NaN(trail_length, 3);
    orient_trail_tangent2  = NaN(trail_length, 3);
end

%% === Nested Function: Animation Loop ===
    function run_animation()

        % Reset trails
        if show_trails
            pos_trail_data(:)        = NaN;
            orient_trail_normal(:)   = NaN;
            orient_trail_tangent1(:) = NaN;
            orient_trail_tangent2(:) = NaN;
        end
        
        fprintf('Starting animation with %d frames...\n', ceil(n_frames/frame_skip));
        
        for i = 1:frame_skip:n_frames
            while paused
                pause(0.05);
            end
            current_time  = t(i);
            current_pos   = [x_pos(i), y_pos(i), z_pos(i)];
            current_phi   = phi(i);
            current_theta = theta(i);
            current_psi   = psi(i);

            % Update marker
            set(current_pos_marker, 'XData', current_pos(1), ...
                                    'YData', current_pos(2), ...
                                    'ZData', current_pos(3));
            
            % Update trail
            if show_trails
                pos_trail_data = [current_pos; pos_trail_data(1:end-1, :)];
                valid_idx = ~isnan(pos_trail_data(:,1));
                set(pos_trail, 'XData', pos_trail_data(valid_idx,1), ...
                               'YData', pos_trail_data(valid_idx,2), ...
                               'ZData', pos_trail_data(valid_idx,3));
            end
            
            % Orientation rotation matrix
            R_x = [1 0 0; 0 cos(current_phi) -sin(current_phi); 0 sin(current_phi) cos(current_phi)];
            R_y = [cos(current_theta) 0 sin(current_theta); 0 1 0; -sin(current_theta) 0 cos(current_theta)];
            R_z = [cos(current_psi) -sin(current_psi) 0; sin(current_psi) cos(current_psi) 0; 0 0 1];
            R_total = R_z * R_y * R_x;

            % Rotate disc
            disc_points = [disc_x(:), disc_y(:), disc_z(:)]';
            disc_rot = R_total * disc_points;
            set(liquid_surf_orient, 'XData', reshape(disc_rot(1,:), size(disc_x)), ...
                                    'YData', reshape(disc_rot(2,:), size(disc_y)), ...
                                    'ZData', reshape(disc_rot(3,:), size(disc_z)));
            
            % Vectors
            normal_rot   = R_total * [0;0;vector_length];
            tangent1_rot = R_total * [vector_length;0;0];
            tangent2_rot = R_total * [0;vector_length;0];
            set(normal_line_orient,   'XData', [0 normal_rot(1)],   'YData', [0 normal_rot(2)],   'ZData', [0 normal_rot(3)]);
            set(tangent1_line_orient, 'XData', [0 tangent1_rot(1)], 'YData', [0 tangent1_rot(2)], 'ZData', [0 tangent1_rot(3)]);
            set(tangent2_line_orient, 'XData', [0 tangent2_rot(1)], 'YData', [0 tangent2_rot(2)], 'ZData', [0 tangent2_rot(3)]);
            
            % Trails for orientation
            if show_trails
                orient_trail_normal   = [normal_rot';   orient_trail_normal(1:end-1,:)];
                orient_trail_tangent1 = [tangent1_rot'; orient_trail_tangent1(1:end-1,:)];
                orient_trail_tangent2 = [tangent2_rot'; orient_trail_tangent2(1:end-1,:)];
                set(trail_normal_orient,   'XData', orient_trail_normal(:,1), 'YData', orient_trail_normal(:,2), 'ZData', orient_trail_normal(:,3));
                set(trail_tangent1_orient, 'XData', orient_trail_tangent1(:,1), 'YData', orient_trail_tangent1(:,2), 'ZData', orient_trail_tangent1(:,3));
                set(trail_tangent2_orient, 'XData', orient_trail_tangent2(:,1), 'YData', orient_trail_tangent2(:,2), 'ZData', orient_trail_tangent2(:,3));
            end
            
            % Text update
            set(pos_text, 'String', sprintf('POSITION (t=%.2fs)\nX: %+.4f\nY: %+.4f\nZ: %+.4f', ...
                                             current_time, current_pos));
            set(angle_text, 'String', sprintf('ORIENTATION (t=%.2fs)\nφ=%.2f° θ=%.2f° ψ=%.2f°', ...
                                              current_time, current_phi*180/pi, current_theta*180/pi, current_psi*180/pi));
            
            drawnow;
            pause(dt/anim_speed);
        end
        fprintf('Animation finished.\n');
    end % end run_animation

%% Run once at start
run_animation();

end % main function
