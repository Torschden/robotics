function main1
% regular, hanging pendulum

    function dxdt = pendulumCT0(x, u)
        % Continuous-time nonlinear dynamic model of liquid sloshing in a vessel
        % Modeled as an equivalent pendulum system with sloshing characteristics
        
        % Physical parameters for liquid sloshing
        M = 1.0;      % vessel mass (kg)
        m_s = 0.5;    % equivalent sloshing mass (kg) - fraction of total liquid
        m_ns = 0.5;   % non-sloshing (rigid) liquid mass (kg)
        g = 9.81;     % gravity (m/s^2)
        l = 0.3;      % equivalent pendulum length (m) - related to liquid height
        
        % Damping coefficients for liquid sloshing
        Kd_vessel = 8;    % vessel damping (N⋅s/m)
        zeta_s = 0.05;    % sloshing damping ratio (typical for water: 0.01-0.1)
        
        % Natural sloshing frequency (rad/s)
        omega_n = sqrt(g/l);
        c_s = 2*zeta_s*omega_n*m_s*l;  % sloshing damping coefficient
        
        % State variables
        z_dot = x(2);      % vessel velocity
        theta = x(3);      % sloshing angle
        theta_dot = x(4);  % sloshing angular velocity
        F = u;             % applied force
        
        % Additional sloshing effects
        % Horizontal acceleration of vessel affects sloshing dynamics
        z_ddot_vessel = (F - Kd_vessel*z_dot - m_s*l*theta_dot^2*sin(theta) + m_s*g*sin(theta)*cos(theta)) / ...
                       (M + m_ns + m_s*sin(theta)^2);
        
        % Sloshing dynamics with vessel acceleration coupling
        % Enhanced model includes:
        % 1. Gravitational restoring force
        % 2. Coupling with vessel acceleration
        % 3. Sloshing-specific damping
        % 4. Nonlinear sloshing effects (for large amplitudes)
        
        theta_ddot = (1/l) * (g*sin(theta) + z_ddot_vessel*cos(theta)) - ...
                     (c_s/(m_s*l^2))*theta_dot - ...
                     0.1*theta_dot*abs(theta_dot); % Additional nonlinear damping for large sloshing
        
        % Vessel dynamics including sloshing reaction forces
        z_ddot = (F - Kd_vessel*z_dot - m_s*l*(theta_ddot*sin(theta) + theta_dot^2*cos(theta))) / ...
                 (M + m_ns + m_s);
        
        dxdt = [
            z_dot;           % vessel velocity
            z_ddot;          % vessel acceleration  
            theta_dot;       % sloshing angular velocity
            theta_ddot;      % sloshing angular acceleration
        ];
    end

% --- pendulumDT0 ---
    function xk1 = pendulumDT0(xk, uk, Ts)
        Nd = 10;
        delta = Ts/Nd;
        xk1 = xk;
        for ct=1:Nd
            xk1 = xk1 + delta*pendulumCT0(xk1,uk);
        end
    end

% --- pendulumStateFcn ---
    function xk1 = pendulumStateFcn(xk, u)
        uk = u(1);
        Ts = u(2);
        xk1 = pendulumDT0(xk, uk, Ts);
    end

% --- pendulumOutputFcn ---
    function y = pendulumOutputFcn(x,u,params)
        y = [x(1); x(3)];
    end

% --- pendulumMeasurementFcn ---
    function y = pendulumMeasurementFcn(xk)
        y = [xk(1);xk(3)];
    end

% --- Main simulation code ---
    nx = 4;
    ny = 2;
    nu = 1;
    nlobj = nlmpc(nx, ny, nu);


    Ts = 0.05;
    nlobj.Ts = Ts;
    nlobj.PredictionHorizon = 16;
    nlobj.ControlHorizon = 5;

    nlobj.Model.StateFcn = @(xk,u,Ts) pendulumDT0(xk,u,Ts);
    nlobj.Model.IsContinuousTime = false;
    nlobj.Model.NumberOfParameters = 1;
    nlobj.Model.OutputFcn = @(x,u,params) pendulumOutputFcn(x,u,params);

    nlobj.Jacobian.OutputFcn = @(x,u,Ts) [1 0 0 0; 0 0 1 0];

    nlobj.Weights.OutputVariables = [5 10];
    nlobj.Weights.ManipulatedVariablesRate = 0.1;

    nlobj.OV(1).Min = -10;
    nlobj.OV(1).Max = 10;

    nlobj.MV.Min = -10;
    nlobj.MV.Max = 10;

    % Initial state: pendulum near downward position (theta ~ 0)
    x0 = [0.1;0.2;0.1;0.3];
    u0 = 0.4;
    validateFcns(nlobj,x0,u0,[],{Ts});

    EKF = extendedKalmanFilter(@pendulumStateFcn, @pendulumMeasurementFcn);

    % Start at downward position
    x = [0;0;0;0];
    y = [x(1);x(3)];
    EKF.State = x;

    mv = 0;

    % Reference: cart at 0, pendulum angle at 0 (downward)
    yref1 = [0 0];
    yref2 = [5 0];

    nloptions = nlmpcmoveopt;
    nloptions.Parameters = {Ts};

    Duration = 50;
    nSteps = Duration/Ts;
    hbar = waitbar(0,'Simulation Progress');
    xHistory = x;
    uHistory = mv; % Store input history
    
    for ct = 1:nSteps
        if ct*Ts<10
            yref = yref1;
        else
            yref = yref2;
        end
        xk = correct(EKF, y);
        [mv,nloptions,info] = nlmpcmove(nlobj,xk,mv,yref,[],nloptions);
        predict(EKF, [mv; Ts]);
        x = pendulumDT0(x,mv,Ts);
        % y = x([1 3]) + randn(2,1)*0.01; 
        y = x([1 3]); 
        xHistory = [xHistory x]; %#ok<*AGROW>
        uHistory = [uHistory mv]; %#ok<*AGROW>
        
        % Update progress bar with proper progress calculation
        progress = ct / nSteps;
        waitbar(progress, hbar, sprintf('Simulation Progress: %d/%d steps (%.1f%%)', ct, nSteps, progress*100));
    end
    close(hbar)

    % Create time vector
    timeVec = 0:Ts:Duration;

    % Plot results including input
    figure('Position', [100, 100, 1200, 800])
    subplot(2,3,1)
    plot(timeVec, xHistory(1,:), 'b-', 'LineWidth', 2)
    xlabel('Time (s)')
    ylabel('Vessel Position (m)')
    title('Liquid Vessel Position')
    grid on

    subplot(2,3,2)
    plot(timeVec, xHistory(2,:), 'r-', 'LineWidth', 2)
    xlabel('Time (s)')
    ylabel('Vessel Velocity (m/s)')
    title('Liquid Vessel Velocity')
    grid on

    subplot(2,3,3)
    plot(timeVec, xHistory(3,:)*180/pi, 'g-', 'LineWidth', 2)
    xlabel('Time (s)')
    ylabel('Sloshing Angle (deg)')
    title('Liquid Sloshing Angle')
    grid on

    subplot(2,3,4)
    plot(timeVec, xHistory(4,:)*180/pi, 'm-', 'LineWidth', 2)
    xlabel('Time (s)')
    ylabel('Sloshing Angular Velocity (deg/s)')
    title('Liquid Sloshing Rate')
    grid on

    subplot(2,3,5)
    plot(timeVec, uHistory, 'k-', 'LineWidth', 2)
    xlabel('Time (s)')
    ylabel('Applied Force (N)')
    title('Control Force on Vessel')
    grid on
    
    % Add reference changes as vertical lines
    hold on
    line([10, 10], [min(uHistory), max(uHistory)], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2)
    text(10.5, max(uHistory)*0.8, 'Reference Change', 'FontSize', 10)
    hold off

    % Create animation
    animFig = figure('Position', [200, 200, 900, 600]);
    l = 0.4; % equivalent pendulum length for sloshing
    vessel_width = 0.8;
    vessel_height = 0.6;
    liquid_height = 0.4;
    
    % Set up the plot limits
    xlim_range = [min(xHistory(1,:))-2, max(xHistory(1,:))+2];
    ylim_range = [-0.3, vessel_height+0.5];
    
    % Add animation controls
    uicontrol('Style', 'pushbutton', 'String', 'Play Animation', ...
              'Position', [20, 20, 100, 40], ...
              'Callback', @(~,~) playAnimation());
    
    % Animation function
    function playAnimation()
        figure(animFig);
        
        for i = 1:5:length(timeVec) % Animate every 5th frame
            if ~ishandle(animFig), break; end % Stop if figure is closed
            
            figure(animFig);
            clf
            
            % Re-add the play button
            uicontrol('Style', 'pushbutton', 'String', 'Play Animation', ...
                      'Position', [20, 20, 100, 40], ...
                      'Callback', @(~,~) playAnimation());
            
            hold on
            
            % Current state
            z = xHistory(1,i);
            theta = xHistory(3,i);
            
            % Draw vessel container (rectangular tank)
            vessel_x = [z-vessel_width/2, z+vessel_width/2, z+vessel_width/2, z-vessel_width/2, z-vessel_width/2];
            vessel_y = [0, 0, vessel_height, vessel_height, 0];
            plot(vessel_x, vessel_y, 'k-', 'LineWidth', 3) % Vessel outline
            
            % Draw liquid surface (sloshing)
            % Create a sinusoidal surface based on sloshing angle
            n_points = 20;
            x_surface = linspace(z-vessel_width/2+0.05, z+vessel_width/2-0.05, n_points);
            sloshing_amplitude = 0.1 * sin(theta); % Surface wave amplitude
            y_base = liquid_height;
            y_surface = y_base + sloshing_amplitude * sin(2*pi*(x_surface-z)/(vessel_width/2));
            
            % Fill liquid area
            liquid_x = [z-vessel_width/2+0.05, x_surface, z+vessel_width/2-0.05];
            liquid_y = [0.05, y_surface, 0.05];
            fill(liquid_x, liquid_y, [0.3 0.6 1], 'EdgeColor', [0.2 0.4 0.8], 'LineWidth', 2, 'FaceAlpha', 0.7)
            
            % Draw sloshing motion indicator (equivalent pendulum)
            sloshing_x = z + l*sin(theta);
            sloshing_y = vessel_height/2 - l*cos(theta);
            plot([z, sloshing_x], [vessel_height/2, sloshing_y], 'r--', 'LineWidth', 2)
            plot(sloshing_x, sloshing_y, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')
            text(sloshing_x+0.1, sloshing_y, 'Sloshing Center', 'FontSize', 8)
            
            % Draw center of mass of vessel
            plot(z, vessel_height/2, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
            
            % Draw ground/platform
            plot(xlim_range, [-0.1, -0.1], 'k-', 'LineWidth', 4)
            
            % Add force arrow if significant
            if abs(uHistory(i)) > 1
                force_scale = 0.015;
                arrow_length = uHistory(i) * force_scale;
                arrow_y = vessel_height/2;
                if arrow_length > 0
                    arrow([z, arrow_y], [z + arrow_length, arrow_y], 'Color', 'g', 'LineWidth', 3)
                    text(z + arrow_length/2, arrow_y + 0.15, sprintf('F=%.1fN', uHistory(i)), ...
                         'HorizontalAlignment', 'center', 'FontSize', 10)
                else
                    arrow([z, arrow_y], [z + arrow_length, arrow_y], 'Color', 'g', 'LineWidth', 3)
                    text(z + arrow_length/2, arrow_y + 0.15, sprintf('F=%.1fN', uHistory(i)), ...
                         'HorizontalAlignment', 'center', 'FontSize', 10)
                end
            end
            
            xlim(xlim_range)
            ylim(ylim_range)
            xlabel('Position (m)')
            ylabel('Height (m)')
            title(sprintf('Liquid Sloshing in Moving Vessel - Time: %.1f s', timeVec(i)))
            grid on
            
            % Add text showing current state
            text(xlim_range(1)+0.5, ylim_range(2)-0.1, sprintf('Vessel Pos: %.2f m', z), ...
                 'FontSize', 12, 'BackgroundColor', 'white')
            text(xlim_range(1)+0.5, ylim_range(2)-0.2, sprintf('Sloshing Angle: %.1f°', theta*180/pi), ...
                 'FontSize', 12, 'BackgroundColor', 'white')
            text(xlim_range(1)+0.5, ylim_range(2)-0.3, sprintf('Applied Force: %.1f N', uHistory(i)), ...
                 'FontSize', 12, 'BackgroundColor', 'white')
            
            % Show reference target
            if timeVec(i) < 10
                ref_text = 'Target: Vessel at 0m, Liquid stable';
            else
                ref_text = 'Target: Vessel at 5m, Liquid stable';
            end
            text(xlim_range(1)+0.5, ylim_range(2)-0.4, ref_text, ...
                 'FontSize', 10, 'BackgroundColor', 'yellow')
            
            % Add sloshing frequency info
            omega_n = sqrt(9.81/l);
            text(xlim_range(2)-3, ylim_range(2)-0.1, sprintf('Natural Freq: %.2f Hz', omega_n/(2*pi)), ...
                 'FontSize', 10, 'BackgroundColor', 'cyan')
            
            drawnow
            pause(0.08) % Adjust for animation speed
        end
        hold off
    end
    
    % Play animation once automatically
    playAnimation();
end

% Helper function to draw arrows (if not available in your MATLAB version)
function arrow(start_point, end_point, varargin)
    % Simple arrow drawing function
    x = [start_point(1), end_point(1)];
    y = [start_point(2), end_point(2)];
    
    % Draw the line
    plot(x, y, varargin{:})
    
    % Draw arrowhead
    dx = end_point(1) - start_point(1);
    dy = end_point(2) - start_point(2);
    
    if abs(dx) > 1e-6 || abs(dy) > 1e-6
        length_arrow = sqrt(dx^2 + dy^2);
        unit_x = dx / length_arrow;
        unit_y = dy / length_arrow;
        
        % Arrowhead size
        head_length = 0.1;
        head_width = 0.05;
        
        % Arrowhead points
        head_x = [end_point(1) - head_length*unit_x + head_width*unit_y, ...
                  end_point(1), ...
                  end_point(1) - head_length*unit_x - head_width*unit_y];
        head_y = [end_point(2) - head_length*unit_y - head_width*unit_x, ...
                  end_point(2), ...
                  end_point(2) - head_length*unit_y + head_width*unit_x];
        
        fill(head_x, head_y, 'g', 'EdgeColor', 'g')
    end
end