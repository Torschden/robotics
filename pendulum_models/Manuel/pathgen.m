function traj = pathgen(p0, p1, p2, t0, t2)
% cubicSpline3Points Generate a cubic spline trajectory through 3 points
% 
% Inputs:
%   p0 : [3x1] start point (x,y,z)
%   p1 : [3x1] intermediate point
%   p2 : [3x1] end point
%   t0 : start time
%   t2 : end time
%
% Output:
%   traj : struct with function handles
%          traj.pos(t) -> position [3x1]
%          traj.vel(t) -> velocity [3x1]
%          traj.acc(t) -> acceleration [3x1]

    % Force intermediate time at midpoint
    t1 = (t0 + t2)/2;

    % Build splines for each coordinate
    pp_x = spline([t0 t1 t2], [p0(1) p1(1) p2(1)]);
    pp_y = spline([t0 t1 t2], [p0(2) p1(2) p2(2)]);
    pp_z = spline([t0 t1 t2], [p0(3) p1(3) p2(3)]);

    % Derivatives
    pp_vx = fnder(pp_x,1);  pp_ax = fnder(pp_x,2);
    pp_vy = fnder(pp_y,1);  pp_ay = fnder(pp_y,2);
    pp_vz = fnder(pp_z,1);  pp_az = fnder(pp_z,2);

    % Function handles
    traj.pos = @(t) [ppval(pp_x,t); ppval(pp_y,t); ppval(pp_z,t)];
    traj.vel = @(t) [ppval(pp_vx,t); ppval(pp_vy,t); ppval(pp_vz,t)];
    traj.acc = @(t) [ppval(pp_ax,t); ppval(pp_ay,t); ppval(pp_az,t)];
end
