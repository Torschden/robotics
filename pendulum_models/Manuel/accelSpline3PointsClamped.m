function [u_pp, traj] = accelSpline3PointsClamped(p0, p1, p2, t0, t2, v0, v2)
% accelSpline3PointsClamped  Per-axis acceleration handle for a 3-point path
% using clamped cubic splines (endpoint velocities specified).
%
% Inputs:
%   p0,p1,p2 : 3x1 position vectors (x;y;z) at times t0, t1, t2 respectively
%              where t1 = (t0+t2)/2 (you can change t1 if needed)
%   t0,t2    : scalars, start and end times (t2 > t0)
%   v0,v2    : 3x1 endpoint velocities at t0 and t2 (optional; default zeros)
%
% Outputs:
%   u_pp(t)  : function handle returning [ax; ay; az] at time t (3xN)
%   traj     : struct with pp-forms for pos/vel/acc per axis and breaks

    if t2 <= t0, error('t2 must be greater than t0'); end
    t1 = 0.5*(t0 + t2);

    p0 = p0(:); p1 = p1(:); p2 = p2(:);
    if nargin < 6 || isempty(v0), v0 = [0;0;0]; else, v0 = v0(:); end
    if nargin < 7 || isempty(v2), v2 = [0;0;0]; else, v2 = v2(:); end

    tx = [t0 t1 t2];

    % --- Build clamped splines by augmenting with endpoint slopes (velocities)
    pp_x = spline(tx, [v0(1) p0(1) p1(1) p2(1) v2(1)]);
    pp_y = spline(tx, [v0(2) p0(2) p1(2) p2(2) v2(2)]);
    pp_z = spline(tx, [v0(3) p0(3) p1(3) p2(3) v2(3)]);

    % --- Derivatives
    apx = fnder(pp_x, 2);
    apy = fnder(pp_y, 2);
    apz = fnder(pp_z, 2);

    % --- Gate to zero outside [t0,t2] to mimic step-like support
    u_pp = @(t) accel_pp_gated(t, apx, apy, apz, t0, t2);

    % --- Return pp forms (optional, useful for debugging/plotting)
    traj.pos.x = pp_x; traj.pos.y = pp_y; traj.pos.z = pp_z;
    traj.acc.x = apx;  traj.acc.y = apy;  traj.acc.z = apz;
    traj.breaks = apx.breaks;   % same breaks for all axes here

end

% ===== helper =====
function U = accel_pp_gated(t, apx, apy, apz, t0, t2)
    t = t(:).';                       % row vector
    ax = ppval(apx, t);
    ay = ppval(apy, t);
    az = ppval(apz, t);
    mask = (t >= t0) & (t <= t2);     % gate outside support
    U = [ax; ay; az] .* [mask; mask; mask];
end