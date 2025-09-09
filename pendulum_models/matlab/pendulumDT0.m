function xk1 = pendulumDT0(xk, uk, Ts, params)
%% Discrete-time nonlinear dynamic model of a pendulum on a cart at time k
%
% 4 states (xk): 
%   cart position (z)
%   cart velocity (z_dot): when positive, cart moves to right
%   angle (theta): when 0, pendulum is at upright position
%   angular velocity (theta_dot): anticlockwise positive
% 
% 1 inputs: (uk)
%   force (F): when positive, force pushes cart to right 
%
% xk1 is the states at time k+1.
%
% Copyright 2018 The MathWorks, Inc.

%#codegen



% M = 1;      % cart mass
% m = 1;      % pendulum mass
% g = 9.81;   % gravity of earth
% l = 0.5;    % pendulum length
% Kd = 10;    % cart damping
% 
% params(1) = 1;      % cart mass
% params(2) = 1;      % pendulum mass
% params(3) = 9.81;   % gravity of earth
% params(4) = 0.5;    % pendulum length
% params(5) = 10;    % cart damping

% Repeat application of Euler method sampled at Ts/Nd.
Nd = 10;
delta = Ts/Nd;
xk1 = xk;
for ct=1:Nd
    xk1 = xk1 + delta*pendulumCT0(xk1,uk, params);
end
% Note that we choose the Euler method (first order Runge-Kutta method)
% because it is more efficient for plant with non-stiff ODEs.  You can
% choose other ODE solvers such as ode23, ode45 for better accuracy or
% ode15s and ode23s for stiff ODEs.  Those solvers are available from
% MATLAB.
