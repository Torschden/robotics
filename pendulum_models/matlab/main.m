clear 
clc

% sim parameters
simTime = 20;
Ts = 0.1;

params.M = 10;     % cart mass
params.m = 1;      % pendulum mass
params.g = 9.81;   % gravity of earth
params.l = 0.5;    % pendulum length
params.Kd = 10;    % cart damping
params.Kp = 1;   % pendulum damping

% limits
cartSpeedLimits = [-10 10];
cartForceLimits = [-500 500];

% MPC parameters
weightGoal = 3;
weightSlosh = 7;
rateLimitation = 0.05;
nStates = 4;
nOutputs = 2;
nInputs = 1;

% trajectory
yref1 = [0 0];
yref2 = [5 0];
yref3 = [0 0];

nlobj = nlmpc(nStates, nOutputs, nInputs);
nlobj.Ts = Ts;
nlobj.PredictionHorizon = 15;
nlobj.ControlHorizon = 7;
nlobj.Model.StateFcn = @(x,u,Ts) pendulumDT0(x,u,Ts,params);
nlobj.Model.IsContinuousTime = false;
nlobj.Model.NumberOfParameters = 1;
nlobj.Model.OutputFcn = 'pendulumOutputFcn';
nlobj.Jacobian.OutputFcn = @(x,u,Ts) [1 0 0 0; 0 0 1 0];
nlobj.Weights.OutputVariables = [weightGoal weightSlosh];
nlobj.Weights.ManipulatedVariablesRate = rateLimitation;
nlobj.OV(1).Min = cartSpeedLimits(1);
nlobj.OV(1).Max = cartSpeedLimits(2);
nlobj.MV.Min = cartForceLimits(1);
nlobj.MV.Max = cartForceLimits(2);

x0 = [0.1;0.2;0.1;0.3];
u0 = 0.4;
validateFcns(nlobj,x0,u0,[],{Ts});

EKF = extendedKalmanFilter(@(x,u) pendulumDT0(x,u,Ts,params), @pendulumMeasurementFcn);

x = [0;0;0;0];
y = [x(1);x(3)];
EKF.State = x;

mv = 0;

nloptions = nlmpcmoveopt;
nloptions.Parameters = {Ts};

hbar = waitbar(0,'Simulation Progress');
sol = zeros(length(x)+1,simTime/Ts+1);
sol(:,1) = [x; mv];
% sim loop
for ct = 1:(20/Ts)
    if ct*Ts<7
        yref = yref1;
    elseif ct*Ts<15
        yref = yref2;
    else
        yref = yref3;
    end
    xk = correct(EKF, y);
    [mv,nloptions,info] = nlmpcmove(nlobj,xk,mv,yref,[],nloptions);
    predict(EKF, mv);
    x = [pendulumDT0(x(1:4),mv,Ts,params); mv];
    % y = x([1 3]) + randn(2,1)*0.01; % white noise
    y = x([1 3]);
    sol(:,ct+1) = x;
    waitbar(ct*Ts/20,hbar);
end
close(hbar)
% plot
plotting(Ts,simTime,sol)
disp(max(abs(sol(3,:)/pi*180)));