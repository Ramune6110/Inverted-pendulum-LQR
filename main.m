clear;
close all;
clc;

% paramater
param.M = 0.7;
param.m = 0.12;
param.J = 0.009;
param.l = 0.3;
param.g = 9.8;
param.alpha = (param.M + param.m) * (param.J + param.m * param.l^2) - param.m^2 * param.l^2;

%time
ts       = 0.0;
dt       = 0.1;
tf       = 25;
param.t  = ts:dt:tf;

% state
param.A   = [0, 0, 1, 0;
             0, 0, 0, 1;
             0, -(param.m^2 * param.l^2 * param.g) / param.alpha, 0, 0;
             0, ((param.M + param.m) * param.m * param.g * param.l) / param.alpha, 0, 0];
param.B   = [0; 0; (param.J * param.m * param.l^2) / param.alpha; -(param.m * param.l) / param.alpha];
param.C   = [1, 0, 0, 0];
Q         = diag([10, 10, 10, 10]);
R         = 1;
x         = [1;0;0;0];
xEst      = [1;0;0;0];
xEst_sec  = [1;0;0;0];
PEst      = diag([0.1, 0.1, 0.1, 0.1]);
u         = 0;
% save data
state_log     = zeros(length(param.t), length(x));
xEst_log      = zeros(length(param.t), length(x));
input_log     = zeros(length(param.t), length(u));
uEst_log      = zeros(length(param.t), length(u));
xEst_sec_log  = zeros(length(param.t), length(x));

i = 1;
tic;
for t = ts:dt:tf
    % input
    [K, S, P] = lqr(param.A, param.B, Q, R);
    u    = -K * x;
    uEst = -K * xEst;
    % kalmanfilter
    [xEst, PEst] =  KalmanFilter(x, u, param, PEst, dt);
    % -- runge-kutta --
    x        = Runge_Kutta(x, u, param, dt);
    xEst_sec = Runge_Kutta(xEst, uEst, param, dt);
  
    % -- save data --
    state_log(i, :)     = x;
    xEst_log(i, :)      = xEst;
    input_log(i, :)     = u;
    uEst_log(i, :)      = uEst;
    xEst_sec_log(i, :)  = xEst_sec;
    i = i + 1;
end
toc

DrowGraph(param, state_log, input_log, xEst_log, uEst_log, xEst_sec_log);

function x = model(x, u, param)
    x = param.A * x + param.B * u;
end

function x = Runge_Kutta(x, u, param, dt)
    % -- runge-kutta --
    k1 = model(x, u, param);
    k2 = model(x + k1 * dt / 2, u, param);
    k3 = model(x + k2 * dt / 2, u, param);
    k4 = model(x + dt * k3, u, param);
    x  = x + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6;
end

function [xEst, PEst] = KalmanFilter(xTrue, u, param, PEst, dt)
            % noise
            Q = 0.00000001;
            R = 0.00000001;
            % Prediction Step
            z     = param.C * xTrue + randn() * sqrt(R);
            xPred = Runge_Kutta(xTrue, u, param, dt) + randn() * sqrt(Q);
            PPred = param.A * PEst * param.A' + param.B * Q * param.B';
            % Filtering Step
            K    = PPred * param.C' / (param.C * PPred * param.C' + R);
            zEst = param.C * xPred;
            xEst = xPred + K * (z - zEst);
            PEst = (eye(size(param.A)) - K * param.C) * PPred;
end

function [] = DrowGraph(param, x, u, xEst, uEst, xEst_sec)
    figure(1);
    for i = 1:4
        subplot(2, 2, i);
        plot(param.t, x(:, i));
        if i == 1
            legend('True');
            title('LQR Control');
        end
    end
    figure(2);
    plot(param.t, u, 'b', param.t, uEst, 'r');
    legend('True', 'Estimate');
    title('Control Inout');
    figure(3);
    for i = 1:4
        subplot(2, 2, i);
        plot(param.t, x(:, i), 'b', param.t, xEst(:, i), 'r');
        if i == 1
            legend('True', 'Estimate');
            title('State estimate');
        end
    end
    figure(4);
    for i = 1:4
        subplot(2, 2, i);
        plot(param.t, x(:, i), 'b', param.t, xEst_sec(:, i), 'r');
        if i == 1
            legend('True', 'Estimate');
            title('LQR Control');
        end
    end
end
