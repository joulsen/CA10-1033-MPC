addpath("..\..\commoncode\")
clear mpc_simulation

data = read_bitzer_data("../../BitzerData/Logdata-d-60/HPLog 026.csv", [1, -1]);
steps = size(data, 1);
d = data.Tamb;

% Defining parameters for the MPC
param = struct("Hu", 20, "Hp", 1000, ...
               "mu1", 1, "mu2", 60, "mu3", 2, "mu4", 0.1, ...
               "x0", [35; 30; 25; -2; 20], ...
               "P", [0.01, 0.02, 0.03, 0.04, 0.05], ...
               "kappa", 30);
param.u0 = d(1) + param.mu3;
param.Q = eye(param.Hp);
param.R = 10*eye(param.Hu);
param = mpc_iteration_initialize(param);
sys = param.sys;

% Defining other options for the MPC object
options = struct("MPC_Enabled", true, "Input_Default", 0);
%% Initialize simulation
% u = zeros(steps, 1);
r = ones(steps, 1)*25;
tic
[u, x, xest, y] = mpc_simulation(r, d, param, options);
toc

clf
stairs(x)