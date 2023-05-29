addpath("..\..\commoncode\")
clear mpc_iteration

data = read_bitzer_data("../../BitzerData/Logdata-d-60/HPLog 026.csv", [1, -1]);
steps = size(data, 1);
u = zeros(steps, 1);
r = ones(steps, 1)*25;
d = data.Tamb;

% Defining parameters for the MPC
param = struct("Hu", 20, "Hp", 1000, ...
               "mu1", 1, "mu2", 60, "mu3", 2, "mu4", 1, ...
               "x0", [35; 30; 25; -2; 20], ...
               "P", [0.01, 0.02, 0.03, 0.04, 0.05]);
param.u0 = d(1) + param.mu3;
param.Q = eye(param.Hp);
param.R = zeros(param.Hu);
param = mpc_iteration_initialize(param);

% Defining other options for the MPC object
options = struct("MPC_Enabled", true, "Input_Default", 0);

[u, x] = mpc_iteration_base(25, d, 20, param, options);
