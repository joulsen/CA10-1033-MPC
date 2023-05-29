clear mpc_iteration_base
clear get_forecast
%% Parameters for MPC
param = struct("Hu", 20, "Hp", 25, ...
               "mu1", 1, "mu2", 60, "mu3", 2, "mu4", 0.1, ...
               "x0", [35; 30; 25; -2; 20], ...
               "P", [0.01, 0.02, 0.03, 0.04, 0.05], ...
               "kappa", 30);
param.u0 = 30; % This cannot be less than ambient for constraint reasons
param.Q = eye(param.Hp);
param.R = 0.1*eye(param.Hu);
% Expand parameters with lifted matrices
param = mpc_initialize(param);
options = struct("MPC_Enabled", true, "Input_Default", 0);
sys = param.sys;

% Initial parameters
u = param.u0;
y = 0;
x = param.x0;
x_est = param.x0;
k = 0;
z = 0;
mpc_period_start = 1;

% Logging variables
T_RUN = 60*24*60*60; % Allocating memory for 60 days for efficiency
u_log = [u; nan(T_RUN, 1)];
y_log = [y; nan(T_RUN, 1)];
x_log = [x'; nan(T_RUN, size(x, 1))];
z_log = [z; nan(T_RUN, 1)];
d_log = [0; nan(T_RUN, 1)];
%% Control loop
% Run this loop below once a minute.
% <- INSERT ->
try
while true
    % Time update
    k = k + 1;
    %% Input
    % Reference value for indoor temperature every minute.
    % Reference must be at least Hp*Kappa in length
    % <- INSERT ->
    r = 22*ones(param.kappa*param.Hp, 1);
    % Disturbance is ambient temperature for every minute in a vector.
    % Length must be at least Kappa.
    d = get_forecast(param.kappa*param.Hp+250);
    %% Observer
    % Measurement of the return temperature
    % <- INSERT ->
    % y = ...
    y_est = sys.Cy * x + sys.Dy * u + sys.Ddy * d(1);
    x_est = sys.A * (x + + param.L*(y - y_est)) + sys.B * u + sys.Bd * d(1);
    %% Controller
    if (~mod(k-1, param.kappa))
        mpc_period_start = k;
        [dbar, Hp_bar] = match_Hp(d, param.Hp, param.kappa);
        Hu_bar = min(param.Hu, Hp_bar);
        rbar = r(1:param.kappa:param.kappa*Hp_bar);
        if Hp_bar == 0
            error("Error: No values in disturbance vector.")
        end
        deltaU = mpc_iteration(x_est, u, rbar, dbar, Hp_bar, Hu_bar, param, options);
        if size(deltaU, 1) > 1
            us = u + deltaU;
            T = (Hu_bar-1)*param.kappa;
            us = interp1(0:param.kappa:T, us, 0:T, "previous")';
        else
            us = zeros(size(us));
            us(k:end) = u + deltaU;
        end
    end
    % Input is selected
    u = us(k - mpc_period_start + 1);
    % Send input u to appropriate method
    % <- INSERT ->
    % Logging
    x_log(k, :) = x';
    u_log(k, :) = u;
    y_log(k, :) = y;
    z_log(k, :) = z;
    d_log(k, :) = d(1);
    % Model
    x = sys.A * x + sys.B * u + sys.Bd * d(1,:)';
    y = sys.Cy * x + sys.Dy * u + sys.Ddy * d(1,:)';
    z = sys.Cz * x + sys.Dz * u + sys.Ddz * d(1,:)';
    if ~mod(k, 12*60*60)
        now = datetime("now", "Format", "yyyy_MM_dd_HH.mm");
        fname = strcat("logs/", char(now), ".mat");
        save(fname, "x_log", "y_log", "z_log", "u_log", "d_log");
    end
end
% Save data in case of error and rethrow
catch ME
    now = datetime("now", "Format", "yyyy_MM_dd_HH.mm");
    fname = strcat("logs/", char(now), "_crash", ".mat");
    save(fname, "x_log", "y_log", "z_log", "u_log", "d_log");
    rethrow(ME);
end