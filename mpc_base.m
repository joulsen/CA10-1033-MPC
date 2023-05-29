addpath("..\..\commoncode\")
data = read_bitzer_data("../../BitzerData/Logdata-d-60/HPLog 026.csv", [1, -1]);
d = data.Tamb(1:steps);
clear mpc_iteration_base
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
param = mpc_iteration_initialize(param);
options = struct("MPC_Enabled", true, "Input_Default", 0);
sys = param.sys;

% Initial parameters
u = param.u0;
y = 0;
x = param.x0;
x_est = param.x0;
k = 0;
z = 0;

% Logging variables
u_log = [u];
y_log = [y];
x_log = [x'];
z_log = [z];
%% Control loop
% Run this loop below once a minute.
% <------------------------------->
while true
    % Time update
    k = k + 1;
    %% Input
    % Reference value for indoor temperature every minute.
    % Reference must be at least Hp*Kappa in length
    r = 25*ones(param.kappa*param.Hp, 1);
    % Disturbance is ambient temperature for every minute in a vector.
    % Length must be at least Kappa.
    % <------------------------>
    d = d(2:end);
    %% Observer
    % Measurement of the return temperature
    % <------------------------>
    % y = ...
    y_est = sys.Cy * x + sys.Dy * u + sys.Ddy * d(1);
    x_est = sys.A * x + sys.B * u + sys.Bd * d(1) + param.L*(y - y_est);
    %% Controller
    if (~mod(k-1, param.kappa))
        [dbar, Hp_bar] = match_Hp(d, param.Hp, param.kappa);
        Hu_bar = min(param.Hu, Hp_bar);
        rbar = r(1:param.kappa:param.kappa*Hp_bar);
        if Hp_bar == 0
            error("Error: No values in disturbance vector.")
        end
        deltaU = mpc_iteration_base(x_est, u, rbar, dbar, Hp_bar, Hu_bar, param, options);
        if size(deltaU, 1) > 1
            us = u + deltaU;
            us = interp1(k + (0:Hu_bar-1)*param.kappa, us, 1:steps, "previous")';
        else
            us = zeros(size(us));
            us(k:end) = u + deltaU;
        end
    end
    % Input is selected
    u = us(k);
    % Send input to appropriate method
    % <--------------------->
    % Logging
    x_log = [x_log; x'];
    u_log = [u_log; u];
    y_log = [y_log; y];
    z_log = [z_log; z];
    % Model
    x = sys.A * x + sys.B * u + sys.Bd * d(1,:)';
    y = sys.Cy * x + sys.Dy * u + sys.Ddy * d(1,:)';
    z = sys.Cz * x + sys.Dz * u + sys.Ddz * d(1,:)';
end


function [mat, Hp_bar] = match_Hp(mat, Hp, kappa)
    v = size(mat, 1);
    Hp_bar = min(Hp, floor(v / kappa));
    mat = mat(1:kappa:Hp_bar*kappa);
end
function deltaU = mpc_iteration_base(x, u, r, d, Hp_bar, Hu_bar, param, options)
    %MPC_ITERATION_BASE Performs a single iteration of the MPC
    if options.MPC_Enabled
        % Cropping to Hp_bar
        dim = struct("p", size(param.msys.Cz, 1), ...
                     "m", size(param.msys.A, 1), ...
                     "n", size(param.msys.B, 2), ...
                     "j", size(param.msys.Bd, 2));
        Psi = param.Psi(1:dim.p*Hp_bar, :);
        Ups = param.Ups(1:dim.p*Hp_bar, :);
        Theta = param.Theta(1:dim.p*Hp_bar, 1:dim.n*Hu_bar);
        Xi = param.Xi(1:dim.p*Hp_bar, 1:dim.j*Hp_bar);
        Q = param.Q(1:Hp_bar, 1:Hp_bar);
        R = param.R(1:Hu_bar, 1:Hu_bar);
        H = Theta' * Q * Theta + R;
        H = (H + H')/2;
        Eps = r - Psi * x - Ups * u - Xi * d;
        G = 2*Theta' * Q * Eps;
        [A_con, b_con] = get_lifted_constraints(Hu_bar, u, d);
        qp_opt = optimoptions('quadprog', 'Algorithm', 'active-set', ...
                              'Display', 'off');
        deltaU = quadprog(2*H, -G', A_con, b_con, ...
                          [],[],[],[],zeros(Hu_bar, 1),qp_opt);
    else
        deltaU = 0;
    end
    function [mat, Hp_bar] = match_Hp(mat, Hp, kappa)
        v = size(mat, 1);
        Hp_bar = min(Hp, floor(v / kappa));
        mat = mat(1:kappa:Hp_bar*kappa);
    end
    function [A_con, b_con] = get_lifted_constraints(Hu, u0, d)
        nu = 1;
        nf = 3;
        % Input constraint F
        F = kron(eye(Hu, Hu), [-1; 1; -1]);
        f1 = kron(ones(Hu, 1), [param.mu1; -param.mu2; param.mu3]);
        f2 = kron(eye(Hu), [0; 0; 1]) * d(1:Hu);
        f = f1+f2;
        F = [F, f];
        % Input rate constrain E = [W w]
        W = [eye(Hu); -eye(Hu)];
        w = [ones(Hu, 1); ones(Hu, 1)] * param.mu4 * param.kappa;
        % Lifting input constraint F
        calF = zeros(Hu*nf, Hu*nu);
        for i=1:Hu
            calFi = zeros(Hu*nf, nu);
            for j=i:Hu
                Fi = F(:,nu*j);
                calFi = calFi + Fi;
            end
            calF(:,i) = calFi;
        end
        A_con = [calF; W];
        b_con = [-calF(:,1:nu) * u0 - f; w];
    end
end