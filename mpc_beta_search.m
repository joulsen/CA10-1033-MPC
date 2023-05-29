addpath("../../commoncode/")
%% Parameters for MPC
param = struct("Hu", 20, "Hp", 35, ...
               "mu1", 1, "mu2", 45, "mu3", 2, "mu4", 0.1, ...
               "mu5", 30, "mu6", 50, ...
               "alpha1", 1, "beta1", 1, ...
               "x0", [35; 30; 25; -2; 20], ...
               "P", 1*[0.01, 0.02, 0.03, 0.04, 0.05], ...
               "kappa", 30);
param.u0 = 30; % This cannot be less than ambient for constraint reasons
param.u_fallback = 32;
param.Q = kron(eye(param.Hp), [1, 0; 0, 0]);
param.R = 0.1*eye(param.Hu);
% Expand parameters with lifted matrices
param = mpc_initialize(param);
options = struct("MPC_Enabled", true, "Input_Default", 0, ...
                 "Simulation", true, ...
                 "SimulationStartDate", datetime(2022,1,1), ...
                 "SimulationEndDate", datetime(2023,1,1), ...
                 "Logging", false, "LoggingInterval", 3*60, ...
                 "LogOnCrash", false, ...
                  "PriceOn", true, ...
                 "PriceNormalisation", "cap", ...
                 "NormMean", 2, ...
                 "NormMin", 0.05, "NormMax", 2.2, ...
                 "qp_opt", optimoptions('quadprog', 'Algorithm', ...
                            'active-set', 'Display', 'off'));
%%
% wait = waitbar(0, "Simulating with different betas");
betas = linspace(0, 3, 25);
logs = cell(length(betas), 1);
for i = 1:length(betas)
    tic
    beta = betas(i);
    param_i = param;
    param_i.beta1 = beta;
%     waitbar(i / length(betas), wait, ...
%         sprintf("Simulating with beta = %f", beta));
    log = mpc_base(param_i, options);
    logs{i} = log;
    toc
end
%%
save("mpc_beta_iteration.mat", "logs", "betas")

% close(wait)

%
%plottime = options.SimulationStartDate:minutes(1):options.SimulationEndDate;
%fig = default_fig();
%plot(plottime,logs{1}.z(:,1),"Color",[0.75, 0.2118, 0.2118]);
%plot(plottime,log.u,"Color",[0.8745, 0.6980, 0.3333]);
%plot(plottime,logs{2}.z(:,1),"linewidth",2,"Color",[0.1059, 0.1059, 0.1059]);
%legend(["Room temp. with no price weight", "Ambient temperature", "Room temp. with excessive weight"],"Location","northwest");
%xlim([options.SimulationStartDate, options.SimulationEndDate-minutes(1)])
%ylabel("Temperature [°C]")
%exportgraphics(fig,"BetaImportance.pdf","ContentType","vector")

%clf
%yyaxis("left")
%hold on
%plot(logs{1}.z(:,1))
%plot(logs{1}.u)
%plot(logs{1}.d)
%yyaxis("right")
%plot(logs{1}.d2)
%plot(logs{1}.d2bar)
%legend(["Room temperature", "Input", "Ambient temperature", "Energy Price", "Normalised Energy Price"])