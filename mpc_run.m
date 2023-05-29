close all
%% Parameters for MPC
param = struct("Hu", 20, "Hp", 72, ...
               "mu1", 1, "mu2", 45, "mu3", 2, "mu4", 0.1, ...
               "mu5", 30, "mu6", 50, ...
               "alpha1", 1, "beta1", 0.75, ...
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

tic
log = mpc_base(param, options);
toc
%%
power = 326.6 * log.z(:,2) / 1000;
price_hour = power .* log.d2 / 60;

%%
plottime = options.SimulationStartDate:minutes(1):options.SimulationEndDate;
fig=figure('Position',[100 100 960 540]);
hold on
stairs(plottime,log.u)
stairs(plottime,log.z(:,1),"LineWidth",2)
stairs(plottime,log.d)
yline(22,"--","Color","black","LineWidth",2)
ax = gca;
hold on
box on
grid on
set(ax,'GridAlpha',0.3)
set(ax,'FontSize',14)
ylabel("Temperature [°C]")
legend(["Input","Room","Ambient","Reference"])
xlim([options.SimulationStartDate+hours(2) options.SimulationEndDate-days(1)])
exportgraphics(fig,"Results_Price_Cont.png","ContentType","vector")
%%
fig1=figure('Position',[100 100 960 540]);
    hold on
    stairs(plottime,log.z(:,1),"LineWidth",2)
    yline(22,"--","Color","black","LineWidth",2)
    ax = gca;
    hold on
    box on
    grid on
    set(ax,'GridAlpha',0.3)
    set(ax,'FontSize',14)
    ylabel("Temperature [°C]")
    legend(["Room","Reference"])
    xlim([options.SimulationStartDate+days(5) options.SimulationStartDate+days(10)])
    exportgraphics(fig1,"Results_Price_Cont_SS.pdf","ContentType","vector")
%%
fig1=figure('Position',[100 100 960 540]);
    hold on
    stairs(plottime,log.z(:,1),"LineWidth",2)
    yline(22,"--","Color","black","LineWidth",2)
    ax = gca;
    hold on
    box on
    grid on
    set(ax,'GridAlpha',0.3)
    set(ax,'FontSize',14)
    ylabel("Temperature [°C]")
    legend(["Room","Reference"])
    xlim([options.SimulationStartDate+hours(3)+minutes(30) options.SimulationStartDate+days(5)])
    ylim([20 22.5])
    exportgraphics(fig1,"Results_Price_Cont_Step.pdf","ContentType","vector")
%%
fig1=figure('Position',[100 100 960 540]);
    hold on
    stairs(plottime,log.u)
    stairs(plottime,log.z(:,1),"LineWidth",2)
    stairs(plottime,log.d)
    yline(22,"--","Color","black","LineWidth",2)
    ax = gca;
    hold on
    box on
    grid on
    set(ax,'GridAlpha',0.3)
    set(ax,'FontSize',14)
    ylabel("Temperature [°C]")
    legend(["Input","Room","Ambient","Reference"],"Location","best")
    xlim([options.SimulationStartDate+minutes(26) options.SimulationStartDate+days(5)])
    exportgraphics(fig1,"Results_Price_Cont_Start.pdf","ContentType","vector")
%% 
    fig2=figure('Position',[100 100 960 540]);
    subplot(311)
    plot(plottime,log.d2)
    ax = gca;
    hold on
    box on
    grid on
    set(ax,'GridAlpha',0.3)
    set(ax,'FontSize',14)
    ylabel("$m_e(k)$ [DKK/kWh]","Interpreter","latex")
    xlim([options.SimulationStartDate+hours(1) options.SimulationEndDate-days(1)])

    subplot(312)
    plot(plottime,power)
    ax = gca;
    hold on
    box on
    grid on
    set(ax,'GridAlpha',0.3)
    set(ax,'FontSize',14)
    ylabel("$P_c(k)$ [kW]","Interpreter","latex")
    xlim([options.SimulationStartDate+hours(1) options.SimulationEndDate-days(1)])

    subplot(313)
    plot(plottime,cumsum(price_hour))
    ax = gca;
    hold on
    box on
    grid on
    set(ax,'GridAlpha',0.3)
    set(ax,'FontSize',14)
    ylabel("$\bar{V}(k)$ [DKK]","Interpreter","latex")
    xlim([options.SimulationStartDate+hours(1) options.SimulationEndDate-days(1)])
    %exportgraphics(fig2,"Results_Price_Econ_Results.pdf","ContentType","vector")