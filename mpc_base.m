clear mpc_iteration_base
clear get_forecast
clear get_forecast_historical
clear get_price_historical
%% Parameters for MPC
param = struct("Hu", 20, "Hp", 25, ...
               "mu1", 1, "mu2", 45, "mu3", 2, "mu4", 0.1, ...
               "x0", [35; 30; 25; -2; 20], ...
               "P", [0.01, 0.02, 0.03, 0.04, 0.05], ...
               "kappa", 30);
param.u0 = 30; % This cannot be less than ambient for constraint reasons
param.Q = eye(param.Hp);
param.R = 0.1*eye(param.Hu);
% Expand parameters with lifted matrices
param = mpc_initialize(param);
options = struct("MPC_Enabled", true, "Input_Default", 0, ...
                 "Simulation", true, ...
                 "SimulationStartDate", datetime(2022,1,1), ...
                 "SimulationEndDate", datetime(2023,1,1), ...
                 "Ploton",1);
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
if options.Simulation
    T_RUN = minutes(options.SimulationEndDate - options.SimulationStartDate);
else
    T_RUN = 60*24*60*60; % Allocating memory for 60 days for efficiency
end
u_log = [u; nan(T_RUN, 1)];
y_log = [y; nan(T_RUN, 1)];
x_log = [x'; nan(T_RUN, size(x, 1))];
z_log = [z; nan(T_RUN, 1)];
d_log = [0; nan(T_RUN, 1)];
Hp_log = [0; nan(T_RUN, 1)];
%% Control loop
% Run this loop below once a minute.
% <- INSERT ->
try
while k < T_RUN
    % Time update
    k = k + 1;
    %% Input
    % Reference value for indoor temperature every minute.
    % Reference must be at least Hp*Kappa in length
    % <- INSERT ->
    r = 22*ones(param.kappa*param.Hp, 1);
    % Disturbance is ambient temperature for every minute in a vector.
    % Length must be at least Kappa.
    if options.Simulation
        time = options.SimulationStartDate + seconds(k*sys.Ts);
        d = get_forecast_historical(1000, time, options.SimulationEndDate);
    else
        d = get_forecast(1000);
    end
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
        Hp_log(k) = Hp_bar;
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

%clear get_price_historical
% Simulation results
% power = 326.6 * (x_log(:,1) - x_log(:, 3)) / 1000; % kW
power = 326.6 * (u_log - x_log(:, 3)) / 1000; % kW
price = get_price_historical(-1, options.SimulationStartDate, options.SimulationEndDate); % DKK/KWh
price_total = power .* price / 60;

%% Plotting spas

if options.Ploton
    plottime = options.SimulationStartDate:minutes(1):options.SimulationEndDate;
    tot_price = cumsum(price_total);

%     figure(1)
%     subplot(211)
%     hold on
%     plot(plottime,x_log)
%     plot(plottime,d_log)
%     ax = gca;
%     hold on
%     box on
%     grid on
%     set(ax,'GridAlpha',0.3)
%     set(ax,'FontSize',14)
%     ylabel("States [°C]")
%     xlim([options.SimulationStartDate options.SimulationEndDate-days(1)])

    fig=figure('Position',[100 100 960 540]);
    hold on
    stairs(plottime,u_log)
    stairs(plottime,z_log,"LineWidth",2)
    stairs(plottime,d_log)
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
    exportgraphics(fig,"Comfort_Room_Temp_Results.png","ContentType","vector")
    
    fig1=figure('Position',[100 100 960 540]);
    hold on
    stairs(plottime,z_log,"LineWidth",2)
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
    ylim([21.9 22.1])
    exportgraphics(fig1,"Comfort_steadystate_results.pdf","ContentType","vector")  

    fig1=figure('Position',[100 100 960 540]);
    hold on
    stairs(plottime,z_log,"LineWidth",2)
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
    exportgraphics(fig1,"Comfort_step_results.pdf","ContentType","vector") 

    fig1=figure('Position',[100 100 960 540]);
    hold on
    stairs(plottime,u_log)
    stairs(plottime,z_log,"LineWidth",2)
    stairs(plottime,d_log)
    yline(22,"--","Color","black","LineWidth",2)
    ax = gca;
    hold on
    box on
    grid on
    set(ax,'GridAlpha',0.3)
    set(ax,'FontSize',14)
    ylabel("Temperature [°C]")
    legend(["Input","Room","Ambient","Reference"])
    xlim([options.SimulationStartDate+minutes(26) options.SimulationStartDate+days(5)])
    exportgraphics(fig1,"Comfort_Startbehave_results.pdf","ContentType","vector")    

    fig2=figure('Position',[100 100 960 540]);
    subplot(311)
    plot(plottime,price)
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
    plot(plottime,tot_price)
    ax = gca;
    hold on
    box on
    grid on
    set(ax,'GridAlpha',0.3)
    set(ax,'FontSize',14)
    ylabel("$\bar{V}(k)$ [DKK]","Interpreter","latex")
    xlim([options.SimulationStartDate+hours(1) options.SimulationEndDate-days(1)])
    exportgraphics(fig2,"Comfort_Econ_Results.pdf","ContentType","vector")
end