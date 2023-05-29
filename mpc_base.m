clear mpc_iteration_base
clear get_forecast
clear get_forecast_historical
clear get_price_v2
%% Parameters for MPC
param = struct("Hu", 20, "Hp", 35, ...
               "mu1", 1, "mu2", 45, "mu3", 2, "mu4", 0.1, ...
               "beta1", 0.01, ...
               "x0", [35; 30; 25; -2; 20], ...
               "P", [0.01, 0.02, 0.03, 0.04, 0.05], ...
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
                 "SimulationEndDate", datetime(2022,1,14), ...
                 "Ploton",0);
sys = param.sys;

% Initial parameters
u = param.u0;
y = 0;
x = param.x0;
x_est = param.x0;
k = 0;
z = [0, 0];
mpc_period_start = 1;

% Logging variables
if options.Simulation
    T_RUN = minutes(options.SimulationEndDate - options.SimulationStartDate);
else
    T_RUN = 60*24*365*5; % Allocating memory for 5 year for efficiency
end
dim = struct("p", size(param.msys.Cz, 1), ...
             "m", size(param.msys.A, 1), ...
             "n", size(param.msys.B, 2), ...
             "j", size(param.msys.Bd, 2));
u_log = [u; nan(T_RUN, 1)];
y_log = [y; nan(T_RUN, 1)];
x_log = [x'; nan(T_RUN, size(x, 1))];
x_est_log = [x_est'; nan(T_RUN, size(x, 1))];
z_log = [z; nan(T_RUN, dim.p)];
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
    r = kron(ones(param.kappa*param.Hp, 1), [22; 0]);
    % Disturbance is ambient temperature for every minute in a vector.
    % Length must be at least Kappa.
    if options.Simulation
        time = options.SimulationStartDate + seconds(k*sys.Ts);
        d = get_forecast_historical(1000, time, options.SimulationEndDate);
        d2 = get_price_v2(time);
    else
        d = get_forecast(1000);
        d2 = get_price_v2(datetime("now"));
    end
    d2 = rmmissing(d2);
    %% Controller
    if (~mod(k-1, param.kappa))
        mpc_period_start = k;
        l1 = floor(size(rmmissing(d), 1));
        l2 = floor(size(rmmissing(d2), 1));
        Hp_bar = min(min(param.Hp, l1  / param.kappa), ...
                     min(param.Hp, l2 / param.kappa));
        Hp_bar=floor(Hp_bar);
        %disp([l1, l2, Hp_bar])
        Hp_log(k) = Hp_bar;
        Hu_bar = min(param.Hu, Hp_bar);
        dbar = d(1:param.kappa:Hp_bar*param.kappa);
        d2bar = d2(1:param.kappa:Hp_bar*param.kappa);
%         [dbar, Hp_bar] = match_Hp(d, param.Hp, param.kappa);
        rmask = [1; zeros(param.kappa-1, 1)];
        rmask = kron(rmask, ones(dim.p, 1));
        rmask = kron(ones(Hp_bar, 1), rmask) == 1;
        rbar = r(rmask);
        price_Q = ones(Hp_bar * dim.p, 1);
        price_Q(2:dim.p:end) = d2bar * param.beta1;
        param.Q = diag(price_Q);
        if Hp_bar == 0
            error("Error: No values in disturbance vector.")
        end
        deltaU = mpc_iteration(x_est, u, rbar, dbar, Hp_bar, Hu_bar, param, options);
        if size(deltaU, 1) > 1
            us = u + deltaU;
            T = (Hu_bar-1)*param.kappa;
            us = interp1(0:param.kappa:T, us, 0:T, "previous")';
        else
            us = zeros(Hu_bar*param.kappa);
            us(k - mpc_period_start + 1:end) = u + deltaU;
        end
    end
    % Input is selected
    u = us(k - mpc_period_start + 1);
    % Send input u to appropriate method
        %% Observer
    % Measurement of the return temperature
    % <- INSERT ->
    % y = ...
    y_est = sys.Cy * x_est + sys.Dy * u + sys.Ddy * d(1);
    x_est = sys.A * x_est + param.L*(y - y_est) + sys.B * u + sys.Bd * d(1);
    % <- INSERT ->
    % Logging
    x_log(k, :) = x';
    x_est_log(k,:) = x_est';
    u_log(k, :) = u;
    y_log(k, :) = y;
    z_log(k, :) = z;
    d_log(k, :) = d(1);
    % Model
    x = sys.A * x + sys.B * u + sys.Bd * d(1,:)';
    y = sys.Cy * x + sys.Dy * u + sys.Ddy * d(1,:)';
    z = sys.Cz * x + sys.Dz * u + sys.Ddz * d(1,:)';
    if ~mod(k, 3*60)
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

clear get_price_historical
% Simulation results
% power = 326.6 * (x_log(:,1) - x_log(:, 3)) / 1000; % kW
power = 326.6 * (u_log - x_log(:, 3)) / 1000; % kW
price = get_price_historical(-1, options.SimulationStartDate, options.SimulationEndDate); % DKK/KWh
price_total = power .* price / 60;

%% Plotting spas

if options.Ploton
    plottime = options.SimulationStartDate:minutes(1):options.SimulationEndDate;
    tot_price = cumsum(price_total);

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
%     exportgraphics(fig,"Comfort_Room_Temp_Results.png","ContentType","vector")
    
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
%     exportgraphics(fig1,"Comfort_steadystate_results.pdf","ContentType","vector")  

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
%     exportgraphics(fig1,"Comfort_step_results.pdf","ContentType","vector") 

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
%     exportgraphics(fig1,"Comfort_Startbehave_results.pdf","ContentType","vector")    

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
%     exportgraphics(fig2,"Comfort_Econ_Results.pdf","ContentType","vector")
end