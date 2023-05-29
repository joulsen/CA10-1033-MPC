addpath("../../commoncode/")
%% Parameters for MPC
param = struct("Hu", 20, "Hp", 35, ...
               "mu1", 1, "mu2", 45, "mu3", 2, "mu4", 0.1, ...
               "mu5", 30, "mu6", 50, ...
               "alpha1", 1, "beta1", 0.375, ...
               "x0", [35; 30; 25; -2; 20], ...
               "P", 1*[0.01, 0.02, 0.03, 0.04, 0.05], ...
               "kappa", 30);
param.u0 = 30; 
param.u_fallback = 32;
param.Q = kron(eye(param.Hp), [1, 0; 0, 0]);
param.R = 0.1*eye(param.Hu);
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

options = {options, options};
options{1}.PriceOn = false;
options{1}.PriceNormalisation = "none";
options{1}.Name = "Comfort Controller";
options{2}.Name = "Economic Controller";
params = {param, param};
params{1}.Hp = 25;

start_year = 2018;
end_year = 2023;
dates = cell(end_year - start_year, 1);
for i=1:end_year - start_year + 1
    dates{i} = datetime(end_year - i + 1, 01, 01);
end
%%
logs = cell(length(dates), length(options));

for i=1:size(dates, 1) - 1
    fprintf("From %s to %s\n", dates{i+1}, dates{i})
    for j=1:length(options)
        tic
        options{j}.SimulationStartDate = dates{i+1};
        options{j}.SimulationEndDate = dates{i};
        log = mpc_base(params{j}, options{j});
        logs{i,j} = log;
%         if options{j}.Name == "Comfort Controller"
%             log.d2 = get_price_historical_all(-1, options{j}.SimulationStartDate, options{j}.SimulationEndDate);
%         end
        power = 326.6 * log.z(:,2) / 1000;
        price = power .* log.d2 / 60;
        duration = toc;
        fprintf("%s: %f DKK, %f °C (%fs)\n", options{j}.Name, ...
                sum(rmmissing(price)), min(rmmissing(log.z(1500:end,1))), duration);
    end
    fprintf("\n")
end
%%
prices = cell(size(logs));
prices_mask = cell(size(logs));
total_prices = zeros(size(logs));
total_prices_mask = zeros(size(logs));
for i=1:size(logs, 1)-1
    for j=1:size(options, 2)
        prices{i,j} = 326.6 * logs{i,j}.z(:,2) / 1000 / 60 .* logs{i,j}.d2;
        total_prices(i,j) = sum(rmmissing(prices{i,j}));
        mask = logs{i,j}.d > 15;
        prices_mask{i,j} = 326.6 * logs{i,j}.z(~mask,2) / 1000 / 60 .* logs{i,j}.d2(~mask);
        total_prices_mask(i,j) = sum(rmmissing(prices_mask{i,j}));
    end
end
reductions = (total_prices(:,1) - total_prices(:,2)) ./ total_prices(:,1) * 100;
reductions_mask = (total_prices_mask(:,1) - total_prices_mask(:,2)) ./ total_prices_mask(:,1) * 100;

%%
time = dates{2}:seconds(60):dates{1};
fig1 = default_fig();
hold on
plot(time, logs{1,1}.z(:,1));
plot(time, logs{1,2}.z(:,1));
legend(["Room temperature (Comfort)", "Room temperature (Economic)"], ...
        "Location", "northwest")
xlabel("Time")
ylabel("Temperature [°C]")
exportgraphics(fig1, "figures/MultiYear2023.pdf", "ContentType", "vector");

time = dates{5}:seconds(60):dates{4};
fig2 = default_fig();
hold on
plot(time, logs{5,1}.z(:,1));
plot(time, logs{5,2}.z(:,1));
legend(["Room temperature (Comfort)", "Room temperature (Economic)"], ...
       "Location", "northwest")
xlabel("Time")
ylabel("Temperature [°C]")
exportgraphics(fig2, "figures/MultiYear2019.pdf", "ContentType", "vector");
