function [log] = mpc_base(param, options)
% MPC_BASE Runs the MPC controller loop and returns the logs of states,
% inputs, outputs etc.
    % Clearing functions with persistents
    clear mpc_iteration_base
    clear get_forecast
    clear get_forecast_historical
    clear get_price
    clear get_price_historical
    sys = param.sys;
    psys = param.psys;
    u = param.u0;
    x = param.x0;
    y = sys.Cy * x;
    x_est = param.x0;
    z = [0, 0];
    k = 0;
    ke = 0;
    mpc_period_start = 1;
    dim = struct("p", size(param.msys.Cz, 1), ...
                 "m", size(param.msys.A, 1), ...
                 "n", size(param.msys.B, 2), ...
                 "j", size(param.msys.Bd, 2));
    % Run-time for pre-allocating memory for speed
    if options.Simulation
        T_RUN = minutes(options.SimulationEndDate - options.SimulationStartDate);
        times = options.SimulationStartDate:seconds(sys.Ts):options.SimulationEndDate;
    else
        T_RUN = 60*24*365*5; % 5 years
        % Connect to device
        javaaddpath(pwd);
        javaaddpath([pwd '/lmtclient']);
        cl = LMTClientClass(REDACTED);
        cl.Connect()
    end
    % Logging
    log = struct("u",     [u; nan(T_RUN, 1)], ...
                 "y",     [y; nan(T_RUN, 1)], ...
                 "x",     [x'; nan(T_RUN, size(x, 1))], ...
                 "x_est", [x_est'; nan(T_RUN, size(x_est, 1))], ...
                 "z",     [z; nan(T_RUN, dim.p)], ...
                 "d",     [0; nan(T_RUN, 1)], ...
                 "d2",    [0; nan(T_RUN, 1)], ...
                 "d2bar", [0; nan(T_RUN, 1)], ...
                 "Hp",    [0; nan(T_RUN, 1)]);
    %% Main control loop
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
            d = get_forecast_historical(2200, times(k), options.SimulationEndDate);
            [d2, ke] = get_price_historical(times(k), options);
        else
            d = get_forecast(2200);
            d2 = get_price(datetime("now", "TimeZone", "UTC"));
        end
        log.d2(k, :) = d2(k);
        if ~options.PriceOn
            d2 = zeros(T_RUN, 1);
            ke = param.Hp * param.kappa;
        end
%         d2 = rmmissing(d2);
%         d2 = d2(~isnan(d2));
        %% Controller
        if (~mod(k-1, param.kappa))
            mpc_period_start = k;
            l1 = floor(size(rmmissing(d), 1));
            if options.Simulation
                Hp_bar = min(param.Hp, ke / param.kappa);
                Hp_bar = min(Hp_bar, l1 / param.kappa);
            else
                l2 = floor(size(rmmissing(d2), 1));
                Hp_bar = min(min(param.Hp, l1  / param.kappa), ...
                             min(param.Hp, l2 / param.kappa));
            end
            Hp_bar=floor(Hp_bar);
            Hu_bar = min(param.Hu, Hp_bar);
            dbar = d(1:param.kappa:Hp_bar*param.kappa);
            if options.Simulation
                iend = min(k+Hp_bar*param.kappa-1, k+ke-1);
                d2bar = d2(k:param.kappa:iend);
            else
                d2bar = d2(1:param.kappa:Hp_bar*param.kappa);
            end
            if (options.PriceNormalisation ~= "none" && options.PriceOn)
                if options.PriceNormalisation == "squish"
                    d2bar = (d2bar - min(d2bar));
                    if ~any(d2bar)
                        d2bar = ones(size(d2bar));
                    else
                        d2bar = d2bar / max(d2bar);
                    end
                elseif options.PriceNormalisation == "mean"
                    d2bar = d2bar - mean(d2bar) + options.NormMean;
                elseif options.PriceNormalisation == "cap"
                    d2bar = min(max(d2bar, options.NormMin), options.NormMax);
                end
            end
            rmask = [1; zeros(param.kappa-1, 1)];
            rmask = kron(rmask, ones(dim.p, 1));
            rmask = kron(ones(Hp_bar, 1), rmask) == 1;
            rbar = r(rmask);
            price_Q = ones(Hp_bar * dim.p, 1);
            price_Q(1:dim.p:end) = param.alpha1;
            price_Q(2:dim.p:end) = d2bar.^2 * param.beta1;
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
                us = zeros(Hu_bar*param.kappa, 1);
                us(k - mpc_period_start + 1:end) = u + deltaU;
            end
        end
        % Input is selected
        u = us(k - mpc_period_start + 1);
        if ~options.Simulation
            set = LimitValues(u, 15, 40); % Limit value just in case
            % In this mode the heat pump always use the min setpoint
            cl.SetVariables({'Heating.SetPointMin'}, Tset)
        end
        %% Observer
        % Measurement of the return temperature
        if ~options.Simulation
            CurrentTemp = cl.GetVariables({'Input.THeatReturn'});
            if(CurrentTemp > 0) % Only update input if get was succesful
                y = CurrentTemp;
            end
        end
        y_est = sys.Cy * x_est + sys.Dy * u + sys.Ddy * d(1);
        x_est = sys.A*x_est + param.L*(y-y_est) + sys.B * u + sys.Bd * d(1);
        % Model
        x = psys.A * x + psys.B * u + psys.Bd * d(1,:)';
        z = psys.Cz * x + psys.Dz * u + psys.Ddz * d(1,:)';
        if options.Simulation
            y = psys.Cy * x + psys.Dy * u + psys.Ddy * d(1,:)';
        end
        % Logging
        log.x(k, :) = x';
        log.x_est(k, :) = x_est';
        log.u(k, :) = u;
        log.y(k, :) = y;
        log.z(k, :) = z;
        log.d(k, :) = d(1);
        log.d2bar(k, :) = d2bar(1);
        log.Hp(k) = Hp_bar;
        if (options.Logging && ~mod(k, options.LoggingInterval))
            now = datetime("now", "Format", "yyyy_MM_dd_HH.mm");
            fname = strcat("logs/", char(now), ".mat");
            save(fname, "log");
        end
    end
    % Save data in case of error and rethrow
    catch ME
        if options.LogOnCrash
            now = datetime("now", "Format", "yyyy_MM_dd_HH.mm");
            fname = strcat("logs/", char(now), "_crash", ".mat");
            save(fname, "log");
        end
        rethrow(ME);
    end
end
