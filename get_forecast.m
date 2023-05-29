function d = get_forecast(limit)
%GET_FORECAST Summary of this function goes here
%   Detailed explanation goes here
    persistent w k
    if (isempty(w) || k > 6*60)
        k = 1;
        weather = webread("https://api.open-meteo.com/v1/forecast?latitude=REDACTED&longitude=REDACTED&hourly=temperature_2m", ...
                          weboptions("Timeout", 20));
        weathertime = datetime(weather.hourly.time, "InputFormat", 'yyyy-MM-dd''T''HH:mm');
        weather.hourly.time = weathertime;
        realtime = datetime("now"):minutes(1):weather.hourly.time(end);
        w = interp1(weather.hourly.time, weather.hourly.temperature_2m, realtime)';
    end
    k = k + 1;
    d = w(k:k + limit);
end

