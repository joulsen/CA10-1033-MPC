function d = get_forecast_historical(limit, startdate, enddate)
%GET_FORECAST_HISTORICAL Summary of this function goes here
%   Detailed explanation goes here
    persistent w k
    if (isempty(w))
        k = 1;
        url = "https://archive-api.open-meteo.com/v1/archive?latitude=REDACTED&longitude=REDACTED&start_date=%s&end_date=%s&hourly=temperature_2m&timezone=Europe/Berlin";
        startdate.Format = 'yyyy-MM-dd';
        enddate.Format = 'yyyy-MM-dd';
        weather = webread(sprintf(url, startdate, enddate), ...
                          weboptions("Timeout", 20));
        weathertime = datetime(weather.hourly.time, "InputFormat", 'yyyy-MM-dd''T''HH:mm');
        weather.hourly.time = weathertime;
        realtime = startdate:minutes(1):enddate;
        w = interp1(weather.hourly.time, weather.hourly.temperature_2m, realtime)';
    end
    limit = min(k+limit, size(w,1));
    d = w(k:limit);
    k = k + 1;
end

