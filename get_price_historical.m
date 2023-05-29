function d = get_price_historical(time, options)
%GET_PRICE_HISTORICAL Summary of this function goes here
%   Detailed explanation goes here
    persistent w k kn webopt pricetime_utc pricetime_dk realtime price
    if (isempty(webopt))
        webopt = weboptions("Timeout", 20);
    end
    if (isempty(w))
        k = 1;
        kn = 0;
        url = "https://api.energidataservice.dk/dataset/Elspotprices?offset=0&start=%sT00:00&end=%sT00:00%s";
        extra = "&filter=%7B%22PriceArea%22:[%22DK1%22]%7D&sort=HourUTC%20DESC&timezone=dk";
        startdate = options.SimulationStartDate;
        enddate = options.SimulationEndDate;
        startdate.Format = 'yyyy-MM-dd';
        enddate.Format = 'yyyy-MM-dd';
        price = webread(sprintf(url, startdate, enddate + days(1), extra), ...
                        webopt);
        price = struct2table(price.records);
        pricetime_utc = datetime(price.HourUTC, "InputFormat", 'yyyy-MM-dd''T''HH:mm:ss');
        pricetime_dk = datetime(price.HourDK, "InputFormat", 'yyyy-MM-dd''T''HH:mm:ss');
        realtime = time:minutes(1):enddate;
        spids_hour = any(hour(pricetime_dk) == [17, 18, 19, 20], 2);
        spids_mon = any(month(pricetime_dk) == [1, 2, 3, 10, 11, 12], 2);
        spids_mask = all([spids_hour, spids_mon], 2);
        lav_mask = ~spids_mask;
        % The price is in DKK/MWh
        price = price.SpotPriceDKK / 1000; % DKK/kWh
        % Tarrifs (see https://n1.dk/priser-og-vilkaar)
        price(lav_mask) = price(lav_mask) + 0.3268;
        price(spids_mask) = price(spids_mask) + 0.8409;
    end
    if (k > kn)
        k = 1;
        timehour = hour(time);
        t0 = time; t1 = time;
        t0.Hour = 13;
        t0.Minute = 0;
        t0.Second = 0;
        t1.Hour = 23;
        t1.Minute = 0;
        t1.Second = 0;
        if timehour < 13
            t0 = t0 - hours(24);
        else
            t1 = t1 + hours(24);
        end
        kn = ceil(minutes(t0 + hours(24) - time));
        sm = all([t0 <= pricetime_dk, pricetime_dk <= t1], 2);
        price_i = price(sm);
        price_i = [price_i; price_i(1)];
        pricetime_utc_i = pricetime_utc(sm);
        pricetime_utc_i = [pricetime_utc_i; pricetime_utc_i(1) + hours(1)];
        w = interp1(pricetime_utc_i, price_i, realtime, "previous")';
    end
    d = w(k:end);
    k = k + 1;
end

