function d = get_price_historical(limit, startdate, enddate)
%GET_PRICE_HISTORICAL Summary of this function goes here
%   Detailed explanation goes here
    persistent w k
    if (isempty(w))
        k = 1;
        url = "https://api.energidataservice.dk/dataset/Elspotprices?offset=0&start=%sT00:00&end=%sT00:00%s";
        extra = "&filter=%7B%22PriceArea%22:[%22DK1%22]%7D&sort=HourUTC%20DESC&timezone=dk";
        startdate.Format = 'yyyy-MM-dd';
        enddate.Format = 'yyyy-MM-dd';
        price = webread(sprintf(url, startdate, enddate + days(1), extra), ...
                        weboptions("Timeout", 20));
        price = struct2table(price.records);
        pricetime_utc = datetime(price.HourUTC, "InputFormat", 'yyyy-MM-dd''T''HH:mm:ss');
        pricetime_dk = datetime(price.HourDK, "InputFormat", 'yyyy-MM-dd''T''HH:mm:ss');
        realtime = startdate:minutes(1):enddate;
        spids_hour = any(hour(pricetime_dk) == [17, 18, 19, 20], 2);
        spids_mon = any(month(pricetime_dk) == [1, 2, 3, 10, 11, 12], 2);
        spids_mask = all([spids_hour, spids_mon], 2);
        lav_mask = ~spids_mask;
        % The price is in DKK/MWh
        price = price.SpotPriceDKK / 1000; % DKK/kWh
        % Tarrifs (see https://n1.dk/priser-og-vilkaar)
        price(lav_mask) = price(lav_mask) + 0.3268;
        price(spids_mask) = price(spids_mask) + 0.8409;
        w = interp1(pricetime_utc, price, realtime, "previous")';
    end
    if limit == -1
        limit = size(w,1);
    else
        limit = min(k+limit, size(w,1));
    end
    d = w(k:limit);
    k = k + 1;
end

