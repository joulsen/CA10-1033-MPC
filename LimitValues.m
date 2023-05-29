function values = LimitValues(values, min_lim, max_lim)

values(values < min_lim) = min_lim;
values(values > max_lim) = max_lim;
