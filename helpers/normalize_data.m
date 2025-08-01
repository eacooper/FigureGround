function x_normed = normalize_data(x,lowerPerc,upperPerc)
%
% take in a vector of data and normalize range to be 0-1, robust to
% outliers

x_min            = quantile(x(:),lowerPerc);
x_max            = quantile(x(:),upperPerc);
x_normed      = (x - x_min) ./ (x_max-x_min);
x_normed(x_normed < 0) = NaN;
x_normed(x_normed > 1) = NaN;

end