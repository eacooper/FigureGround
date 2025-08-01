function out = compute_mean_and_quantiles(data)
%
% compute the mean, std, 25th, 50th, and 75th quantiles for some data

% out.m = mean(data(~isnan(data)));
% out.std = std(data(~isnan(data)));
% out.CI = 1.96*std(data(~isnan(data)))./sqrt(numel(data(~isnan(data))));
% out.q1= quantile(data(~isnan(data)),0.25);
% out.q2= quantile(data(~isnan(data)),0.50);
% out.q3 = quantile(data(~isnan(data)),0.75);

out.m = mean(data, 'omitnan');
out.std = std(data, 'omitnan');
out.CI = 1.96*std(data, 'omitnan')./sqrt(numel(data(~isnan(data))));

out.q1 = [];
out.q2 = [];
out.q3 = [];

for n = 1:size(data,2)

    out.q1(n)= quantile(data(~isnan(data(:, n)), n),0.25);
    out.q2(n)= quantile(data(~isnan(data(:, n)), n),0.50);
    out.q3(n)= quantile(data(~isnan(data(:, n)), n),0.75);

end


