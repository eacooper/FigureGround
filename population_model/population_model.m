clear all; close all;
addpath(genpath('.'));
rng(1986);

% Parameters for the Gaussian tuning curves
num_neurons = 50;           % Number of neurons in the population
sigma = 1;                  % Standard deviation in log10 space
A = 1;                      % Amplitude

% Tile the preferred speeds (mu_values) from log10(0.0001) to log10(10000)
mu_values = linspace(log10(0.0001), log10(10000), num_neurons);

% Define the Gaussian tuning curve function (log10 of speed as input)
gaussian_tuning = @(x, mu, sigma, A) A * exp(-(x - mu).^2 / (2 * sigma^2));

% Define the example figure speeds in deg/sec
figure_speeds = logspace(log10(0.5),log10(8),8);
figure_speeds_log = log10(figure_speeds);

% Define the example ground speeds, as a scale factor on figure speeds
ground_scales = rand(1,8)*0.5 + 0.3;
ground_speeds_log = log10(ground_scales.*figure_speeds);

% Preallocate vectors to store decoded speeds (in log10 space)
decoded_speeds_prioritizefaster = zeros(size(figure_speeds_log));
decoded_speeds_average = zeros(size(figure_speeds_log));
decoded_speeds_prioritizeslower = zeros(size(figure_speeds_log));

% Create output folder for EPS files
output_folder = 'population_model_plots';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end


% Loop through each figure speed and decode it
for i = 1:length(figure_speeds_log)

    f_log = figure_speeds_log(i);  % Current figure speed in log10 space
    g_log = ground_speeds_log(i);  % Current ground speed in log10 space

    % Compute the population response for the faster stimulus
    faster_log = max([f_log g_log]);
    responses_faster = gaussian_tuning(faster_log, mu_values, sigma, A);
    decoded_speeds_prioritizefaster(i) = sum(mu_values .* responses_faster) / sum(responses_faster);

    % Compute the population response for the slower stimulus
    slower_log = min([f_log g_log]);
    responses_slower = gaussian_tuning(slower_log, mu_values, sigma, A);
    decoded_speeds_prioritizeslower(i) = sum(mu_values .* responses_slower) / sum(responses_slower);

    % Compute the population response for the average
    responses_average = mean([responses_faster ; responses_slower]);
    decoded_speeds_average(i) = sum(mu_values .* responses_average) / sum(responses_average);

    if i == 3
        % Create figures for the population responses
        figure; setupfig(2,0.5,10); hold on;
        bar(mu_values, responses_faster, 'FaceColor', [0 0.5 0]);
        set(gca, 'XTick', [], 'YTick', []);
        plot(decoded_speeds_prioritizefaster(i),1.1*max(responses_faster),'ko');
        ylim([0 1.2*max(responses_faster)]);

        % Save the figure
        saveas(gcf, fullfile(output_folder, 'population_responses_faster.eps'), 'epsc');

        figure; setupfig(2,0.5,10); hold on;
        bar(mu_values, responses_average, 'FaceColor', [0.5 0.5 0]);
        set(gca, 'XTick', [], 'YTick', []);
        plot(decoded_speeds_average(i),1.1*max(responses_average),'ko');
        ylim([0 1.2*max(responses_average)]);

        % Save the figure
        saveas(gcf, fullfile(output_folder, 'population_responses_avg.eps'), 'epsc');

        figure; setupfig(2,0.5,10); hold on;
        bar(mu_values, responses_slower, 'FaceColor', [0.5 0 0.5]);
        set(gca, 'XTick', [], 'YTick', []);
        plot(decoded_speeds_prioritizeslower(i),1.1*max(responses_slower),'ko');
        ylim([0 1.2*max(responses_slower)]);

        % Save the figure
        saveas(gcf, fullfile(output_folder, 'population_responses_slow.eps'), 'epsc');

    end
end


% Plot the decoded results with log-log axes and custom tick labels
figure; setupfig(2,1.2,5); 
%set(gca, 'FontSize', 12 * 2);

% Plot decoded speeds with log-log scaling
plot(figure_speeds, 10.^decoded_speeds_prioritizefaster, 'o-', 'LineWidth', 1, 'Color', [0 0.5 0], 'DisplayName', 'Faster','markersize',2);  % Green
hold on;
plot(figure_speeds, 10.^decoded_speeds_average, 'o-', 'LineWidth', 1, 'Color', [0.5 0.5 0], 'DisplayName', 'Average','markersize',2);  % Olive
plot(figure_speeds, 10.^decoded_speeds_prioritizeslower, 'o-', 'LineWidth', 1, 'Color', [0.5 0 0.5], 'DisplayName', 'Slower','markersize',2);  % Purple
plot([0.2, 10], [0.2, 10], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Identity Line');  % Identity line

% Set log-log axes, custom limits, and custom tick labels
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([0.2, 10]);
ylim([0.2, 10]);

% Set custom tick labels to match figure speeds
set(gca, 'XTick', [0.5, 1, 2, 4, 8], 'YTick', [0.5, 1, 2, 4, 8]);

% Save the figure
saveas(gcf, fullfile(output_folder, 'decoded_results.eps'), 'epsc');



% Plot all tuning curves with restricted x-axis limits
x_log = linspace(log10(0.0001), log10(10000), 100);
tuning_curves = zeros(num_neurons, length(x_log));
for i = 1:num_neurons
    tuning_curves(i, :) = gaussian_tuning(x_log, mu_values(i), sigma, A);
end

figure;
plot(10.^x_log, tuning_curves', 'LineWidth', 1.2);
set(gca, 'XScale', 'log', 'FontSize', 12 * 2);
xlim([10^-2, 10^2]);
set(gca, 'XTick', [], 'YTick', []);  % Remove x and y tick labels

% Save the figure
saveas(gcf, fullfile(output_folder, 'tuning_curves.eps'), 'epsc');

% Generate and plot overlapping histograms with filled red and blue bars
figure;
data_red = normrnd(6, 1, [1, 1000]);  % Red distribution with higher mean
data_blue = normrnd(4, 1, [1, 1000]); % Blue distribution with lower mean
histogram(data_red, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 1);  % Opaque red bars
hold on;
histogram(data_blue, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 1);  % Opaque blue bars

% Remove axis labels and ticks
set(gca, 'XTick', [], 'YTick', [], 'FontSize', 12 * 2);

% Save the histogram plot
saveas(gcf, fullfile(output_folder, 'speed_histograms.eps'), 'epsc');

% Display decoded results
disp('Figure Speed (deg/sec)   Decoded (Original)   Decoded (75%)   Decoded (50%)');
disp([figure_speeds(:), 10.^decoded_speeds_prioritizefaster(:), 10.^decoded_speeds_average(:), 10.^decoded_speeds_prioritizeslower(:)]);
