% Plots and analysis of simulated receptive fields
close all;

% use this flag to toggle whether analyzing world motion and simulated
% pursuit motion
%fixes_flag = 'stationary';
fixes_flag = 'pursuit';

% set path and load data
save_dir = fullfile(resultsPath, 'motion_RF_analysis');
load(fullfile(save_dir,'RF_motion_measurements.mat'));
display('loaded RF_data');

% if we're analyzing pursuit, grab the appropriate variables, otherwise the
% default is to use the stationary gaze simulation
if strcmp(fixes_flag,'pursuit')

    RF_hist_mag     = RF_hist_mag_pursuit;
    RF_hist_mag_log = RF_hist_mag_log_pursuit;
    RF_hist_ori     = RF_hist_ori_pursuit;

    RF_hist_mag_log_fig = RF_hist_mag_log_fig_pursuit;
    RF_hist_ori_fig = RF_hist_ori_fig_pursuit;

    RF_hist_mag_log_bg = RF_hist_mag_log_bg_pursuit;
    RF_hist_ori_bg = RF_hist_ori_bg_pursuit;

    RF_mag_mean_fig = RF_mag_mean_fig_pursuit;
    RF_ori_mean_fig = RF_ori_mean_fig_pursuit;
    RF_ori_var_fig = RF_ori_var_fig_pursuit;

    RF_mag_mean_bg = RF_mag_mean_bg_pursuit;
    RF_ori_mean_bg = RF_ori_mean_bg_pursuit;
    RF_ori_var_bg = RF_ori_var_bg_pursuit;

    RF_mag_mean_all = RF_mag_mean_all_pursuit;
    RF_ori_mean_all = RF_ori_mean_all_pursuit;
    RF_ori_var_all = RF_ori_var_all_pursuit;

end

% we'll analyze all data together and by eccentricity
cond = {'all','2pt5_deg','5_deg','10_deg','15_deg'};

% for each eccentricity condition
for pl = 1:numel(cond)

    display(['analyzing ' cond{pl} '...']);

    % get correct indices
    switch cond{pl}

        % get indices of all desired RFs (R)
        case 'all'
            sim_indsFG = 1:numel(RF_diameters);
        case '2pt5_deg'
            sim_indsFG = find(RF_diameters == RF_sizes(1));
        case '5_deg'
            sim_indsFG = find(RF_diameters == RF_sizes(2));
        case '10_deg'
            sim_indsFG = find(RF_diameters == RF_sizes(3));
        case '15_deg'
            sim_indsFG = find(RF_diameters == RF_sizes(4));
    end


    %% basic description of RFs and fixations

    if strcmp(cond{pl},'all')

        % report number of RFs
        display(['total  RFs = ' num2str(numel(sim_indsFG))]);

        % mean percent of each RF that is fg and bg
        display(['-- mean perc fig = ' num2str(mean(RF_perc_fig))]);
        display(['-- mean perc bg = ' num2str(mean(RF_perc_bg))]);

        % percent of fixations that fall on figures
        these_fixation_figs = RF_sim_fixation_fig(sim_indsFG);
        display(['-- perc of fixations that fall on figures ' num2str(100*sum(these_fixation_figs)/numel(these_fixation_figs))]);

        % plot basic statistics of the simulation
        figure; hold on; setupfig(12,12,12);
        subplot(2,2,1); hold on; title('scene indices'); histogram(RF_scenes(sim_indsFG)); box on;
        subplot(2,2,2); hold on; title('RF diameters in deg'); histogram(RF_diameters(sim_indsFG),0:20); box on;
        subplot(2,2,3); hold on; title('percent fig/bg')
        hf(1) = histogram(RF_perc_fig(sim_indsFG)); box on;
        hf(2) = histogram(RF_perc_bg(sim_indsFG));
        legend(hf,'perc fig','perc bg')
        subplot(2,2,4); hold on; title('RF locations'); scatter(RF_center_pix(sim_indsFG,2),RF_center_pix(sim_indsFG,1),'r');
        xlabel('pixel location x'); ylabel('pixel location y');

        saveas(gcf,fullfile(save_dir,['RF_stats_fixations-' fixes_flag '.png']));

        % smooth pursuit stats for FG RFs
        if strcmp(fixes_flag,'pursuit')

            % grab relevant data
            these_pursuit_mags = RF_sim_fixation_mag(sim_indsFG);
            these_pursuit_oris = RF_sim_fixation_ori(sim_indsFG);
            these_fixations = RF_sim_fixation_pix(sim_indsFG,:);

            % pursuit characteristics
            display(['mean / median speed pursuit = ' num2str(mean(these_pursuit_mags)) ' | ' num2str(median(these_pursuit_mags))]);
            display(['min / max speed pursuit = ' num2str(min(these_pursuit_mags)) ' | ' num2str(max(these_pursuit_mags))]);
            display(['prop zero speed pursuit = ' num2str(sum(these_pursuit_mags == 0)/numel(these_pursuit_mags))]);

            figure; hold on; setupfig(8,5,12);
            subplot(1,3,1); hold on; title('pursuit speeds'); histogram(these_pursuit_mags); box on;
            subplot(1,3,2); hold on; title('pursuit directions'); ph(1) = histogram(these_pursuit_oris); box on;
            subplot(1,3,3); hold on; title('fixation locations'); scatter(these_fixations(:,2),these_fixations(:,1),'k','filled');
            xlabel('pixel location x'); ylabel('pixel location y');

            saveas(gcf,fullfile(save_dir,'RF_pursuit_stats.png'));

        end

    end


    %% compute histogram of speed difference between figure and ground in border RFs

    FG_diff_mag{pl}     = RF_mag_mean_fig(sim_indsFG) - RF_mag_mean_bg(sim_indsFG);
    FG_diff_mag_log{pl} = log10(RF_mag_mean_fig(sim_indsFG)) - log10(RF_mag_mean_bg(sim_indsFG));

    % prop diff and binomial CI
    [phat_mag_log_prop_faster(pl),pci_mag_log_prop_faster(pl,:)] = binofit(sum(FG_diff_mag_log{pl} > 0),numel(FG_diff_mag_log{pl}));


    % compute means
    mean_diff_mag(pl) = mean(FG_diff_mag{pl});
    mean_diff_mag_log(pl) = mean(FG_diff_mag_log{pl});

    % t-test
    [~,p_mag_log_FG_diff(pl),~,stats_tmp] = ttest(FG_diff_mag_log{pl},0,'Alpha',0.05,'tail','both');
    stats_mag_log_FG_diff.tstat(pl) = stats_tmp.tstat;
    stats_mag_log_FG_diff.df(pl) = stats_tmp.df;

    % cohen's D
    eff_size_mag_log_FG_diff(pl) = mean_diff_mag_log(pl) / std(FG_diff_mag_log{pl});


    %% circular variance difference between fig and ground

    FG_diff_circ_var{pl}  = RF_ori_var_fig(sim_indsFG) - RF_ori_var_bg(sim_indsFG);

    % proportion more coherent
    [phat_circ_var_prop_fig_less(pl),pci_circ_var_prop_fig_less(pl,:)] = binofit(sum(FG_diff_circ_var{pl} < 0),numel(FG_diff_circ_var{pl}));

    % compute means
    mean_diff_circ_var(pl) = mean(FG_diff_circ_var{pl});

    % t-test
    [~,p_circ_var_FG_diff(pl),~,stats_tmp] = ttest(FG_diff_circ_var{pl},0,'Alpha',0.05,'tail','both');
    stats_circ_var_FG_diff.tstat(pl) = stats_tmp.tstat;
    stats_circ_var_FG_diff.df(pl) = stats_tmp.df;

    % cohen's D
    eff_size_circ_var_FG_diff(pl) = mean_diff_circ_var(pl) / std(FG_diff_circ_var{pl});


    %% average over the histograms -

    % get appropriate values

    RF_hist_mag_logFG0      = RF_hist_mag_log(sim_indsFG,:);
    RF_hist_oriFG0          = RF_hist_ori(sim_indsFG,:);

    RF_hist_mag_logFG_fig0  = RF_hist_mag_log_fig(sim_indsFG,:);
    RF_hist_oriFG_fig0      = RF_hist_ori_fig(sim_indsFG,:);

    RF_hist_mag_logFG_bg0   = RF_hist_mag_log_bg(sim_indsFG,:);
    RF_hist_oriFG_bg0       = RF_hist_ori_bg(sim_indsFG,:);

    % convert raw counts to probability density

    % totals
    total_mag_logFG0 = sum(RF_hist_mag_logFG0,2);
    total_oriFG0 = sum(RF_hist_oriFG0,2);

    total_mag_logFG_fig0 = sum(RF_hist_mag_logFG_fig0,2);
    total_oriFG_fig0 = sum(RF_hist_oriFG_fig0,2);

    total_mag_logFG_bg0 = sum(RF_hist_mag_logFG_bg0,2);
    total_oriFG_bg0 = sum(RF_hist_oriFG_bg0,2);

    % norm
    RF_pd_mag_logFG0 = RF_hist_mag_logFG0./total_mag_logFG0;
    RF_pd_oriFG0 = RF_hist_oriFG0./total_oriFG0;

    RF_pd_mag_logFG_fig0 = RF_hist_mag_logFG_fig0./total_mag_logFG_fig0;
    RF_pd_oriFG_fig0 = RF_hist_oriFG_fig0./total_oriFG_fig0;

    RF_pd_mag_logFG_bg0 = RF_hist_mag_logFG_bg0./total_mag_logFG_bg0;
    RF_pd_oriFG_bg0 = RF_hist_oriFG_bg0./total_oriFG_bg0;

    % compute averages and ranges - raw histograms

    % FG RFs, all points
    mag_logFG0   = compute_mean_and_quantiles(RF_hist_mag_logFG0); % normed mag, whole RF - fig/bg RFs
    oriFG0   = compute_mean_and_quantiles(RF_hist_oriFG0); % normed ori, whole RF - fig/bg RFs

    % FG RFs, figure regions
    mag_logFG_fig0   = compute_mean_and_quantiles(RF_hist_mag_logFG_fig0); % normed mag, figure region - fig/bg RFs
    oriFG_fig0   = compute_mean_and_quantiles(RF_hist_oriFG_fig0); % normed ori, figure region - fig/bg RFs

    % FG RFs, ground regions
    mag_logFG_bg0   = compute_mean_and_quantiles(RF_hist_mag_logFG_bg0); % normed mag, background region - fig/bg RFs
    oriFG_bg0   = compute_mean_and_quantiles(RF_hist_oriFG_bg0); % normed ori, background region - fig/bg RFs

    % averages - probability density for ratios

    % FG RFs, figure regions
    mag_logFG_fig0_pd   = compute_mean_and_quantiles(RF_pd_mag_logFG_fig0); % normed mag, figure region - fig/bg RFs
    oriFG_fig0_pd   = compute_mean_and_quantiles(RF_pd_oriFG_fig0); % normed ori, figure region - fig/bg RFs

    % FG RFs, ground regions
    mag_logFG_bg0_pd   = compute_mean_and_quantiles(RF_pd_mag_logFG_bg0); % normed mag, background region - fig/bg RFs
    oriFG_bg0_pd   = compute_mean_and_quantiles(RF_pd_oriFG_bg0); % normed ori, background region - fig/bg RFs

    this_line_width = 2;

    if strcmp(cond{pl},'all')

        % Figure and ground probabilities:

        figure; hold on; setupfig(4,4,12)
        a = area(cntr_mag,mag_logFG0.m);
        a.FaceColor = [0.75 0.75 0.75];
        a.EdgeColor = [0.75 0.75 0.75];
        hf(1) = plot(cntr_mag,mag_logFG_fig0.m,'-', 'Color',figColour,'LineWidth', this_line_width);
        plot(cntr_mag, [mag_logFG_fig0.m - mag_logFG_fig0.CI; mag_logFG_fig0.m + mag_logFG_fig0.CI], ':', 'Color',figColour, 'LineWidth', this_line_width);
        hf(2) = plot(cntr_mag, mag_logFG_bg0.m,'-', 'Color',bgColour,'LineWidth', this_line_width);
        plot(cntr_mag, [mag_logFG_bg0.m - mag_logFG_bg0.CI; mag_logFG_bg0.m + mag_logFG_bg0.CI], ':', 'Color',bgColour, 'LineWidth', this_line_width);
        xlabel('speed (norm)'); axis square; box on;
        ylabel('frequency'); axis square; box on;
        ylim([0 3000]);
        legend(hf,'figure','ground');
        saveas(gcf,fullfile(save_dir,['RF_fg_speed_probs-' fixes_flag '.png']));

        figure; setupfig(4,4,12)
        polarplot(edges_ori, [oriFG0.m, oriFG0.m(1)],'-', 'Color',[0.75 0.75 0.75],'LineWidth', this_line_width); hold on;
        polarplot(edges_ori, [oriFG_fig0.m, oriFG_fig0.m(1)],'-', 'Color',figColour,'LineWidth', this_line_width); hold on;
        polarplot(edges_ori, [[oriFG_fig0.m - oriFG_fig0.CI, oriFG_fig0.m(1) - oriFG_fig0.CI(1)]; [oriFG_fig0.m + oriFG_fig0.CI, oriFG_fig0.m(1) + oriFG_fig0.CI(1)]], ':', 'Color',figColour, 'LineWidth', this_line_width);
        polarplot(edges_ori, [oriFG_bg0.m, oriFG_bg0.m(1)],'-', 'Color',bgColour,'LineWidth', this_line_width);
        polarplot(edges_ori, [[oriFG_bg0.m - oriFG_bg0.CI, oriFG_bg0.m(1) - oriFG_bg0.CI(1)]; [oriFG_bg0.m + oriFG_bg0.CI, oriFG_bg0.m(1) + oriFG_bg0.CI(1)]], ':', 'Color',bgColour, 'LineWidth', this_line_width);
        rlim([0 10000]); rticks([]);
        saveas(gcf,fullfile(save_dir,['RF_fg_direction_probs-' fixes_flag '.png']));
    end

    % figure probability ratio
    ratio_mag_log(pl,:) = mag_logFG_fig0_pd.m./mag_logFG_bg0_pd.m;
    ratio_ori(pl,:)     = oriFG_fig0_pd.m./oriFG_bg0_pd.m;

    % bootstrap the 95% CI for these ratios
    doboot = @(f,g) (mean(f, 'omitnan'))./(mean(g, 'omitnan'));
    ratio_mag_log_ci(:,:,pl) = bootci(1000,doboot,RF_pd_mag_logFG_fig0,RF_pd_mag_logFG_bg0);
    ratio_ori_ci(:,:,pl) = bootci(1000,doboot,RF_pd_oriFG_fig0,RF_pd_oriFG_bg0);

end

%% plot histogram of speed difference between figure and ground

all_ind = find(strcmp(cond,'all'));

figure; hold on; setupfig(4,4,12)
histogram(FG_diff_mag_log{all_ind},linspace(-2,2,21),'facecolor',[0.75 0.75 0.75],'facealpha',1)
xlabel('speed difference between figure and ground (log deg/s)'); axis square; box on;
ylabel('frequency')
saveas(gcf,fullfile(save_dir,'RF_fg_speed_differences.png'));

% each eccentricity separately
ecc_inds = find(~strcmp(cond,'all'));
ecc_dark = [0.9 0.7 0.5 0.3]; % color mapping

% each eccentricity separately - line plot
figure; hold on; setupfig(4,4,12);

for ex = 1:numel(ecc_inds)
    % Extract data for the current eccentricity
    this_ecc = ecc_inds(ex);
    data = FG_diff_mag_log{this_ecc};

    % Compute mean and 95% confidence interval (CI)
    mean_val = mean(data);
    std_err = std(data) / sqrt(length(data)); % Standard Error
    ci_95 = 1.96 * std_err; % 95% CI using z-score for normal distribution

    % Plot mean as a point and CI as error bars
    errorbar(ex, mean_val, ci_95, 'o-', ...
        'Color', repmat(ecc_dark(ex), 1, 3), ...
        'LineWidth', 1.5, 'MarkerSize', 8,'MarkerFaceColor', repmat(ecc_dark(ex), 1, 3)); hold on;
end

% Label axes and set plot properties
xlabel('Eccentricity (deg)');
ylabel('Mean speed difference (log deg/s)');
axis square; box on;
xlim([0.5 4.5]); ylim([-0.1 0.6]);
plot([0.5 4.5],[0 0],'k-')
xticks(1:numel(ecc_inds));
xticklabels({'2.5','5','10','15'});

if strcmp(fixes_flag,'stationary')
    ylim([-0.1 1.2])
end

saveas(gcf,fullfile(save_dir,['RF_fg_speed_differences_by_ecc-' fixes_flag '.png']));

% speed difference stats

T_DIFF_MAG = table(cond',phat_mag_log_prop_faster',pci_mag_log_prop_faster(:,1),pci_mag_log_prop_faster(:,2),...
    mean_diff_mag',mean_diff_mag_log',stats_mag_log_FG_diff.tstat',stats_mag_log_FG_diff.df',...
    eff_size_mag_log_FG_diff',p_mag_log_FG_diff',...
    'VariableNames',{'ecc','prop_fig_faster','prop_fig_faster_LCI','prop_fig_faster_HCI',...
    'mean_Diff_degsec','mean_Diff_LOGdegsec','tstat','df',...
    'cohensD','pval'});

display(T_DIFF_MAG);

save(fullfile(save_dir,['RFs_stats_speed_diffs-' fixes_flag '.mat']),'T_DIFF_MAG');


% run a speed difference anova on Eccentricity
% data appear roughly normal with good homogeneity of variances

% Initialize containers for combined data and group labels
allData = [];
group = [];

% Loop through each group to extract data and create group labels
for ex = 1:numel(ecc_inds)
    this_ecc    = ecc_inds(ex);
    groupData   = FG_diff_mag_log{this_ecc}';          % Data for the current group
    allData     = [allData; groupData];  % Append the group's data to allData
    group       = [group; repmat(cond(this_ecc), numel(groupData), 1)];  % Create group labels
end

% Run one-way ANOVA
[p_mag_log, tbl_mag_log, stats_mag_log] = anova1(allData, group);

% Display the ANOVA table
display('Ecc ANOVA for Speed diff')
disp(tbl_mag_log);

SS_between = cell2mat(tbl_mag_log(2, 2));  % Sum of Squares between groups (SSB)
SS_total = cell2mat(tbl_mag_log(4, 2));    % Total Sum of Squares (SST)
eta_squared = SS_between / SS_total;
fprintf('Effect size (η²): %.4f\n', eta_squared);

% Perform post-hoc multiple comparisons
[comparison_mag_log,means_mag_log,h_mag_log,gnames_mag_log] = multcompare(stats_mag_log, "CriticalValueType","hsd");

resultsTable = create_eccentricity_comparison_table(stats_mag_log, comparison_mag_log, 'Mean Diff (log deg/sec)',allData,group);

% Display the Table
disp(resultsTable);

% Save the Table to CSV
writetable(resultsTable, fullfile(save_dir,['RFs_stats_speed_diffs_pairwise-' fixes_flag '.csv']));


%% coherence difference

figure; hold on; setupfig(4,4,12)
histogram(FG_diff_circ_var{all_ind},linspace(-1,1,21),'facecolor',[0.75 0.75 0.75],'facealpha',1)
xlabel('coherence difference between figure and ground'); axis square; box on;
ylabel('freqency'); set(gca,'xdir','reverse')
saveas(gcf,fullfile(save_dir,'RF_fg_orivar_differences.png'));


% each eccentricity separately - line plot
figure; hold on; setupfig(4,4,12);

for ex = 1:numel(ecc_inds)
    % Extract data for the current eccentricity
    this_ecc = ecc_inds(ex);
    data = FG_diff_circ_var{this_ecc};

    % Compute mean and 95% confidence interval (CI)
    mean_val = mean(data);
    std_err = std(data) / sqrt(length(data)); % Standard Error
    ci_95 = 1.96 * std_err; % 95% CI using z-score for normal distribution

    % Plot mean as a point and CI as error bars
    errorbar(ex, mean_val, ci_95, 'o-', ...
        'Color', repmat(ecc_dark(ex), 1, 3), ...
        'LineWidth', 1.5, 'MarkerSize', 8,'MarkerFaceColor', repmat(ecc_dark(ex), 1, 3)); hold on;
end

% Label axes and set plot properties
xlabel('Eccentricity (deg)');
ylabel('Mean variance difference');
axis square; box on;
set(gca,'ydir','reverse')
xlim([0.5 4.5]); ylim([-0.15 0.05]);
plot([0.5 4.5],[0 0],'k-')
xticks(1:numel(ecc_inds));
xticklabels({'2.5','5','10','15'});

if strcmp(fixes_flag,'stationary')
    ylim([-0.2 0.05])
end

saveas(gcf,fullfile(save_dir,['RF_fg_orivar_differences_by_ecc-' fixes_flag '.png']));

% coherence difference stats

T_DIFF_CIRC_VAR = table(cond',phat_circ_var_prop_fig_less',pci_circ_var_prop_fig_less(:,1),pci_circ_var_prop_fig_less(:,2),...
    mean_diff_circ_var',stats_circ_var_FG_diff.tstat',stats_circ_var_FG_diff.df',...
    eff_size_circ_var_FG_diff',p_circ_var_FG_diff',...
    'VariableNames',{'ecc','prop_fig_more_coh','prop_fig_more_coh_LCI','prop_fig_more_coh_HCI',...
    'mean_Diff_circ_var','tstat','df',...
    'cohensD','pval'});


display(T_DIFF_CIRC_VAR);

save(fullfile(save_dir,['RFs_stats_orivar_diffs-' fixes_flag '.mat']),'T_DIFF_CIRC_VAR');

% run an orientation difference anova on Eccentricity
% data appear roughly normal with good homogeneity of variances

% Initialize containers for combined data and group labels
allData = [];
group = [];

% Loop through each group to extract data and create group labels
for ex = 1:numel(ecc_inds)
    this_ecc    = ecc_inds(ex);
    groupData   = FG_diff_circ_var{this_ecc}';          % Data for the current group
    allData     = [allData; groupData];  % Append the group's data to allData
    group       = [group; repmat(cond(this_ecc), numel(groupData), 1)];  % Create group labels
end

% Run one-way ANOVA
[p_circ_var, tbl_circ_var, stats_circ_var] = anova1(allData, group);

% Display the ANOVA table
display('Ecc ANOVA for Circ Var diff')
disp(tbl_circ_var);

SS_between = cell2mat(tbl_circ_var(2, 2));  % Sum of Squares between groups (SSB)
SS_total = cell2mat(tbl_circ_var(4, 2));    % Total Sum of Squares (SST)
eta_squared = SS_between / SS_total;
fprintf('Effect size (η²): %.4f\n', eta_squared);

% Perform post-hoc multiple comparisons
[comparison_circ_var,means_circ_var,h_circ_var,gnames_circ_var] = multcompare(stats_circ_var, "CriticalValueType","hsd");

resultsTable = create_eccentricity_comparison_table(stats_circ_var, comparison_circ_var, 'Mean Diff (deg2)',allData,group);

% Display the Table
disp(resultsTable);

% Save the Table to CSV
writetable(resultsTable, fullfile(save_dir,['RFs_stats_orivar_diffs_pairwise-' fixes_flag '.csv']));


% figure/ground probability ratio

% all ecc

figure; hold on; setupfig(4,4,12)
plot(cntr_mag,ratio_mag_log(all_ind,:),'k-', 'LineWidth', this_line_width)
plot(cntr_mag,ratio_mag_log_ci(1,:,all_ind),'k:', 'LineWidth', this_line_width)
plot(cntr_mag,ratio_mag_log_ci(2,:,all_ind),'k:','LineWidth', this_line_width)
plot([cntr_mag(1) cntr_mag(end)],[1 1],'k-')
ax = gca;
ax.Color = 'none';
xlabel('speed (norm)');
ylabel('figure/ground probability'); axis square; box on;
saveas(gcf,fullfile(save_dir,['RF_speed_probability_ratios-' fixes_flag '.png']));

figure; setupfig(4,4,12)
polarplot(edges_ori, [ratio_ori(all_ind,:), ratio_ori(all_ind,1)],'k-','LineWidth', this_line_width); hold on;
polarplot(edges_ori, [ratio_ori_ci(1,:,all_ind), ratio_ori_ci(1,1,all_ind)],'k:', 'LineWidth', this_line_width); hold on;
polarplot(edges_ori, [ratio_ori_ci(2,:,all_ind), ratio_ori_ci(2,1,all_ind)],'k:','LineWidth', this_line_width); hold on;
rticks([0.5 1 1.5]);
saveas(gcf,fullfile(save_dir,['RF_direction_probability_ratios-' fixes_flag '.png']));
