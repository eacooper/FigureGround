% Plots and analysis of simulated receptive fields
close all;

% set path and load data
save_dir = fullfile(resultsPath, 'distance_RF_analysis');
load(fullfile(save_dir,'RF_distance_measurements.mat'));
display('loaded RF_data');

% flip Dioptric values left/right so that 0 = near and 1 = far
RF_hist_Ds = fliplr(RF_hist_Ds);
RF_hist_Ds_fig = fliplr(RF_hist_Ds_fig);
RF_hist_Ds_bg = fliplr(RF_hist_Ds_bg);

% we'll analyze all data together and by eccentricity
cond = {'all','2pt5_deg','5_deg','10_deg','15_deg'};

% for each condition
for pl = 1:numel(cond)

    display(['analyzing ' cond{pl} '...']);

    % get correct indices
    switch cond{pl}

        % get indices of all desired RFs (R)
        % and the subset for which there'sat least 25% figure and 25% background (FG)
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

    %% basic description

    if strcmp(cond{pl},'all')

        % report number of RFs and number of F/G RFs
        display(['total RFs = ' num2str(numel(sim_indsFG))]);

        % mean percent fg and bg
        display(['-- mean perc fig = ' num2str(mean(RF_perc_fig))]);
        display(['-- mean perc bg = ' num2str(mean(RF_perc_bg))]);

        display(['mean / median fixation distance = ' num2str(mean(RF_sim_fixation(sim_indsFG,3))) ' | ' num2str(median(RF_sim_fixation(sim_indsFG,3)))]);

        % plot basic statistics of the simulation
        figure; hold on; setupfig(8,8,12);
        subplot(2,2,1); hold on; title('scene indices'); histogram(RF_scenes(sim_indsFG)); box on;
        subplot(2,2,2); hold on; title('RF diameters in deg'); histogram(RF_diameters(sim_indsFG),0:20); box on;
        subplot(2,2,3); hold on; title('percent fig/bg')
        hf(1) = histogram(RF_perc_fig(sim_indsFG)); box on;
        hf(2) = histogram(RF_perc_bg(sim_indsFG));
        legend(hf,'perc fig','perc_bg')
        subplot(2,2,4); hold on; title('RF locations'); scatter(RF_center_pix(sim_indsFG,2),RF_center_pix(sim_indsFG,1),'r');
        xlabel('pixel location x'); ylabel('pixel location y');

        saveas(gcf,fullfile(save_dir,'RF_stats.png'));

        figure; hold on; setupfig(8,5,12); sgtitle(cond{pl})
        subplot(1,2,1); hold on; title('simulated fixation point distance'); histogram(RF_sim_fixation(sim_indsFG,3)); box on;
        subplot(1,2,2); hold on; title('fixation locations'); scatter(RF_sim_fixation_r_c(:,2),RF_sim_fixation_r_c(:,1),'k','filled');
        xlabel('pixel location x'); ylabel('pixel location y');

        saveas(gcf,fullfile(save_dir,'RF_fixation_distance_stats.png'));

    end

    %% compute histogram of distance difference between figure and ground in border RFs
    FG_diff_Ds{pl}     = (1./RF_dist_mean_fig(sim_indsFG)) - (1./RF_dist_mean_bg(sim_indsFG));
    %FG_diff_meters{pl} = RF_dist_mean_fig(sim_indsFG) - RF_dist_mean_bg(sim_indsFG);

    % disparity - flip sign so that figure nearer = positive
    FG_diff_disparity{pl} = (-RF_disparity_mean_fig(sim_indsFG)) - (-RF_disparity_mean_bg(sim_indsFG));

    % proportion closer
    [phat_prop_fig_closer(pl),pci_prop_fig_closer(pl,:)] = binofit(sum(FG_diff_disparity{pl} > 0),numel(FG_diff_disparity{pl}));

    % compute means
    mean_diff_Ds(pl) = mean(FG_diff_Ds{pl});
    mean_diff_disparity(pl) = mean(FG_diff_disparity{pl});

    % rank sum test on disparity
    [~,p_FG_disparity_diff(pl),~,stats_tmp] = ttest(FG_diff_disparity{pl},0,'Alpha',0.05,'tail','both');
    stats_FG_disparity_diff.tstat(pl) = stats_tmp.tstat;
    stats_FG_disparity_diff.df(pl) = stats_tmp.df;

    % effect size
    % https://rpkgs.datanovia.com/rstatix/reference/wilcox_effsize.html
    eff_size_FG_disparity_diff(pl) = mean_diff_disparity(pl) / std(FG_diff_disparity{pl});


    %% average over the histograms

    % grab appropriate histograms

    % border RFs
    %RF_hist_distancesFG0 = RF_hist_distances(sim_indsFG,:);
    %RF_hist_DsFG0 = RF_hist_Ds(sim_indsFG,:);
    RF_hist_disparitiesFG0 = fliplr(RF_hist_disparities(sim_indsFG,:)); % disparity - flip sign so that figure nearer = positive

    % figures
    %RF_hist_distancesFG_fig0 = RF_hist_distances_fig(sim_indsFG,:);
    %RF_hist_DsFG_fig0 = RF_hist_Ds_fig(sim_indsFG,:);
    RF_hist_disparitiesFG_fig0 = fliplr(RF_hist_disparities_fig(sim_indsFG,:)); % disparity - flip sign so that figure nearer = positive

    % grounds
    %RF_hist_distancesFG_bg0 = RF_hist_distances_bg(sim_indsFG,:);
    %RF_hist_DsFG_bg0 = RF_hist_Ds_bg(sim_indsFG,:);
    RF_hist_disparitiesFG_bg0 = fliplr(RF_hist_disparities_bg(sim_indsFG,:)); % disparity - flip sign so that figure nearer = positive


    % convert raw counts to probability density

    % totals
    %total_distancesFG0 = sum(RF_hist_distancesFG0,2);
    %total_DsFG0 = sum(RF_hist_DsFG0,2);
    total_disparitiesFG0 = sum(RF_hist_disparitiesFG0,2);

    %total_distancesFG_fig0 = sum(RF_hist_distancesFG_fig0,2);
    %total_DsFG_fig0 = sum(RF_hist_DsFG_fig0,2);
    total_disparitiesFG_fig0 = sum(RF_hist_disparitiesFG_fig0,2);

    %total_distancesFG_bg0 = sum(RF_hist_distancesFG_bg0,2);
    %total_DsFG_bg0 = sum(RF_hist_DsFG_bg0,2);
    total_disparitiesFG_bg0 = sum(RF_hist_disparitiesFG_bg0,2);

    % norm

    %RF_pd_distancesFG0 = RF_hist_distancesFG0./total_distancesFG0;
    %RF_pd_DsFG0 = RF_hist_DsFG0./total_DsFG0;
    RF_pd_disparitiesFG0 = RF_hist_disparitiesFG0./total_disparitiesFG0;

    %RF_pd_distancesFG_fig0 = RF_hist_distancesFG_fig0./total_distancesFG_fig0;
    %RF_pd_DsFG_fig0 = RF_hist_DsFG_fig0./total_DsFG_fig0;
    RF_pd_disparitiesFG_fig0 = RF_hist_disparitiesFG_fig0./total_disparitiesFG_fig0;

    %RF_pd_distancesFG_bg0 = RF_hist_distancesFG_bg0./total_distancesFG_bg0;
    %RF_pd_DsFG_bg0 = RF_hist_DsFG_bg0./total_DsFG_bg0;
    RF_pd_disparitiesFG_bg0 = RF_hist_disparitiesFG_bg0./total_disparitiesFG_bg0;

    % averages - raw histograms

    % FG RFs, all points
    %distanceFG   = compute_mean_and_quantiles(RF_hist_distancesFG0); % normed distance m - all RFs
    %DFG          = compute_mean_and_quantiles(RF_hist_DsFG0); % normed distance m - all RFs
    dispFG       = compute_mean_and_quantiles(RF_hist_disparitiesFG0); % normed distance m - all RFs

    % FG RFs, figure regions
    %distance_fig   = compute_mean_and_quantiles(RF_hist_distancesFG_fig0); % normed distance m - all RFs
    %D_fig          = compute_mean_and_quantiles(RF_hist_DsFG_fig0); % normed distance m - all RFs
    disp_fig       = compute_mean_and_quantiles(RF_hist_disparitiesFG_fig0); % normed distance m - all RFs

    % FG RFs, ground regions
    %distance_bg   = compute_mean_and_quantiles(RF_hist_distancesFG_bg0); % normed distance m - all RFs
    %D_bg          = compute_mean_and_quantiles(RF_hist_DsFG_bg0); % normed distance m - all RFs
    disp_bg       = compute_mean_and_quantiles(RF_hist_disparitiesFG_bg0); % normed distance m - all RFs

    % averages - probability density for ratios

    % FG RFs, figure regions
    disp_fig_pd       = compute_mean_and_quantiles(RF_pd_disparitiesFG_fig0); % normed distance m - all RFs

    % FG RFs, ground regions
    disp_bg_pd       = compute_mean_and_quantiles(RF_pd_disparitiesFG_bg0); % normed distance m - all RFs


    % figure/ground probability plot
    if strcmp(cond{pl},'all')

        this_line_width = 2;

        figure; hold on; setupfig(4,4,12)
        a = area(cntr_disparity,dispFG.m);
        a.FaceColor = [0.75 0.75 0.75];
        a.EdgeColor = [0.75 0.75 0.75];
        plot(cntr_disparity,disp_fig.m,'-', 'Color',figColour,'LineWidth', this_line_width);
        plot(cntr_disparity, [disp_fig.m - disp_fig.CI; disp_fig.m + disp_fig.CI], ':', 'Color',figColour, 'LineWidth', this_line_width);
        plot(cntr_disparity,disp_bg.m,'-', 'Color',bgColour,'LineWidth', this_line_width);
        plot(cntr_disparity, [disp_bg.m - disp_bg.CI; disp_bg.m + disp_bg.CI], ':', 'Color',bgColour, 'LineWidth', this_line_width);
        xlabel('disparity (deg)'); axis square; box on;
        xlim([-1.5 1.5])
        ylabel('frequency');
        saveas(gcf,fullfile(save_dir,'RF_fg_disparity_probs.png'));

    end

    % figure probability ratio

    ratio_disp(pl,:)  = disp_fig_pd.m./disp_bg_pd.m;

    % bootstrap the 95% CI for these ratios

    doboot_disp = @(f,g) (mean(f, 'omitnan'))./(mean(g, 'omitnan')+eps); % we add eps to the denominator because the disparity histograms often contain zeros
    ratio_disp_ci(:,:,pl) = bootci(1000,doboot_disp,RF_pd_disparitiesFG_fig0,RF_pd_disparitiesFG_bg0);

end

%% plot histogram of distance difference between figure and ground

all_ind = find(strcmp(cond,'all'));

figure; hold on; setupfig(4,4,12)
h = histogram(FG_diff_disparity{all_ind},linspace(-1.5,1.5,21),'facecolor',[0.75 0.75 0.75],'facealpha',1);
xlabel('relative disparity between figure and ground (deg)'); axis square; box on;
ylabel('frequency')
h.FaceAlpha = 1;
xlim([-1.5 1.5])
saveas(gcf,fullfile(save_dir,'RF_fg_disparity_differences.png'));

% each eccentricity separately
ecc_inds = find(~strcmp(cond,'all'));
ecc_dark = [0.9 0.7 0.5 0.3]; % color mapping

figure; hold on; setupfig(4,4,12);

for ex = 1:numel(ecc_inds)
    % Extract data for the current eccentricity
    this_ecc = ecc_inds(ex);
    data = FG_diff_disparity{this_ecc};

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
ylabel('relative disparity between figure and ground (deg)');
axis square; box on;
xlim([0.5 4.5]); ylim([-0.1 0.6]);
plot([0.5 4.5],[0 0],'k-')
xticks(1:numel(ecc_inds));
xticklabels({'2.5','5','10','15'});

saveas(gcf,fullfile(save_dir,'RF_fg_diparity_differences_by_ecc.png'));

% store and save in table

T_DIFF_DEP = table(cond',phat_prop_fig_closer',pci_prop_fig_closer(:,1),pci_prop_fig_closer(:,2), ...
    mean_diff_Ds',mean_diff_disparity',stats_FG_disparity_diff.tstat',stats_FG_disparity_diff.df',eff_size_FG_disparity_diff',p_FG_disparity_diff',...
    'VariableNames',{'ecc','prop_fig_closer','prop_fig_closer_LCI','prop_fig_closer_HCI',...
    'mean_Diff_D','mean_Diff_disparity','tstat','df','cohensD','pval'});

display(T_DIFF_DEP);

save(fullfile(save_dir,['RFs_stats_disparity_diffs.mat']),'T_DIFF_DEP');


% run a depth difference anova on Eccentricity
% data appear roughly normal with good homogeneity of variances

% Initialize containers for combined data and group labels
allData = [];
group = [];

% Loop through each group to extract data and create group labels
for ex = 1:numel(ecc_inds)
    this_ecc    = ecc_inds(ex);
    groupData   = FG_diff_disparity{this_ecc}';          % Data for the current group
    allData     = [allData; groupData];  % Append the group's data to allData
    group       = [group; repmat(cond(this_ecc), numel(groupData), 1)];  % Create group labels
end

% Run one-way ANOVA
[p_disparity, tbl_disparity, stats_disparity] = anova1(allData, group);

% Display the ANOVA table
display('Ecc ANOVA for disparity diff')
disp(tbl_disparity);

SS_between = cell2mat(tbl_disparity(2, 2));  % Sum of Squares between groups (SSB)
SS_total = cell2mat(tbl_disparity(4, 2));    % Total Sum of Squares (SST)
eta_squared = SS_between / SS_total;
fprintf('Effect size (η²): %.4f\n', eta_squared);

% (Optional) Perform post-hoc multiple comparisons
[comparison_disparity,means_disparity,h_disparity,gnames_disparity] = multcompare(stats_disparity, "CriticalValueType","hsd");

resultsTable = create_eccentricity_comparison_table(stats_disparity, comparison_disparity, 'Mean Diff (deg)',allData,group);

% Display the Table
disp(resultsTable);

% Save the Table to CSV
writetable(resultsTable, fullfile(save_dir,['RFs_stats_disparity_diffs_pairwise.csv']));

% figure/ground probability ratio

figure; hold on; setupfig(4,4,12)
plot(cntr_disparity,ratio_disp(all_ind,:),'k-', 'LineWidth', this_line_width)
plot(cntr_disparity,ratio_disp_ci(1,:,all_ind),'k:', 'LineWidth', this_line_width)
plot(cntr_disparity,ratio_disp_ci(2,:,all_ind),'k:', 'LineWidth', this_line_width)
plot([cntr_disparity(1) cntr_disparity(end)],[1 1],'k-')
set(gca, 'YScale', 'log');
ax = gca;
ax.Color = 'none';
ylim([10^-4 10^4])
xlim([-1.5 1.5])
xlabel('disparity (deg)');
ylabel('figure/ground probability'); axis square; box on;
saveas(gcf,fullfile(save_dir,'RF_disparity_probability_ratios.png'));