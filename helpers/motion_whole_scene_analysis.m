% Motion across the whole scenes
close all;

% Set save path for visuals
save_dir = fullfile(resultsPath, 'motion_whole_scene_analysis');
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

% load annotation labels
[~, ~, ~, dict] = convertrgb2labels([]);

%% Histogram and percentages of labels

% report percent sky/invalid
display(['percent sky/invalid labels = ' num2str(100*sum(all_set_labels(:)==1)/numel(all_set_labels))])

% report perc fig/bg of valid pixels
display(['percent fig labels = ' num2str(100*sum(all_set_figs(:))/(sum(all_set_bgs(:)) + sum(all_set_figs(:))))])
display(['percent bg labels = ' num2str(100*sum(all_set_bgs(:))/(sum(all_set_bgs(:)) + sum(all_set_figs(:))))])

figure; hold on;
histogram(all_set_labels(all_set_labels > 3));
set(gca,'xtick',3:15,'xticklabel',dict.name(3:15));
saveas(gcf,fullfile(save_dir,'figure_categories.png'));


%% set up histograms

% Histogram bin edges
edges_mag_lin   = linspace(0, 150, 50);
edges_mag_log   = logspace(log10(magthreshold),log10(150),50);
edges_ori       = linspace(0,2*pi,51);

% Histogram bin centres
cntr_mag_lin   = edges_mag_lin(1:end-1) + (edges_mag_lin(2) - edges_mag_lin(1))/2;
cntr_mag_log   = edges_mag_log(1:end-1) + (edges_mag_log(2) - edges_mag_log(1))/2;
cntr_ori       = edges_ori(1:end-1) + (edges_ori(2) - edges_ori(1))/2;

%% analysis

% clean up speed and orientation data
all_set_mag_clean  = all_set_mag.*(FPS*actual_FOV_deg/horiz); %converts from pix/frame to visdeg/s

% report proportion less than magthreshold
prop_thresh = sum(all_set_mag_clean(:) < magthreshold)/numel(all_set_mag_clean);
display(['percent at or below speed threshold = ' num2str(100*prop_thresh)])

all_set_mag_clean(all_set_mag_clean <= magthreshold)    = NaN;
all_set_ori(isnan(all_set_mag_clean))                   = NaN;
all_set_ori                                             = wrapTo2Pi(all_set_ori);

all_set_figs = logical(all_set_figs);
all_set_bgs  = logical(all_set_bgs);

%% Fill up histograms
hist_mag_log        = histcounts(all_set_mag_clean(:),edges_mag_log,'normalization','count');
hist_mag_log_fig    = histcounts(all_set_mag_clean(all_set_figs),edges_mag_log,'normalization','count');
hist_mag_log_bg     = histcounts(all_set_mag_clean(all_set_bgs),edges_mag_log,'normalization','count');

hist_ori            = histcounts(all_set_ori(:),edges_ori,'normalization','count');
hist_ori_fig        = histcounts(all_set_ori(all_set_figs),edges_ori,'normalization','count');
hist_ori_bg         = histcounts(all_set_ori(all_set_bgs),edges_ori,'normalization','count');


%% Plot them
grey_colour = [0.7 0.7 0.7];

%% speed

% logX linY
f2 = figure; hold on; setupfig(4,3,12);

a = area(cntr_mag_log,hist_mag_log, 'FaceColor', grey_colour, 'EdgeColor', grey_colour);
plot(cntr_mag_log,hist_mag_log_fig,'-','color',figColour,'linewidth',2)
plot(cntr_mag_log,hist_mag_log_bg,'-','color',bgColour,'linewidth',2)
set(gca,'ytick',[])
xlim([.5 130]); box on;
set(gca,'xscale','log');
set(gca,'xtick',[.5 1 3 10 30 100],'xticklabel',{'0.5','1.0','3.0','10','30','100'})
xlabel('speed (deg/s)');
ylabel('frequency');

exportgraphics(f2, fullfile(save_dir, strcat('whole_scene_logXlinY_hist_mag.png')));


%% orientation polar plot

p1 = figure;
polarplot(edges_ori,hist_ori([1:end 1]),'-','color',grey_colour,'linewidth',2);
hold on;
setupfig(4,3,12);
polarplot(edges_ori,hist_ori_fig([1:end 1]),'-','color', figColour,'linewidth',2);
polarplot(edges_ori,hist_ori_bg([1:end 1]),'-','color', bgColour,'linewidth',2);

thetaticks(0:30:330);
rticks([]);
exportgraphics(p1, fullfile(save_dir, strcat('whole_scene_hist_ori.png')));


%% stats

% speed distributions

% grab values to test
mag_fig = all_set_mag_clean(all_set_figs);
mag_bg  = all_set_mag_clean(all_set_bgs);

% compute medians
median_mag_all = median(all_set_mag_clean(~isnan(all_set_mag_clean)), 'all');
median_mag_fig = median(mag_fig, 'all', 'omitnan');
median_mag_bg = median(mag_bg, 'all', 'omitnan');

mean_mag_all = mean(all_set_mag_clean(~isnan(all_set_mag_clean)), 'all');
mean_mag_fig = mean(mag_fig, 'all', 'omitnan');
mean_mag_bg = mean(mag_bg, 'all', 'omitnan');

% rank sum test
[p_mag,h_mag,stats_mag] = ranksum(mag_fig,mag_bg);

% effect size
% https://rpkgs.datanovia.com/rstatix/reference/wilcox_effsize.html
eff_size_mag = stats_mag.zval/sqrt(numel(mag_fig) + numel(mag_bg));

% store and save in table
T_whole_scene_MAG = table(median_mag_all,median_mag_fig,median_mag_bg,mean_mag_all, mean_mag_fig,mean_mag_bg, stats_mag.zval,stats_mag.ranksum,eff_size_mag,p_mag,'VariableNames',{'median_mag_all','median_mag_fig','median_mag_bg','mean_mag_all','mean_mag_fig','mean_mag_bg','zval_mag','ranksum_mag','eff_size_mag','pval_mag'});
display(T_whole_scene_MAG);

% orientation distributions

% grab values to test
ori_fig = intersect(all_set_ori(all_set_figs), all_set_ori(~isnan(all_set_mag_clean)));
ori_bg = intersect(all_set_ori(all_set_bgs), all_set_ori(~isnan(all_set_mag_clean)));

% compute medians
mean_ori_all = circ_mean(all_set_ori(~isnan(all_set_mag_clean)));
mean_ori_fig = circ_mean(ori_fig(~isnan(ori_fig)));
mean_ori_bg = circ_mean(ori_bg(~isnan(ori_bg)));

% kuiper test
[p_ori,k_ori,K_ori] = circ_kuipertest_new(ori_fig,ori_bg);

% effect size
% https://rpkgs.datanovia.com/rstatix/reference/wilcox_effsize.html
%eff_size = stats.zval/sqrt(numel(ori_fig) + numel(ori_bg));
display("Find how to do effect size for kuiper test vals");

% store and save in table
T_whole_scene_ORI = table(mean_ori_all,mean_ori_fig,mean_ori_bg,p_ori, k_ori, 'VariableNames',{'mean_ori_all','mean_ori_fig','mean_ori_bg','pval_ori', 'test stat'});
display(T_whole_scene_ORI);

% orientations versus isotropic

% rayleigh test
[p_ori_fig,z_ori_fig] = circ_rtest(cntr_ori,hist_ori_fig,cntr_ori(3) - cntr_ori(2));
[p_ori_bg,z_ori_bg] = circ_rtest(cntr_ori,hist_ori_bg,cntr_ori(3) - cntr_ori(2));
T_diff_from_zero_ORI = table(p_ori_fig,z_ori_fig,p_ori_bg,z_ori_bg, 'VariableNames',{'p_ori_fig','z_ori_fig','p_ori_bg','z_ori_bg'});
display(T_diff_from_zero_ORI);

% store in table
save(fullfile(save_dir, strcat('whole_scene_stats.mat')),'T_whole_scene_MAG', 'T_whole_scene_ORI');


