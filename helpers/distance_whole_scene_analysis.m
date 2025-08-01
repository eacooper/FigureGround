% Absolute distance across the whole scenes
close all;

% Set save path for visuals
save_dir = fullfile(resultsPath, 'distance_whole_scene_analysis');
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


%% Global histograms

% set up histogram bins for linear and log distance (m)
edges_log    = logspace(log10(2),log10(100),51);

% histogram bin centers
cntr_log       = edges_log(1:end-1) + (edges_log(2) - edges_log(1))/2;

% histogram bin heights 
hAll_log    = histcounts(all_set_distance(:),edges_log,'normalization','count');
hFg_log     = histcounts(all_set_distance(logical(all_set_figs)),edges_log,'normalization','count');
hBg_log     = histcounts(all_set_distance(logical(all_set_bgs)),edges_log,'normalization','count');



% distance in meters (log)
figure; hold on; setupfig(4,3,12); %title('Distance')
a = area(cntr_log(2:end-1),hAll_log(2:end-1));
h(1) = plot(cntr_log(2:end-1),hFg_log(2:end-1),'-','color',figColour,'linewidth',2);
h(2) = plot(cntr_log(2:end-1),hBg_log(2:end-1),'-','color',bgColour,'linewidth',2);
set(gca,'xscale','log');
a.FaceColor = [0.75 0.75 0.75];
a.EdgeColor = [0.75 0.75 0.75];
xlim([2 85]); box on;
xlabel('distance (m)'); ylabel('frequency');
set(gca,'ytick',[]);
set(gca,'xtick',[2.5 5 10 20 40 80]);
legend(h,'figure','ground');

exportgraphics(gcf, fullfile(save_dir, strcat('whole_scene_distance.png')));

% stats comparing figure/ground median distance

% grab values to test
figDs = all_set_distance(logical(all_set_figs));
bgDs  = all_set_distance(logical(all_set_bgs));

% compute medians
median_all = median(all_set_distance(~isnan(all_set_distance)));
median_fig = median(figDs);
median_bg = median(bgDs);

% rank sum test
[p,h,stats] = ranksum(figDs,bgDs);

% effect size
% https://rpkgs.datanovia.com/rstatix/reference/wilcox_effsize.html
eff_size = stats.zval/sqrt(numel(figDs) + numel(bgDs));

% store and save in table
T_whole_scene_DEP = table(median_all,median_fig,median_bg,stats.zval,stats.ranksum,eff_size,p,'VariableNames',{'median_all','median_fig','median_bg','zval','ranksum','eff_size','pval'});
display(T_whole_scene_DEP);

save(fullfile(save_dir, strcat('whole_scene_stats.mat')),'T_whole_scene_DEP');