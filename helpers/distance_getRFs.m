close all;

% assert a seed for randomization so that RFs are reproducible
rng(1986);

%% set up for simulations

% Set save path for results
save_dir = fullfile(resultsPath, 'distance_RF_analysis');
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

% get indices for rows and columns of each frame
% col = x, row = y
[col_vals,row_vals] = meshgrid(1:dims(1),1:dims(2));

% scope of simulations to run
% number of times to sample a random pixel (this will de divided by 4 for 4 RF sizes)
num_simulations     = 800; 

% set up histogram bins for motion analysis (note depth will be in normalized unit from 0-1
edges_distance  = linspace(0,1,51);
edges_disparity = linspace(-1.5,1.5,51);

% histogram bin centers
cntr_distance    = edges_distance(1:end-1) + (edges_distance(2) - edges_distance(1))/2;
cntr_disparity   = edges_disparity(1:end-1) + (edges_disparity(2) - edges_disparity(1))/2;

% counter for number of simulations
sim_count = 1;

%% run simulations

% run total number of simulations
while sim_count <= num_simulations

    display(['running random RF simulation ' num2str(sim_count) '...']);

    % determine which eccentricity / receptive field diameter to simulate based on count
    if sim_count <= num_simulations/4;          Rind = 1;
    elseif sim_count <= 2*num_simulations/4;    Rind = 2;
    elseif sim_count <= 3*num_simulations/4;    Rind = 3;
    else;                                       Rind = 4;
    end

    % store diameter in visual degrees and pixels
    RF_r        = RF_sizes(Rind);
    RF_r_pix    = RF_sizes_pix(Rind);

    % compute buffer around image edge beyond with RF centers cannot go
    % this is just is 1/2 the RF diameter
    buffer = ceil(RF_r_pix/2) + 1;

    % boolean to indicate when a valid RF is located (we will discard RFs
    % that are mostly sky/invalid samples
    ii_r  = 0;

    % counter for invalid trials
    invalid_count = 0;

    % while a valid RF hasn't been found, keep searching
    while ii_r < 1

        %% random scene

        % pick a random scene
        rand_frame_ind = randi(num_scenes);

        % get the distance map and xyz coords for this scene
        scene_distance  = all_set_distance(:, :, rand_frame_ind);
        scene_x         = all_set_x(:, :, rand_frame_ind);
        scene_y         = all_set_y(:, :, rand_frame_ind);
        scene_z         = all_set_z(:, :, rand_frame_ind);

        % get saliency
        scene_saliency    = all_set_most_salient(:, :, rand_frame_ind);

        % convert distance in meters to diopters
        scene_D         = 1./scene_distance;

        %% random RF center

        % random row/column for center of RF, outside of buffer region
        r_rand = randi([buffer dims(2)-buffer]);
        c_rand = randi([buffer dims(1)-buffer]);

        % binary mask of this RF
        RFmask = create_RF_mask(r_rand,c_rand,RF_r_pix,col_vals,row_vals);

        % indices of pixels in RF
        pixels_in_RF = find(RFmask);

        % restart if there's too much sky/nans (50% of pixels)
        if sum(isnan(scene_distance(pixels_in_RF)))/numel(scene_distance(pixels_in_RF)) >= 0.5
            invalid_count = invalid_count + 1;
            %display('too many NaNs, moving on')
            continue
        end

        % check that at least 25% are figure and 25% are ground
        % figure out which pixels in the RF are in the figure, which in the ground
        figs    = all_set_figs(:,:,rand_frame_ind);
        bgs     = all_set_bgs(:,:,rand_frame_ind);

        pixels_in_fig   = find(figs & RFmask);
        pixels_in_bg    = find(bgs & RFmask);

        % compute and display the figure/ground proportions
        perc_fig    = 100*numel(pixels_in_fig)/numel(pixels_in_RF);
        perc_bg     = 100*numel(pixels_in_bg)/numel(pixels_in_RF);

        % restart if there's not at least 25% of each
        if perc_fig < 25 || perc_bg < 25
            invalid_count = invalid_count + 1;
            continue
        end

        % restart if any region has all the same values, because then we
        % cant compute a normalized histogram (this is really rare)
        if range(scene_distance(pixels_in_RF)) < 0.001 || range(scene_distance(pixels_in_fig)) < eps || range(scene_distance(pixels_in_bg)) < eps
            invalid_count = invalid_count + 1;
            display('too low variance')
            continue
        end
        
        display(['percent fig: ' num2str(perc_fig) ' | percent bg: ' num2str(perc_bg)])

        % simulate a random fixation point in meters for computing disparity in visual degrees
        fix = generate_random_fixation_point_distance(r_rand,c_rand,row_vals,col_vals,RF_r_pix,scene_x,scene_y,scene_z,scene_distance,scene_saliency);

        % restart if there's not a valid fixation point
        if isnan(fix.z)
            invalid_count = invalid_count + 1;
            continue
        end

        ii_r = 1;
        display('...')

    end

    % now that we have a good RF...

    % get indices of figure and ground pixels within pixels_in_RF vector,
    % so we can subset the analysis
    fig_indices = find(ismember(pixels_in_RF, pixels_in_fig));
    bg_indices = find(ismember(pixels_in_RF, pixels_in_bg));

    % grab vector of all point locations and distances in RF
    RF_distance = scene_distance(pixels_in_RF);
    RF_D        = scene_D(pixels_in_RF);
    RF_x        = scene_x(pixels_in_RF);
    RF_y        = scene_y(pixels_in_RF);
    RF_z        = scene_z(pixels_in_RF);

    % compute disparity of all points in the RF, based on this fixation point
    RF_disparity = calc_disparity([fix.x, fix.y, fix.z], [RF_x RF_y RF_z], eye_1, eye_2);

    % normalize distance data
    RF_distance_norm    = normalize_data(RF_distance,lowerPerc,upperPerc);
    RF_D_norm           = normalize_data(RF_D,lowerPerc,upperPerc);

    % compute all histograms
    hist_rf_distance    = histcounts(RF_distance_norm,edges_distance,'Normalization','count');
    hist_rf_D           = histcounts(RF_D_norm,edges_distance,'Normalization','count');
    hist_rf_disparity   = histcounts(RF_disparity,edges_disparity,'Normalization','count');

    % figure histos
    if ~isempty(fig_indices)
        hist_rf_distance_fig    = histcounts(RF_distance_norm(fig_indices),edges_distance,'Normalization','count');
        hist_rf_D_fig           = histcounts(RF_D_norm(fig_indices),edges_distance,'Normalization','count');
        hist_rf_disparity_fig   = histcounts(RF_disparity(fig_indices),edges_disparity,'Normalization','count');
    else
        hist_rf_distance_fig    = nan(1,numel(edges_distance)-1);
        hist_rf_D_fig           = nan(1,numel(edges_distance)-1);
        hist_rf_disparity_fig   = nan(1,numel(edges_distance)-1);
    end

    % ground histos
    if ~isempty(bg_indices)
        hist_rf_distance_bg    = histcounts(RF_distance_norm(bg_indices),edges_distance,'Normalization','count');
        hist_rf_D_bg           = histcounts(RF_D_norm(bg_indices),edges_distance,'Normalization','count');
        hist_rf_disparity_bg   = histcounts(RF_disparity(bg_indices),edges_disparity,'Normalization','count');
    else
        hist_rf_distance_bg    = nan(1,numel(edges_distance)-1);
        hist_rf_D_bg           = nan(1,numel(edges_distance)-1);
        hist_rf_disparity_bg   = nan(1,numel(edges_distance)-1);
    end

    % aggregate RF info
    RF_scenes(sim_count)             = rand_frame_ind;
    RF_diameters(sim_count)          = RF_r;
    RF_center_pix(sim_count,:)       = [r_rand c_rand];
    RF_sim_fixation_r_c(sim_count,:) = [fix.r fix.c];
    RF_sim_fixation(sim_count,:)     = [fix.x fix.y fix.z];
    RF_perc_fig(sim_count)           = perc_fig;
    RF_perc_bg(sim_count)            = perc_bg;

    RF_hist_distances(sim_count,:)   = hist_rf_distance;
    RF_hist_Ds(sim_count,:)          = hist_rf_D;
    RF_hist_disparities(sim_count,:) = hist_rf_disparity;

    RF_hist_distances_fig(sim_count,:)   = hist_rf_distance_fig;
    RF_hist_Ds_fig(sim_count,:)          = hist_rf_D_fig;
    RF_hist_disparities_fig(sim_count,:) = hist_rf_disparity_fig;
    
    RF_hist_distances_bg(sim_count,:)   = hist_rf_distance_bg;
    RF_hist_Ds_bg(sim_count,:)          = hist_rf_D_bg;
    RF_hist_disparities_bg(sim_count,:) = hist_rf_disparity_bg;

    % mean/median distances in whole RF and in figure/ground regions
    RF_dist_med_whole(sim_count)        = median(RF_distance,'omitnan');
    RF_dist_mean_whole(sim_count)        = mean(RF_distance);

    RF_disparity_med_whole(sim_count)        = median(RF_disparity,'omitnan');
    RF_disparity_mean_whole(sim_count)        = mean(RF_disparity);
    
    if isempty(pixels_in_fig)
        RF_dist_med_fig(sim_count)       = NaN;
        RF_dist_mean_fig(sim_count)      = NaN;
        RF_disparity_med_fig(sim_count)  = NaN;
        RF_disparity_mean_fig(sim_count) = NaN;
    else
        RF_dist_med_fig(sim_count)         = median(RF_distance(fig_indices),'omitnan');
        RF_dist_mean_fig(sim_count)        = mean(RF_distance(fig_indices),'omitnan');
        RF_disparity_med_fig(sim_count)    = median(RF_disparity(fig_indices),'omitnan');
        RF_disparity_mean_fig(sim_count)   = mean(RF_disparity(fig_indices),'omitnan');
    end

    if isempty(pixels_in_bg)
        RF_dist_med_bg(sim_count)      = NaN;
        RF_dist_mean_bg(sim_count)      = NaN;
        RF_disparity_med_bg(sim_count)      = NaN;
        RF_disparity_mean_bg(sim_count)      = NaN;
    else
        RF_dist_med_bg(sim_count)        = median(RF_distance(bg_indices),'omitnan');
        RF_dist_mean_bg(sim_count)        = mean(RF_distance(bg_indices),'omitnan');
        RF_disparity_med_bg(sim_count)        = median(RF_disparity(bg_indices),'omitnan');
        RF_disparity_mean_bg(sim_count)        = mean(RF_disparity(bg_indices),'omitnan');
    end

    sim_count = sim_count + 1;

end

% save files for next analysis step
save(fullfile(save_dir,'RF_distance_measurements.mat'),...
    'cntr_distance','cntr_disparity', ...
    'RF_scenes','RF_diameters','RF_center_pix','RF_sim_fixation','RF_sim_fixation_r_c','RF_perc_fig','RF_perc_bg',...
    'RF_dist_med_whole','RF_dist_med_fig','RF_dist_med_bg',...
    'RF_dist_mean_whole','RF_dist_mean_fig','RF_dist_mean_bg',...
    'RF_disparity_med_whole','RF_disparity_med_fig','RF_disparity_med_bg',...
    'RF_disparity_mean_whole','RF_disparity_mean_fig','RF_disparity_mean_bg',...
    'RF_hist_distances','RF_hist_Ds','RF_hist_disparities',...
    'RF_hist_distances_fig','RF_hist_Ds_fig','RF_hist_disparities_fig',...
    'RF_hist_distances_bg','RF_hist_Ds_bg','RF_hist_disparities_bg','fix');

