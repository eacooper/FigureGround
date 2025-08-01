close all;

% assert a seed for randomization so that RFs are reproducible
rng(1999);

% clean up speed and orientation data
all_set_mag                                 = all_set_mag.*(FPS*actual_FOV_deg/horiz); %converts from pix/frame to visdeg/s
all_set_mag(all_set_mag <= magthreshold)    = magthreshold;
all_set_ori(all_set_mag == magthreshold)    = NaN;
all_set_ori                                 = wrapTo2Pi(all_set_ori);

%% set up for simulations

% Set save path for results
save_dir = fullfile(resultsPath, 'motion_RF_analysis');
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

% Number of frames in the video dataset
num_scenes = size(all_set_mag, 3);

% get indices for rows and columns of each frame
% col = x, row = y
[col_vals,row_vals] = meshgrid(1:dims(1),1:dims(2));

% scope of simulations to run
% number of times to sample a random pixel (this will de divided by 4 for 4 RF sizes)
num_simulations     = 800;

% set up histogram bins for motion analysis (note depth will be in normalized unit from 0-1
edges_mag   = linspace(0,1,51);
edges_ori   = linspace(0,2*pi,51);

% histogram bin centers
cntr_mag    = edges_mag(1:end-1) + (edges_mag(2) - edges_mag(1))/2;
cntr_ori    = edges_ori(1:end-1) + (edges_ori(2) - edges_ori(1))/2;
cntr_ori    = circshift(cntr_ori, 1);

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

        % Get all the info pertaining to this specific frame
        scene_mag        = all_set_mag(:, :, rand_frame_ind);
        scene_ori        = all_set_ori(:, :, rand_frame_ind);
        scene_label      = all_set_labels(:, :, rand_frame_ind);
        scene_saliency   = all_set_most_salient(:, :, rand_frame_ind);

        %% random RF center

        % random row/column for center of RF, outside of buffer region
        r_rand = randi([buffer dims(2)-buffer]);
        c_rand = randi([buffer dims(1)-buffer]);

        % binary mask of this RF
        RFmask = create_RF_mask(r_rand,c_rand,RF_r_pix,col_vals,row_vals);

        % indices of pixels in RF
        pixels_in_RF = find(RFmask);

        % restart if there's too much sky/nans (50% of pixels)
        if sum(isnan(scene_ori(pixels_in_RF)))/numel(scene_ori(pixels_in_RF)) >= 0.5
            invalid_count = invalid_count + 1;
            %display('too many NaNs, moving on')
            continue
        end

        % check that at least 25% are figure and 25% are ground
        % figure out which pixels in the RF are in the figure, which in the ground
        figs                = all_set_figs(:,:,rand_frame_ind);
        bgs                 = all_set_bgs(:,:,rand_frame_ind);

        pixels_in_fig   = find(figs & RFmask);
        pixels_in_bg    = find(bgs & RFmask);

        % compute and display the figure/ground proportions
        perc_fig    = 100*numel(pixels_in_fig)/numel(pixels_in_RF);
        perc_bg     = 100*numel(pixels_in_bg)/numel(pixels_in_RF);

        % restart if there's not at least 25% of each
        if perc_fig < 25 || perc_bg < 25
            invalid_count = invalid_count + 1;
            %display('not enough FG pixels, moving on')
            continue
        end

        % restart if any region has all the same values, because then we
        % cant compute a normalized histogram (this is really rare)
        if range(scene_mag(pixels_in_RF)) < 0.001 || range(scene_mag(pixels_in_fig)) < eps || range(scene_mag(pixels_in_bg)) < eps ||...
                range(scene_ori(pixels_in_RF)) < 0.001 || range(scene_ori(pixels_in_fig)) < eps || range(scene_ori(pixels_in_bg)) < eps
            invalid_count = invalid_count + 1;
            display('too low variance')
            continue
        end

        display(['percent fig: ' num2str(perc_fig) ' | percent bg: ' num2str(perc_bg)])

        % simulate a random fixation point for subtracting out smooth tracking motion
        fix = generate_random_fixation_point_motion(r_rand,c_rand,row_vals,col_vals,RF_r_pix,scene_mag,scene_ori,scene_saliency,magthreshold);

        % restart if there's not a valid fixation point
        if isnan(fix.mag)
            invalid_count = invalid_count + 1;
            %display('not a valid fixation, moving on')
            continue
        end

        % if all is good, keep going
        ii_r = 1;
        display('...')

    end

    % now that we have a good RF...

    % store epoque info
    scene_ep_num        = all_set_epoque(rand_frame_ind);

    % get indices of figure and ground pixels within pixels_in_RF vector, so we can subset the analysis
    fig_indices = find(ismember(pixels_in_RF, pixels_in_fig));
    bg_indices = find(ismember(pixels_in_RF, pixels_in_bg));

    % grab vector of all point speeds and orientations in RF
    RF_mag = scene_mag(pixels_in_RF);
    RF_ori = scene_ori(pixels_in_RF);

    % perform correction to model smooth pursuit of fixation point

    % convert scene motion to cartesian
    scene_x = scene_mag.*cos(scene_ori);
    scene_y = scene_mag.*sin(scene_ori);

    % convert pursuit motion to cartesian
    fix.x = fix.mag.*cos(fix.ori);
    fix.y = fix.mag.*sin(fix.ori);

    % correct for pursuit
    scene_x_pursuit = scene_x - fix.x;
    scene_y_pursuit = scene_y - fix.y;

    % convert back to mag/ori
    [scene_mag_pursuit, scene_ori_pursuit] = calc_mag_oris(scene_x_pursuit, scene_y_pursuit);
    scene_ori_pursuit = wrapTo2Pi(scene_ori_pursuit);

    % determine whether the fixation point is on a figure-region
    fix.fig = figs(fix.r,fix.c);

    % grab values in RF
    RF_mag_pursuit = scene_mag_pursuit(pixels_in_RF);
    RF_ori_pursuit = scene_ori_pursuit(pixels_in_RF);

    % convert speed to log speed
    RF_mag_log = log10(RF_mag);
    RF_mag_log_pursuit = log10(RF_mag_pursuit);

    % normalize speed data
    RF_mag_norm                  = normalize_data(RF_mag,lowerPerc,upperPerc);
    RF_mag_norm_log              = normalize_data(RF_mag_log,lowerPerc,upperPerc);
    RF_ori(isnan(RF_mag_norm))   = NaN;

    RF_mag_norm_pursuit                  = normalize_data(RF_mag_pursuit,lowerPerc,upperPerc);
    RF_mag_norm_log_pursuit              = normalize_data(RF_mag_log_pursuit,lowerPerc,upperPerc);
    RF_ori_pursuit(isnan(RF_mag_norm_pursuit))   = NaN;

    % compute all histograms
    hist_rf_mag           = histcounts(RF_mag_norm,edges_mag,'Normalization','count');
    hist_rf_mag_log       = histcounts(RF_mag_norm_log,edges_mag,'Normalization','count');
    hist_rf_ori_unnorm    = histcounts(RF_ori,edges_ori,'Normalization','count');

    hist_rf_mag_pursuit           = histcounts(RF_mag_norm_pursuit,edges_mag,'Normalization','count');
    hist_rf_mag_log_pursuit       = histcounts(RF_mag_norm_log_pursuit,edges_mag,'Normalization','count');
    hist_rf_ori_unnorm_pursuit    = histcounts(RF_ori_pursuit,edges_ori,'Normalization','count');

    % normalize direction data, return number of bins to displace for normalization
    [hist_rf_ori, max_ind, ~] = normalize_ori_data(hist_rf_ori_unnorm, cntr_ori);
    [hist_rf_ori_pursuit, max_ind_pursuit, ~] = normalize_ori_data(hist_rf_ori_unnorm_pursuit, cntr_ori);

    % figure histos
    if ~isempty(fig_indices)
        hist_rf_mag_fig           = histcounts(RF_mag_norm(fig_indices),edges_mag,'Normalization','count');
        hist_rf_mag_log_fig       = histcounts(RF_mag_norm_log(fig_indices),edges_mag,'Normalization','count');
        hist_rf_ori_fig_unnorm    = histcounts(RF_ori(fig_indices), edges_ori, 'Normalization','count');
        hist_rf_ori_fig           = normalize_ori_data(hist_rf_ori_fig_unnorm, cntr_ori, max_ind);

        hist_rf_mag_fig_pursuit           = histcounts(RF_mag_norm_pursuit(fig_indices),edges_mag,'Normalization','count');
        hist_rf_mag_log_fig_pursuit       = histcounts(RF_mag_norm_log_pursuit(fig_indices),edges_mag,'Normalization','count');
        hist_rf_ori_fig_unnorm_pursuit    = histcounts(RF_ori_pursuit(fig_indices), edges_ori, 'Normalization','count');
        hist_rf_ori_fig_pursuit           = normalize_ori_data(hist_rf_ori_fig_unnorm_pursuit, cntr_ori, max_ind_pursuit);
    else
        hist_rf_mag_fig     = NaN(1,numel(edges_mag)-1);
        hist_rf_mag_log_fig = NaN(1,numel(edges_mag)-1);
        hist_rf_ori_fig     = NaN(1,numel(edges_ori)-1);
        hist_rf_mag_fig_pursuit     = NaN(1,numel(edges_mag)-1);
        hist_rf_mag_log_fig_pursuit = NaN(1,numel(edges_mag)-1);
        hist_rf_ori_fig_pursuit     = NaN(1,numel(edges_ori)-1);
    end

    % ground histos
    if ~isempty(bg_indices)
        hist_rf_mag_bg          = histcounts(RF_mag_norm(bg_indices),edges_mag,'Normalization','count');
        hist_rf_mag_log_bg      = histcounts(RF_mag_norm_log(bg_indices),edges_mag,'Normalization','count');
        hist_rf_ori_bg_unnorm    = histcounts(RF_ori(bg_indices), edges_ori, 'Normalization','count');
        hist_rf_ori_bg           = normalize_ori_data(hist_rf_ori_bg_unnorm, cntr_ori, max_ind);

        hist_rf_mag_bg_pursuit          = histcounts(RF_mag_norm_pursuit(bg_indices),edges_mag,'Normalization','count');
        hist_rf_mag_log_bg_pursuit      = histcounts(RF_mag_norm_log_pursuit(bg_indices),edges_mag,'Normalization','count');
        hist_rf_ori_bg_unnorm_pursuit    = histcounts(RF_ori_pursuit(bg_indices), edges_ori, 'Normalization','count');
        hist_rf_ori_bg_pursuit           = normalize_ori_data(hist_rf_ori_bg_unnorm_pursuit, cntr_ori, max_ind_pursuit);
    else
        hist_rf_mag_bg    = NaN(1,numel(edges_mag)-1);
        hist_rf_mag_log_bg    = NaN(1,numel(edges_mag)-1);
        hist_rf_ori_bg    = NaN(1,numel(edges_ori)-1);
        hist_rf_mag_bg_pursuit    = NaN(1,numel(edges_mag)-1);
        hist_rf_mag_log_bg_pursuit    = NaN(1,numel(edges_mag)-1);
        hist_rf_ori_bg_pursuit    = NaN(1,numel(edges_ori)-1);
    end

    %% aggregate RF info
    RF_scenes(sim_count)             = rand_frame_ind;
    RF_epoque(sim_count)             = scene_ep_num;
    RF_diameters(sim_count)          = RF_r;
    RF_center_pix(sim_count,:)       = [r_rand c_rand];
    RF_perc_fig(sim_count)           = perc_fig;
    RF_perc_bg(sim_count)            = perc_bg;

    % fixation info
    RF_sim_fixation_pix(sim_count,:) = [fix.r fix.c];
    RF_sim_fixation_mag(sim_count) = fix.mag;
    RF_sim_fixation_ori(sim_count) = fix.ori;
    RF_sim_fixation_fig(sim_count,:) = [fix.fig];

    % histograms
    RF_hist_mag(sim_count,:)         = hist_rf_mag;
    RF_hist_mag_log(sim_count,:)     = hist_rf_mag_log;
    RF_hist_ori(sim_count,:)         = hist_rf_ori;

    RF_hist_mag_fig(sim_count,:)     = hist_rf_mag_fig;
    RF_hist_mag_log_fig(sim_count,:) = hist_rf_mag_log_fig;
    RF_hist_ori_fig(sim_count,:)     = hist_rf_ori_fig;

    RF_hist_mag_bg(sim_count,:)      = hist_rf_mag_bg;
    RF_hist_mag_log_bg(sim_count,:)  = hist_rf_mag_log_bg;
    RF_hist_ori_bg(sim_count,:)      = hist_rf_ori_bg;

    % histograms - pursuit
    RF_hist_mag_pursuit(sim_count,:)         = hist_rf_mag_pursuit;
    RF_hist_mag_log_pursuit(sim_count,:)     = hist_rf_mag_log_pursuit;
    RF_hist_ori_pursuit(sim_count,:)         = hist_rf_ori_pursuit;

    RF_hist_mag_fig_pursuit(sim_count,:)     = hist_rf_mag_fig_pursuit;
    RF_hist_mag_log_fig_pursuit(sim_count,:) = hist_rf_mag_log_fig_pursuit;
    RF_hist_ori_fig_pursuit(sim_count,:)     = hist_rf_ori_fig_pursuit;

    RF_hist_mag_bg_pursuit(sim_count,:)      = hist_rf_mag_bg_pursuit;
    RF_hist_mag_log_bg_pursuit(sim_count,:)  = hist_rf_mag_log_bg_pursuit;
    RF_hist_ori_bg_pursuit(sim_count,:)      = hist_rf_ori_bg_pursuit;

    % central tendencies and variance for speed

    RF_mag_med_all(sim_count)        = median(scene_mag(pixels_in_RF),'omitnan');
    RF_mag_mean_all(sim_count)       = mean(scene_mag(pixels_in_RF),'omitnan');

    RF_mag_med_all_pursuit(sim_count)  = median(scene_mag_pursuit(pixels_in_RF),'omitnan');
    RF_mag_mean_all_pursuit(sim_count) = mean(scene_mag_pursuit(pixels_in_RF),'omitnan');

    % central tendencies and variance for orientation
    % for orientation, remove nans manually because circ stats won't do it
    RF_ori_mean_all(sim_count) = circ_mean_omit_nan(scene_ori(pixels_in_RF));
    RF_ori_var_all(sim_count)  = circ_var_omit_nan(scene_ori(pixels_in_RF));

    RF_ori_mean_all_pursuit(sim_count) = circ_mean_omit_nan(scene_ori_pursuit(pixels_in_RF));
    RF_ori_var_all_pursuit(sim_count)  = circ_var_omit_nan(scene_ori_pursuit(pixels_in_RF));

    if ~isempty(pixels_in_fig)

        RF_mag_med_fig(sim_count)        = median(scene_mag(pixels_in_fig),'omitnan');
        RF_mag_mean_fig(sim_count)       = mean(scene_mag(pixels_in_fig),'omitnan');

        RF_mag_med_fig_pursuit(sim_count)        = median(scene_mag_pursuit(pixels_in_fig),'omitnan');
        RF_mag_mean_fig_pursuit(sim_count)       = mean(scene_mag_pursuit(pixels_in_fig),'omitnan');

        RF_ori_mean_fig(sim_count) = circ_mean_omit_nan(scene_ori(pixels_in_fig));
        RF_ori_var_fig(sim_count)  = circ_var_omit_nan(scene_ori(pixels_in_fig));

        RF_ori_mean_fig_pursuit(sim_count) = circ_mean_omit_nan(scene_ori_pursuit(pixels_in_fig));
        RF_ori_var_fig_pursuit(sim_count)  = circ_var_omit_nan(scene_ori_pursuit(pixels_in_fig));

    else
        RF_mag_med_fig(sim_count) = NaN; RF_mag_mean_fig(sim_count) = NaN;
        RF_ori_mean_fig(sim_count) = NaN; RF_ori_var_fig(sim_count) = NaN;
        RF_mag_med_fig_pursuit(sim_count) = NaN; RF_mag_mean_fig_pursuit(sim_count) = NaN;
        RF_ori_mean_fig_pursuit(sim_count) = NaN; RF_ori_var_fig_pursuit(sim_count) = NaN;
    end

    if ~isempty(pixels_in_bg)

        RF_mag_med_bg(sim_count)        = median(scene_mag(pixels_in_bg),'omitnan');
        RF_mag_mean_bg(sim_count)        = mean(scene_mag(pixels_in_bg),'omitnan');

        RF_mag_med_bg_pursuit(sim_count)        = median(scene_mag_pursuit(pixels_in_bg),'omitnan');
        RF_mag_mean_bg_pursuit(sim_count)        = mean(scene_mag_pursuit(pixels_in_bg),'omitnan');

        RF_ori_mean_bg(sim_count) = circ_mean_omit_nan(scene_ori(pixels_in_bg));
        RF_ori_var_bg(sim_count)  = circ_var_omit_nan(scene_ori(pixels_in_bg));

        RF_ori_mean_bg_pursuit(sim_count) = circ_mean_omit_nan(scene_ori_pursuit(pixels_in_bg));
        RF_ori_var_bg_pursuit(sim_count)  = circ_var_omit_nan(scene_ori_pursuit(pixels_in_bg));

    else
        RF_mag_med_bg(sim_count) = NaN; RF_mag_mean_bg(sim_count) = NaN;
        RF_ori_mean_bg(sim_count) = NaN; RF_ori_var_bg(sim_count) = NaN;
        RF_mag_med_bg_pursuit(sim_count) = NaN; RF_mag_mean_bg_pursuit(sim_count) = NaN;
        RF_ori_mean_bg_pursuit(sim_count) = NaN; RF_ori_var_bg_pursuit(sim_count) = NaN;
    end

    sim_count = sim_count + 1;

    close all;

end

% save files for next analysis step
save(fullfile(save_dir,'RF_motion_measurements.mat'),...
    'cntr_mag','cntr_ori', 'edges_ori', ...
    'RF_scenes','RF_diameters','RF_center_pix','RF_perc_fig','RF_perc_bg',...
    'RF_hist_mag','RF_hist_mag_log','RF_hist_ori',...
    'RF_hist_mag_fig','RF_hist_mag_log_fig','RF_hist_ori_fig',...
    'RF_hist_mag_bg','RF_hist_mag_log_bg','RF_hist_ori_bg', ...
    'RF_mag_med_fig', 'RF_mag_mean_fig','RF_ori_mean_fig', 'RF_ori_var_fig',...
    'RF_mag_med_bg', 'RF_mag_mean_bg', 'RF_ori_mean_bg', 'RF_ori_var_bg',...
    'RF_mag_med_all', 'RF_mag_mean_all', 'RF_ori_mean_all', 'RF_ori_var_all',...
    'RF_sim_fixation_pix','RF_sim_fixation_mag','RF_sim_fixation_ori','RF_sim_fixation_fig',...
    'RF_hist_mag_pursuit','RF_hist_mag_log_pursuit','RF_hist_ori_pursuit',...
    'RF_hist_mag_fig_pursuit','RF_hist_mag_log_fig_pursuit','RF_hist_ori_fig_pursuit',...
    'RF_hist_mag_bg_pursuit','RF_hist_mag_log_bg_pursuit','RF_hist_ori_bg_pursuit', ...
    'RF_mag_med_fig_pursuit', 'RF_mag_mean_fig_pursuit','RF_ori_mean_fig_pursuit', 'RF_ori_var_fig_pursuit',...
    'RF_mag_med_bg_pursuit', 'RF_mag_mean_bg_pursuit', 'RF_ori_mean_bg_pursuit', 'RF_ori_var_bg_pursuit',...
    'RF_mag_med_all_pursuit', 'RF_mag_mean_all_pursuit', 'RF_ori_mean_all_pursuit', 'RF_ori_var_all_pursuit');
