% load in, preprocess, and save the annotations

% paths
saveDir = fullfile(resultsPath,'motion_perVideo_info');

% type of frames
frame_types = {"ep1a", ...
    "ep1b", ...
    "ep1c", ...
    "ep2a", ...
    "ep2b", ...
    "ep2c"};

% For each video set
for v = 1:length(vidset_IDs)

    display(['OF array compilation for ' vidset_IDs{v}]);
    vidset_ID   = vidset_IDs{v};

    % Video number counter
    vc = 0;

    % initialize
    all_set_labels = []; all_set_figs = []; all_set_bgs = []; all_set_border = [];

    for w = 1:vid_nums(v)
        
        vid_num  = w;
        vid_name = strcat(vidset_ID, num2str(vid_num));
        display(vid_name);

        % skip the bad ones
        if ismember(w, ALL_ex{v})
            display('bad video, skipping')
            continue;
        end

        vc = vc + 1;

        for u = 1:num_epoques*3

            annotations_dir = fullfile(dataPath, 'Annotations', vidset_ID);
        
            % load per-dataset structures
            filePattern = fullfile(annotations_dir,strcat('*', frame_types{u},'*', vidset_ID, num2str(w), '-*.png'));
            files_list = dir(filePattern);
            loadoName = files_list(1).name;
            this_mask = imread(fullfile(annotations_dir,loadoName));

            % convert RGB to categorical labels
            [labels, figs, bgs, dict] = convertrgb2labels(this_mask);

            % make a boundary mask
            bnd         = boundarymask(labels);           % boundary between labelled regions
            bnd_f       = fspecial('average',bnd_buff);   % filter to thicken boundary
            bnd_thick   = imfilter(double(bnd),bnd_f);    % apply filter
            bnd_thick(bnd_thick > 0) = 1;
        
            % apply boundary mask to fig/ground masks
            figs(bnd_thick > 0) = 0;
            bgs(bnd_thick > 0)  = 0;
        
            % combine for visual
            fgbg = figs;
            fgbg(bgs == 1) = 0.5;
        
            n = (vc-1)*num_epoques*3 + u;
            display(n)

            %% store for analysis
            all_set_labels(:,:,n)       = labels;
            all_set_figs(:,:,n)         = figs;
            all_set_bgs(:,:,n)          = bgs;
            all_set_border(:,:,n)       = bnd;

        end

    end

    % Save the annotations
    saveName = fullfile(saveDir, strcat(vidset_ID,'_anns_allTOG.mat'));
    save(saveName, '-v7.3', 'all_set_labels', 'all_set_figs', 'all_set_bgs', 'all_set_border');

    %% use annotations to remove sky

    % optic flow
    filePattern     = fullfile(saveDir, strcat(vidset_ID, '_MagOri_allTOG.mat'));
    load(filePattern);

    all_set_mag(all_set_labels == 1) = NaN;
    all_set_ori(all_set_labels == 1) = NaN;

    % Save
    saveName = fullfile(saveDir, strcat(vidset_ID,'_MagOri_allTOG_skygone.mat'));
    save(saveName, '-v7.3',   'all_set_mag', 'all_set_ori');

end

