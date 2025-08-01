close all;

% set paths
epoques_dir     = fullfile(codePath,'helpers','motion_epoque_inds');
raw_frames_dir  = fullfile(dataPath,'DataSets','Motion');

save_dir        = fullfile(resultsPath,'motion_perVideo_info');
if ~exist(save_dir,'dir')
    mkdir(save_dir);  % Create the directory if it does not exist
end

% motion types
motion_types = {'start','peak','stop'};

% For each video set
for v = 1:length(vidset_IDs)

    display(['OF calculation for ' vidset_IDs{v}]);

    % initialize empty matrices
    all_set_vidset = []; all_set_vidnum = []; all_set_epoque = [];
    all_set_motion_type = []; all_set_Vx = []; all_set_Vy = []; 
    all_set_ims = []; all_set_mag = []; all_set_ori = [];

    % counter for all frames
    fc = 1;

    % for each video files in this video set
    for w = 1:vid_nums(v)

        % select this frameset
        vidset_ID   = vidset_IDs{v};
        vid_num     = w;
        vid_name    = strcat(vidset_ID, num2str(vid_num));
        display(vid_name);

        % skip the bad ones
        if ismember(w, ALL_ex{v})
            display('bad video, skipping')
            continue;
        end

        % location of directory containing the selected epoque frames of this video
        frames_dir = fullfile(raw_frames_dir, vidset_ID);

        %Loads in the file with pre-selected epoque frame indices
        load(fullfile(epoques_dir, strcat(vid_name,'_max', num2str(num_epoques), '.mat')));

        % for each epoque
        for epc = 1:length(maxima_inds)

            display(['epoque ' num2str(epc)]);

            % for each motion type
            for m = 1:length(motion_types)

                typeName = motion_types{m};
                display(['analyzing motion at ' typeName]);

                % get the frame
                switch(typeName)
                    case 'start'
                        frame_num = start_inds(epc);  % grab the frame number
                    case 'peak'
                        frame_num = maxima_inds(epc);
                    case 'stop'
                        frame_num = stop_inds(epc);
                    otherwise
                end

                % calculate the optic flow
                the_works_OF;                 

                % store the info
                all_set_vidset{fc}              = vidset_ID;
                all_set_vidnum(fc)              = vid_num;
                all_set_epoque(fc)              = epc;
                all_set_motion_type{fc}         = typeName;

                all_set_Vx(:, :, fc)       = Vx;
                all_set_Vy(:, :, fc)       = Vy;
                all_set_ims(:, :, :, fc)   = vid_segment; % images will be a 1024:1024:taps:total_num_frames array now I think
                all_set_mag(:, :, fc)      = this_mag;
                all_set_ori(:, :, fc)      = this_ori;

                fc = fc+1;                     % increment counter for the all together matrix

            end

        end

    end

    %% Save all-dataset parameters
   
    % Save the metadata
    saveName = fullfile(save_dir, strcat(vidset_ID,'_MetaData_allTOG.mat'));
    save(saveName, '-v7.3',   'all_set_vidset', 'all_set_vidnum','all_set_epoque','all_set_motion_type');

    % Save the speed/direction data
    saveName = fullfile(save_dir, strcat(vidset_ID,'_MagOri_allTOG.mat'));
    save(saveName, '-v7.3',   'all_set_mag', 'all_set_ori');

    % Save the images
    saveName = fullfile(save_dir, strcat(vidset_ID,'_ims_allTOG.mat'));
    save(saveName, '-v7.3', 'all_set_ims');

    % Save the OF vectors
    saveName = fullfile(save_dir, strcat(vidset_ID,'_Vec_allTOG.mat'));
    save(saveName, '-v7.3', 'all_set_Vx', 'all_set_Vy');

end
