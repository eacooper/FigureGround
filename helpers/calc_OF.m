function[Vx, Vy, vid_segment] = calc_OF(vidset_ID, vid_num, frame_num, taps, scale, nlevels, im_folder, yesUndistort)
%
% calculate optic flow using matlab deformable image registration

%% set up
if yesUndistort
    load('cameraParams_checkerchecker.mat');
end

%Calculating OF imreg
vid_segment     = [];
all_Vx          = [];
all_Vy          = [];

% deal with even and odd taps
if mod(taps,2)
    inds = -floor(taps/2):taps-ceil(taps/2);
else
    inds = -floor(taps/2)+1:taps-ceil(taps/2);
end

%% Get and prep the frames in the tap
for k = inds
    
    % Get name of file
    fullFileName = get_full_videoframe_filename(vidset_ID, vid_num, frame_num+k, im_folder, true);
    
    % load file
    this_frame = load_raw_file(fullFileName, [1024, 1024])./2^16; % divide 2^16  to get grayscale
    
    if yesUndistort
        this_frame = undistortImage(this_frame,cameraParams_checkerchecker);
    end
    
    % apply median filter to reduce noise
    this_frame = medfilt2( this_frame, [3 3] );

    % resize if running coarse flow
    if scale < 1
        this_frame = imresize(this_frame, scale);
    end
    
    vid_segment = cat(3, vid_segment, this_frame);
    
end

%% Calc OFs between each pair of subsequent frames
for k = 1:taps-1

    [this_Vx, this_Vy] = frame_pair_imreg(vid_segment(:, :, k), vid_segment(:, :, k+1), nlevels);
    all_Vx = cat(3, all_Vx, this_Vx);
    all_Vy = cat(3, all_Vy, this_Vy);

end

Vx = mean(all_Vx, 3);%.*(FPS*FOV/horiz);
Vy = mean(all_Vy, 3);%.*(FPS*FOV/horiz);


end

function [Vx, Vy] = frame_pair_imreg(i1, i2, levels)

    % histogram equalization
    x2 = histeq(i1);
    x1 = histeq(i2);
    
    % histogram match
    x2 = imhistmatch(x2,x1);
    x1 = imhistmatch(x1,x2);

    %GridSpacing must be a two-element vector, and its default value is [4 4]. 
    % Smaller values of GridSpacing specify a finer grid resolution.
    grid_spacing = [4 4];
    
    % Weighing factor for grid displacement regularization, specified as a nonnegative scalar.
    % A large value for GridRegularization can create a smooth output displacement field, whereas a small value can create more localized displacements.
    grid_regularization = 1;

    [DisplacementField, ~] = imregdeform(x2,x1,NumPyramidLevels=levels,GridSpacing=grid_spacing,GridRegularization=grid_regularization, DisplayProgress=false);

    % magnitude of motion in pixels/frame
    Vx = DisplacementField(:,:,1);
    Vy = DisplacementField(:,:,2);

end