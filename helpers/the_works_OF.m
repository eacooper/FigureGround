% run full resolution optic flow

% Get the full path filename
fullFrameName = get_full_videoframe_filename(vidset_ID, vid_num, frame_num, frames_dir, true);

% Calculate the final optic flow
[Vx, Vy, vid_segment] = calc_OF(vidset_ID, vid_num, frame_num, taps, scale, pyr_levels, frames_dir, yesUndistort);

% Calculate optic flow speed and direction
[this_mag, this_ori] = calc_mag_oris(Vx, Vy); %pixels/frame

% Convert to singles
Vx          = single(Vx);
Vy          = single(Vy);
this_mag    = single(this_mag);
this_ori    = single(this_ori);
vid_segment = single(vid_segment./2^16); % center frame in video
