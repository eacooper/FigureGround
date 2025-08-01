function [FPName] =  get_full_videoframe_filename(vidset_ID, vid_num, im_num, im_folder, isFull)

    vid_name = strcat(vidset_ID, num2str(vid_num));

    filePattern = fullfile(im_folder, strcat(vid_name, '*-', num2str(im_num), '.raw'));

    images_list = dir(filePattern);
    
    %checks if the file doesn't exist
    if isempty(images_list)
        FPName = 0;
        return
    end

    FPName = fullfile(images_list(1).folder, images_list(1).name);

    if ~isFull
        FPName = images_list(1).name;
    end

end