function [im] = load_raw_file(FileName, im_dims)

%FileName = ex. 'Raw16_Mode0_raw.Raw' (with path, if needed)
%im_dims = 2 values for image size/resolution (ex. [1920, 780])

    if isnan(im_dims)
        im_dims = [1024, 1024]; % default - based on camera parameters
    end

    % open the file
    fid = fopen(FileName, 'r');
    if fid == -1
      error('Cannot open file: %s', FileName);
    end

    % read imfile as uint16
    im = fread(fid, im_dims, 'uint16');

    % for some reason, it's tilted 90deg, so rotate
    im = (im)';

    % close the file
    fclose(fid);

end