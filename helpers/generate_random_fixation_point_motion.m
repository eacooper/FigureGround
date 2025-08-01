function fix = generate_random_fixation_point_motion(r_rand,c_rand,row_vals,col_vals,RF_r_pix,scene_mag,scene_ori,scene_saliency,magthreshold)
%
% simulate a random fixation point for computing disparity in visual degrees
% fixation point is constrained so that the RF eccentricity = the RF
% diameter

% get radial distance of all pixels from the RF center
pix_dist    = sqrt((r_rand - row_vals).^2 + (c_rand - col_vals).^2);

% indices of all pixels at this distance where speed is not a NaN (not sky)
% and pixel is salient
dist_inds   = find(pix_dist < RF_r_pix + 1 & pix_dist > RF_r_pix - 1 &...
    ~isnan(scene_mag) & scene_saliency == 1);

if numel(dist_inds) > 0

    % pick one at random
    fix_ind     = randsample(dist_inds,1);

    % grab fixation location info
    [fix.r,fix.c] = ind2sub(size(scene_mag),fix_ind);
    fix.mag         = scene_mag(fix.r,fix.c);
    fix.ori         = scene_ori(fix.r,fix.c);

    % if we've hit a point where speed == magthreshold, we assume
    % stationary fixation
    if fix.mag == magthreshold
        fix.mag         = 0;
        fix.ori         = 0;
    end

else

    fix.r = NaN;
    fix.c = NaN;
    fix.mag = NaN;
    fix.ori = NaN;

end