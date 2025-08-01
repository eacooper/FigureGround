function fix = generate_random_fixation_point_distance(r_rand,c_rand,row_vals,col_vals,RF_r_pix,scene_x,scene_y,scene_z,scene_distance,scene_saliency)
%
% simulate a random fixation point for computing disparity in visual degrees
% fixation point is constrained so that the RF eccentricity = the RF
% diameter

% get radial distance of all pixels from the RF center
pix_dist    = sqrt((r_rand - row_vals).^2 + (c_rand - col_vals).^2);

% indices of all pixels at this distance that aren't NaNs
dist_inds   = find(pix_dist < RF_r_pix + 1 & pix_dist > RF_r_pix - 1 & ...
    ~isnan(scene_distance) & scene_saliency == 1);

if numel(dist_inds) > 0

    % pick one at random
    fix_ind     = randsample(dist_inds,1);

    % grab fixation location info
    [fix.r,fix.c] = ind2sub(size(scene_distance),fix_ind);
    fix.x         = scene_x(fix.r,fix.c);
    fix.y         = scene_y(fix.r,fix.c);
    fix.z         = scene_z(fix.r,fix.c);

else

    fix.r = NaN;
    fix.c = NaN;
    fix.x = NaN;
    fix.y = NaN;
    fix.z = NaN;
end