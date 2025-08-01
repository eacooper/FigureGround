function disparity = calc_disparity(fixation, RF_points, eye_1, eye_2)
% fixation = [f_x, f_y, f_z], position of fixation point
% point    = [p_x, p_y, p_z], positions of point of interest

f1 = fixation - eye_1;  % fixation point, right eye origin
f2 = fixation - eye_2;  % fixation point, left eye origin

p1 = RF_points - eye_1;  % interest point, right eye origin
p2 = RF_points - eye_2;  % interest point, left eye origin

% angular position of both points w.r.t both eyes ()
angle_f1 = atan2d(f1(3), f1(1));
angle_f2 = atan2d(f2(3), f2(1));
angle_p1 = atan2d(p1(:, 3), p1(:, 1));
angle_p2 = atan2d(p2(:, 3), p2(:, 1));

% disparity is the "difference of differences"
beta_1    = angle_f1 - angle_p1;
beta_2    = angle_f2 - angle_p2;

disparity = beta_2 - beta_1;