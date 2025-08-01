function [mags, oris] = calc_mag_oris(Vx, Vy)
%
% get speed in pixels per frame and motion direction from the estimation
% X/Y motion

    mags = zeros(size(Vx));
    oris = zeros(size(Vx));

     for n = 1:size(Vx, 3)
    
        opflow = opticalFlow(Vx(:, :, n),Vy(:, :, n));
        mags(:, :, n) = sqrt(Vx(:, :, n).^2 + Vy(:, :, n).^2);%.*(FPS*actual_FOV_deg/horiz);
        oris(:, :, n) = opflow.Orientation;

     end

end