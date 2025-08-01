function [hist_ori_norm, max_ind, max_ori] = normalize_ori_data(hist_ori, cntr_ori, max_ind)

%Normalize direction distribution so that the histogram bin with the
%highest number of pixels = bin 1 to get relative direction distribution within-RF

switch nargin %only calculate bin displacement index if not provided
    case 2
        max_ind = find(hist_ori == max(hist_ori));

        % in the case where multiple bins have the same number of pixels
        if numel(max_ind > 1)
            rand_max = randi(length(max_ind));
            max_ind = max_ind(rand_max);
        end

        max_ori = cntr_ori(max_ind);
end

hist_ori_norm = circshift(hist_ori, -(max_ind-1));

end