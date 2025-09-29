% set up path
saveDir      = fullfile(resultsPath,'motion_perVideo_info');

% load and concatenate
UCBim = load(fullfile(saveDir, 'UCB_ims_allTOG.mat'));
display(strcat('loaded images for UCB'));

ZOOim = load(fullfile(saveDir, 'ZOO_ims_allTOG.mat'));
display(strcat('loaded images for ZOO'));

all_set_ims = cat(4, UCBim.all_set_ims, ZOOim.all_set_ims); clear UCBim; clear ZOOim;

% for each image
for x = 1:size(all_set_ims,4)

    % load image and gamma correct (secnd frame in dim 3 is target image
    im = (double(all_set_ims(:, :, 2, x))./255).^(1/gamma_val);

    % compute the default bottom-up SUN saliency
    [Saliency_Map, Feature_Maps, ICA_Maps, ~] = Run_SUN(im, []);

    % scale saliency map from 0-1
    Saliency_Map = (Saliency_Map - min(Saliency_Map(:)))/range(Saliency_Map(:));

    % resize to match original image resultion
    Saliency_Map = imresize(Saliency_Map,'OutputSize',size(im));

    % store in matrix
    all_set_saliency(:,:,x) = Saliency_Map;

    % find 50% most salient pixels
    median_sal = median(Saliency_Map(:));

    most_salient = Saliency_Map;
    most_salient(Saliency_Map <= median_sal) = 0;
    most_salient(Saliency_Map > median_sal) = 1;

    % store these as well
    all_set_most_salient(:,:,x) = most_salient;

end


% Save the saliency info
saveName = fullfile(saveDir, 'ALL_Saliency_allTOG.mat');
save(saveName, '-v7.3', 'all_set_saliency','all_set_most_salient');

