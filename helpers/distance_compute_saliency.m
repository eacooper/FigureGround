% compute saliency maps for images in distance dataset

% set up paths
dmapPath = fullfile(dataPath, 'DataSets', 'DISTANCE');
outPath = fullfile(resultsPath,'distance_perImage_info');

if ~exist(outPath,'dir')
    mkdir(outPath);
end

% load metadata telling us which scenes to analyze from range dataset
load(fullfile(dmapPath,'dmap_indices.mat'));

% remove indices of scenes that have annotation or alignment quality issues
nlist           = 1:numel(dmap_indices);
nlist           = setdiff(nlist,exclude);
dmap_indices    = dmap_indices(nlist);

% initialize matrices for images, dmaps, and annotations
num_scenes       = numel(dmap_indices);

% for each scene
for x = 1:num_scenes

    display(['loading scene ' num2str(x) '...']);

    % rgb image
    im = imread(fullfile(dmapPath, 'images',['rImage' sprintf('%03d',dmap_indices(x)) '.png']));

    % preprocess image
    im      = double(im);
    im      = im./max(im(:));
    im      = im.^0.45;

    % compute the default bottom-up SUN saliency
    [Saliency_Map, Feature_Maps, ICA_Maps, ~] = Run_SUN(im, []);

    % scale saliency map from 0-1
    Saliency_Map = (Saliency_Map - min(Saliency_Map(:)))/range(Saliency_Map(:));

    % resize to match original image resultion
    Saliency_Map = imresize(Saliency_Map,'OutputSize',size(im(:,:,1)));

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
save(fullfile(outPath,'ALL_Saliency.mat'),'all_set_most_salient');

