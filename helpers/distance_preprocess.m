% compile images, distance maps, and annotations into a single mat file
close all;

% set paths
dmapPath = fullfile(dataPath, 'DataSets', 'DISTANCE');
annPath = fullfile(dataPath, 'Annotations', 'DISTANCE');
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
dims             = [1920 1080];

for n = 1:num_scenes

    display(['loading scene ' num2str(n) '...']);

    %% load 

    % rgb image
    im = imread(fullfile(dmapPath, 'images',['rImage' sprintf('%03d',dmap_indices(n)) '.png']));
    
    % distance map
    dmap = load(fullfile(dmapPath, 'distancemaps',['rRange' sprintf('%03d',dmap_indices(n)) '.mat']));
    
    % annotations
    findAnn = dir(fullfile(annPath, ['_depth_seg_masks_*_' sprintf('%02d',nlist(n)) '_seg.png']));
    ann = imread(fullfile(annPath, findAnn.name));
    
    %% image

    % preprocess image
    im      = double(im);
    im      = im./max(im(:));
    im      = im.^0.45;
    imgr    = rgb2gray(im);

    %% distance

    % preprocess distance map
    X        = dmap.rangeMap(:,:,2); % - left, + right
    Y        = dmap.rangeMap(:,:,3); % - down, + up
    Z        = dmap.rangeMap(:,:,1); % - away
    distance = dmap.rangeImg;        % This is a matrix of the euclidian distance to every pixel of the image
    
    % set sky pixels and invalid measurements to NaNs
    X(Z > -2)           = NaN; 
    Y(Z > -2)           = NaN;
    distance(Z > -2)    = NaN;
    Z(Z > -2)           = NaN;

    % make away Z = positive
    Z = -Z; 

    %% labels

    % preprocess labels

    % convert RGB to categorical labels
    [labels, figs, bgs, dict] = convertrgb2labels(ann);

    % augment sky labels with information from the distance maps
    labels(isnan(Z))    = 1;
    figs(isnan(Z))      = 0;
    bgs(isnan(Z))       = 0;

    % make a boundary mask
    bnd         = boundarymask(labels);             % boundary between labelled regions
    bnd_f       = fspecial('average',bnd_buff);   % filter to thicken boundary
    bnd_thick   = imfilter(double(bnd),bnd_f);% apply filter
    bnd_thick(bnd_thick > 0) = 1;

    % apply boundary mask to fig/ground masks
    figs(bnd_thick > 0) = 0;
    bgs(bnd_thick > 0)  = 0;

    % combine for visual
    fgbg = figs;
    fgbg(bgs == 1) = 0.5;

    %% store for analysis
    all_set_ims(:,:,n)          = imgr;
    all_set_x(:,:,n)            = X;
    all_set_y(:,:,n)            = Y;
    all_set_z(:,:,n)            = Z;
    all_set_distance(:,:,n)     = distance;
    all_set_labels(:,:,n)       = labels;
    all_set_figs(:,:,n)         = figs;
    all_set_bgs(:,:,n)          = bgs;
    all_set_border(:,:,n)       = bnd;

end

% save files for next analysis step
save(fullfile(outPath,'ALL_DistancesImagesAnnotations.mat'),...
    'dims','num_scenes','dict','all_set_ims','all_set_x','all_set_y','all_set_z','all_set_distance',...
    'all_set_labels','all_set_figs','all_set_bgs','all_set_border')

