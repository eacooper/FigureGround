function [labels, figs, bgs, dict] = convertrgb2labels(im)
%
% convert rgb values to numeric labels and a figure/ground binary masks
% if input is empty, this just returns the label dictionary

% label dictionary

% rgb values
dict.rgb = {[159, 113, 181], [0, 66, 207],[245, 99, 222],...
    [172, 0, 255],[255, 148, 116],[47, 255, 0],[255, 161, 20],...
    [187, 141, 58],[139, 68, 108],[125, 250, 242],[241, 226, 87],...
    [169, 161, 192],[73, 243, 151],[255, 34, 0],[255, 60, 130]};

% associated categories
dict.name = {'sky','background','plant','bird','small_infrastructure',...
    'pavement','ground','lawn','misc','person','car','bike',...
    'person_on_bike','building','animal'};

% numerical indices for categories
dict.num = 1:15;    % all
dict.fig = 3:15;    % indices of figure regions
dict.bg   = 2;      % indices for background

if ~isempty(im)

    % initialize output matrices
    labels  = zeros(size(im,1),size(im,2));
    figs    = zeros(size(im,1),size(im,2));
    bgs     = zeros(size(im,1),size(im,2));

    % Loop through each category in the dictionary
    for i = 1:numel(dict.rgb)

        % Get the RGB values for the current category
        current_rgb = dict.rgb{i};

        % Find pixels in the input image with the current RGB values
        mask = all(bsxfun(@eq, im, reshape(current_rgb, [1, 1, 3])), 3);

        % Assign the corresponding numerical label to pixels in the mask
        labels(mask) = dict.num(i);

        % If the category is a figure, set the corresponding pixels in figs to 1
        if ismember(dict.num(i), dict.fig)
            figs(mask) = 1;
        end
    end

    % Set background pixels in bgs to 1
    bgs(labels == dict.bg) = 1;

else
    labels = []; figs = []; bgs = [];
end
