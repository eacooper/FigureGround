% Load all data needed for analysis

% Load and combine all motion data

load_dir = fullfile(resultsPath,'motion_perVideo_info');

% Images
UCBim = load(fullfile(load_dir, 'UCB_ims_allTOG.mat'));
display(strcat('loaded images for UCB'));

ZOOim = load(fullfile(load_dir, 'ZOO_ims_allTOG.mat'));
display(strcat('loaded images for ZOO'));

all_set_ims = cat(4, UCBim.all_set_ims, ZOOim.all_set_ims); 
clear UCBim; clear ZOOim;


% OF 
UCBmagori = load(fullfile(load_dir, 'UCB_MagOri_allTOG_skygone.mat'));
display(strcat('loaded optic flow data for UCB'));

ZOOmagori = load(fullfile(load_dir, 'ZOO_MagOri_allTOG_skygone.mat'));
display(strcat('loaded optic flow data for ZOO'));

all_set_mag = cat(3, UCBmagori.all_set_mag, ZOOmagori.all_set_mag);
all_set_ori = cat(3, UCBmagori.all_set_ori, ZOOmagori.all_set_ori); 
clear UCBmagori; clear ZOOmagori;

% Annotations
UCBann = load(fullfile(load_dir, 'UCB_anns_allTOG.mat'));
display(strcat('loaded annotations for UCB'));

ZOOann = load(fullfile(load_dir, 'ZOO_anns_allTOG.mat'));
display(strcat('loaded annotations for ZOO'));

all_set_figs = cat(3, UCBann.all_set_figs, ZOOann.all_set_figs);
all_set_bgs = cat(3, UCBann.all_set_bgs, ZOOann.all_set_bgs);
all_set_border = cat(3, UCBann.all_set_border, ZOOann.all_set_border);
all_set_labels = cat(3, UCBann.all_set_labels, ZOOann.all_set_labels); 
clear UCBann; clear ZOOann;


% Meta Data
UCBmeta = load(fullfile(load_dir, 'UCB_MetaData_allTOG.mat'));
display(strcat('loaded metadata for UCB')); 

ZOOmeta = load(fullfile(load_dir, 'ZOO_MetaData_allTOG.mat'));
display(strcat('loaded metadata for ZOO')); 

all_set_epoque = cat(2, UCBmeta.all_set_epoque, ZOOmeta.all_set_epoque);
all_set_motion_type = cat(2, UCBmeta.all_set_motion_type, ZOOmeta.all_set_motion_type);
all_set_vidnum = cat(2, UCBmeta.all_set_vidnum, ZOOmeta.all_set_vidnum);
all_set_vidset = cat(2, UCBmeta.all_set_vidset, ZOOmeta.all_set_vidset); 
clear UCBmeta; clear ZOOmeta;


% Saliency maps for fixation modelling
load(fullfile(load_dir,'ALL_Saliency_allTOG.mat'));
display(strcat('loaded saliency maps')); 