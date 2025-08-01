clear all; close all;

%% Step 1: set up paths & parameters

%%%%%%%%%% EDIT THESE BASED ON YOUR SYSTEM & RUN BEFORE STARTING ANALYSIS %%%%%%%%%%%%%%
% working directory containing this code
codePath = '.';

% add helper functions to your path
addpath(genpath(fullfile(codePath,'helpers')));

% path to location of DataSets and Annotations folders downloaded from Zenodo Repository
dataPath = '/Users/emily/Desktop/FGRepo';

% path where you'd like to store interim analysis files and results
resultsPath = fullfile(codePath,'results');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set analysis parameters
distance_set_analysis_params;

%% Step 2: preprocess distance map data
% Note: you can comment this section out and just run:
% distance_load_files if you just want to do Step 3

% organize distance maps and images
distance_preprocess;

% compute saliency maps for each image
distance_compute_saliency;

% 
% %% preprocess and plot all the distance map data and annotations
% % can also just load preprocessed data
% if(0);      distance_preprocess_and_plot;
% else;       load(strjoin({userPath,'Code','interim_analysis_files','alldata_DEP.mat'},pathJoiner));
% end
% 
% %% compute saliency maps for each image
% if(0)
%     distance_compute_saliency;
% end

%% Step 3: data analysis

% compute whole scene histograms of absolute distance, make plots and run stats
distance_whole_scene_analysis;

% find random RFs, plot each one, and aggregate results
distance_getRFs;

% plot figures and run statistical tests for these RFs
% Note: does not require running distance_load_files so long as
% distance_getRFs has been run
distance_analyzeRFs;


