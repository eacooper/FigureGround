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
motion_set_analysis_params;

%% Step 2: compute and clean up optic flow
% Note: you can comment this section out and just run:
% motion_load_files if you just want to do Step 3

% compute high quality optic flow maps and store them for analysis
% Note: this step may take a few hours
motion_calc_final_OF;

% compile together annotations, preprocess/save masks
motion_annotations_compiler;

% compute saliency maps for each frame
motion_compute_saliency;

%% Step 3: data analysis

% compute whole scene histograms of absolute speed and direction
motion_whole_scene_analysis;

% find random RFs, plot each one, and aggregate results
motion_getRFs;

% plot figures and run statistical tests for these RFs
% Note: does not require running motion_load_files so long as
% motion_getRFs has been run
motion_analyzeRFs;
