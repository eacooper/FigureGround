% parameters used to process and analyze the motion dataset

% names for each recording site
vidset_IDs  = {'UCB','ZOO'};

% number of videos from each recording site
vid_nums = [24, 25];

% camera parameters
dims            = [1024, 1024]; % pixels
horiz           = dims(1);
FPS             = 120;
actual_FOV_deg  = 54*horiz/1280;
degreePerPixel  = actual_FOV_deg/horiz;
pixelPerDegree  = horiz/actual_FOV_deg;
gamma_val       = 2.2;

% based on noise analysis, we only consider speeds that are at least 0.5 deg/sec
magthreshold = 0.5;

% quantiles used for normalization
lowerPerc           = 0.01;
upperPerc           = 0.99;

% OF calculation parameters - for COARSE passthrough
taps_coarse         = 2;
scale_coarse        = 0.125;
pyr_levels_coarse   = 1;
yesUndistort_coarse = false;

% OF calculation parameters - for FINAL calculation on selected frames
taps            = 3;
scale           = 1;
pyr_levels      = 6;
yesUndistort    = true;

% Epoque selection parameters
perc_cutoff     = 0.5; 
window_size     = 120; % in number of frames over which to smooth the mean speed data with a Gaussian
num_epoques     = 2;

secTimeSep          = 3; %seconds
epoqueTimeLength    = secTimeSep-0.5;
secFrameSep         = FPS*secTimeSep;
changeNeighbourhood = FPS*epoqueTimeLength;

% the start and stop of motion must be at least this number of frames before or after the peak
central_buff = 20;

% Number of frames to copy over
num_copy_frames = 11;

% eccentricities / diameters of RFs to simulate in visual degrees 
% (note = in MT the eccentricity is about equal to the typical RF diameter)
RF_sizes        = [2.5, 5, 10, 15];
RF_sizes_pix    = RF_sizes/degreePerPixel;

% buffers
edge_buff   = 5;   % buffer around image border in pixels
bnd_buff    = 5;   % buffer around fig/gr border in pixels

% Excluded frames:
UCB_ex = [6 12];
% UCB6: Too dark, motion noise
% UCB12: Motion overflow from figure to ground
ZOO_ex = [7 10 11 17]; 
% ZOO7: Poor annotation
% ZOO10: Too dark, motion noise
% ZOO11: motion in lower visual field due to camera moving
% ZOO17: Too dark, motion noise

ALL_ex = {UCB_ex, ZOO_ex};

% colors for visualization
figColour   = [ 228 30 38 ]./255;
bgColour    = [ 51 127 186 ]./255;
