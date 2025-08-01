% parameters used to process and analyze the DEP dataset

% camera/data format parameters
dims            = [1920 1080];          % matrix dims in pixels
actual_FOV      = 35;                   % horizontal camera FOV
degreePerPixel  = actual_FOV/dims(1);
pixelPerDegree  = dims(1)/actual_FOV;

% buffer around fig/gr border in pixels
bnd_buff = 5;

% define eye positions in meters for disparity calculation
% IPD = 6.2cm. eye1 = left eye
eye_1  = [-0.031, 0, 0];
eye_2  = [0.031, 0, 0];

% quantiles used for normalization
lowerPerc   = 0.01;
upperPerc   = 0.99;

% eccentricities / diameters of RFs to simulate in visual degrees 
% (note = in MT the eccentricity is about equal to the typical RF diameter)
RF_sizes        = [2.5, 5, 10, 15]; 
RF_sizes_pix    = RF_sizes/degreePerPixel;

% indices of scenes that have annotation or alignment quality issues and
% will be excluded
% 22, 26, 51 & 64 -- poor annotations
% 39 & 74 -- poor image/distance map alignment
exclude = [22 26 51 64 39 74];

% colors for visualization
figColour   = [ 228 30 38 ]./255;
bgColour    = [ 51 127 186 ]./255;