%% Load MAT files (in case variables are not in the environment)
load('/Users/vinicius/Documents/MATLAB/for image df/matlab/2025-08-26/20250624_m0055_00003_00122_mcor_timeSeriesFromFijiROIs.mat')

%% Load Files
testFile = '/Users/vinicius/Documents/MATLAB/for image df/20250624_m0055_00005_mcor.tif';
stack = single(read_file(testFile));

%% Compute Relevant Frames

frameOdorOnset = round(odor_onset_s * frame_rate_hz);
frameOdorOffset = round(odor_offset_s * frame_rate_hz);

% Number of frames between onset and offset (including both ends)
frameOdorDuration = frameOdorOffset - frameOdorOnset + 1;

% Baseline starts one window duration before the onset
frameBaselineStart = frameOdorOnset - frameOdorDuration;

% Baseline ends just before the odor onset
frameBaselineEnd = frameOdorOnset - 1;

%% dF/F computation

% Baseline average
baselineAvgImage = mean( ...
    stack(:,:,frameBaselineStart:frameBaselineEnd), ndims(stack) ...
);

odorPresentationAvgImage = mean( ...
    stack(:,:,frameOdorOnset:frameOdorOffset), ndims(stack) ...
);

SBR = (odorPresentationAvgImage - baselineAvgImage) ./ baselineAvgImage;

%% Plot results

figure
lowQuantile = double(quantile(SBR(:), 0.05));
highQuantile = double(quantile(SBR(:), 0.95));
quantiles = [lowQuantile, highQuantile];
imshow(SBR, quantiles);
colormap turbo;
c = colorbar('southoutside');
ylabel(c,'$dF/F$', 'interpreter','latex', 'FontSize', 20);
