%% Load MAT files (in case variables are not in the environment)
% Use the MAT file from timeSeriesFromFijiROIs
load(['../for image df/matlab/2025-08-26/' ... 
    '20250624_m0055_00003_00122_mcor_timeSeriesFromFijiROIs.mat'])

%% Creates Diverging Colormap
% Define the colors
REDDISH_PURPLE = [.8000 .4745 .6549];
BLUE           = [    0 .4471 .6980];

% Creates 512 color steps between the two colors (512 is arbitrary)
% The number just need to be high enough for the gradient to be smooth
colorpoints = linspace(0.0, 1.0, 512);

% This function is from the following website:
% https://www.kennethmoreland.com/color-maps/
divergingGradient = divergingMap(colorpoints, BLUE, REDDISH_PURPLE);

%% Compute Relevant Frames

% ASSUMPTIONS: 
%   1 - Baseline and Signal intervals have the same number of frames
%   2 - Signal interval is the odor presentation interval

frameOdorOnset = round(odor_onset_s * frame_rate_hz);
frameOdorOffset = round(odor_offset_s * frame_rate_hz);

% Number of frames between onset and offset (including both ends)
frameOdorDuration = frameOdorOffset - frameOdorOnset + 1;

% Baseline starts one window duration before the onset
frameBaselineStart = frameOdorOnset - frameOdorDuration;

%% Iterate over Data

nPlots = 0;

% Iterate through programs
for program_name = string(fieldnames(s_olfactometer))'
    program = s_olfactometer.(program_name);

    % Skip programs that don't have summaries
    if ~isfield(program, 'summary_by_trial')
        continue
    end

    % Get program number 
    programSplit = split(program_name, '_');
    programNumber = str2double(programSplit(2));
        
    summaryByTrial = program.summary_by_trial;
    programOdors = unique(summaryByTrial.odor);
    
    % Iterating through odors in the program
    for odor = programOdors'
        % Select and iterate through acqIdx corresponding to the odor
        odorRows = summaryByTrial.odor == odor;

        % Initialize variables that contain the current average images
        % MATLAB will expand those to matrices of the correct size later
        baselineAvgImage = 0;
        signalAvgImage = 0;

        % nImages counts how many images where used in the average
        nImages = 0;

        for acqIdx = summaryByTrial(odorRows, :).acqIdx'
            % Skip trials without acquisition
            if isnan(acqIdx)
                continue
            end

            filename = imgsToAnalyzeDirs(acqIdx).name;
            fileDir = imgsToAnalyzeDirs(acqIdx).folder;
            filepath = fullfile(fileDir, filename);
            
            % Skips computation and prints warning if file not found
            if ~isfile(filepath)
                % fprintf("WARNING: File %s not found.\n", filename)
                continue
            end

            % Since there is a file, we increase the counter
            nImages = nImages + 1;

            % Load and computes the mean of the relevant frames      
            % We use a trick to compute the average on the fly:
            % https://stackoverflow.com/a/23493727
            baselineFrames = single(read_file( ...
                filepath, frameBaselineStart, frameOdorDuration));
            baselineAvgImage = baselineAvgImage * (nImages-1)/nImages + ...
                mean(baselineFrames, ndims(baselineFrames)) / nImages;

            signalFrames = single(read_file( ...
                filepath, frameOdorOnset, frameOdorDuration));
            signalAvgImage = signalAvgImage * (nImages-1)/nImages + ...
                mean(signalFrames, ndims(baselineFrames)) / nImages;
        end

        % Don't plot the image if there are no acquisition files
        if nImages == 0
            continue
        end

        % Computes Signal to Baseline Ratio
        signalBaselineRatioImage = ...
            (signalAvgImage - baselineAvgImage) ./ baselineAvgImage;

        % There is one more SBR image to plot
        nPlots = nPlots + 1;

        % Add current program, odor, and SBR image to plots
        plots(nPlots).program = programNumber;
        plots(nPlots).odor = odor;
        plots(nPlots).image = signalBaselineRatioImage;
    end
end

%% Plots Data on SBRPlots

programs = sort(unique([plots.program]));
odors = sort(unique([plots.odor]));

nPrograms = length(programs);
nOdors = length(odors);

% Find the position and quantiles for each image
for i = 1:length(plots)
    row = find(programs == plots(i).program, 1);
    column = find(odors == plots(i).odor, 1);

    plots(i).position = column + nOdors * (row - 1);

    % Extreme values of the signal-to-baseline ratio can be so high that it
    % is hard to see anything else in the figure. Thus, we use quantiles to
    % exclude the outliers. The specific values (5% and 95%) are arbitrary.
    plots(i).lowQuantile = ...
        double(quantile(plots(i).image(:), 0.05));
    plots(i).highQuantile = ...
        double(quantile(plots(i).image(:), 0.95));    
end

% Make limits uniform across plots
lowerLimit = min([plots.lowQuantile]);
upperLimit = max([plots.highQuantile]);
plotRange = [lowerLimit upperLimit];

% Make range symmetrical so white means 0
absoluteLimit = max(abs(plotRange));
plotRange = [-absoluteLimit absoluteLimit];

% Start figure and make it big!
fig = figure;
set(fig, 'Position', [200 100 1300 800]);

for i = 1:length(plots)
    subplot(nPrograms, nOdors, plots(i).position);
    imshow(plots(i).image, plotRange);

    t = sprintf("Program %d Odor %d", plots(i).program, plots(i).odor);
    title(t, 'interpreter','latex', 'FontSize', 20)
end

% Creating ax here makes sure the colorbar below encopasses all plots and
% not just the last subplot in the for loop above.
ax = axes(fig);
t = 'Average Signal-to-Baseline Ratio Plots';
title(ax, t, 'interpreter','latex', 'FontSize', 24);

axis off;
axis tight;

% Uses the colormap create at the start
colormap(divergingGradient);
clim(plotRange);

c = colorbar('southoutside');
ylabel(c,'$dF/F$', 'interpreter','latex', 'FontSize', 20);

drawnow;

%% Save figure
savePath = 'SignalBaselineRatio.png';
saveas(fig, savePath);