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

% There is one figure for each combination of program and odor
% iFigure is the index of the figure currently being computed
iFigure = 0;

% Iterate through programs
for program_name = string(fieldnames(s_olfactometer))'
    program = s_olfactometer.(program_name);

    % Skip programs that don't have summaries
    if ~isfield(program, 'summary_by_trial')
        continue
    end

    % Get program number (number after 'program_')
    programSplit = split(program_name, '_');
    programNumber = str2double(programSplit(2));
        
    summaryByTrial = program.summary_by_trial;
    programOdors = unique(summaryByTrial.odor);
    
    % Iterating through odors in the program
    for odor = programOdors'
        % Get table with data only for the current odor
        odorRows = summaryByTrial.odor == odor;
        odorTable = summaryByTrial(odorRows, :);

        % Initialize variables that contain the current average images
        % MATLAB will expand those to matrices of the correct size later
        % ASSUMPTION: There is only one baseline for all outcomes!
        avgBaseline = 0;
        avgHit = 0;
        avgMiss = 0;
        avgFalse = 0;

        % Initialize counts of how many images where used in the average
        nBaseline = 0;
        nHit = 0;
        nMiss = 0;
        nFalse = 0;

        for row = 1:height(odorTable)
            acqIdx = odorTable{row, 'acqIdx'};

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

            % Since there is a file, we increase the baseline counter
            nBaseline = nBaseline + 1;

            % Load and computes the mean of the relevant frames      
            % We use a trick to compute the average on the fly:
            % https://stackoverflow.com/a/23493727
            frames = single(read_file( ...
                filepath, frameBaselineStart, frameOdorDuration));

            avgBaseline = ...
                avgBaseline * (nBaseline-1) / nBaseline + ...
                mean(frames, ndims(frames)) / nBaseline;
            
            % Do the same as above, but for the corresponding outcome
            % The array 'frames' is the same for all possible outcomes
            frames = single(read_file( ...
                        filepath, frameOdorOnset, frameOdorDuration));

            switch odorTable{row, 'outcome'}
                case "hit"
                    nHit = nHit + 1;
                    avgHit = ...
                        avgHit * (nHit - 1) / nHit + ...
                        mean(frames, ndims(frames)) / nHit;
                
                case "miss"
                    nMiss = nMiss + 1;
                    avgMiss = ...
                        avgMiss * (nMiss - 1) / nMiss + ...
                        mean(frames, ndims(frames)) / nMiss;

                case "false choice"
                    nFalse = nFalse + 1;
                    avgFalse = ...
                        avgFalse * (nFalse - 1) / nFalse + ...
                        mean(frames, ndims(frames)) / nFalse;
            
            end
        end

        % Don't compute the ratio if there are no acquisition files
        % If there were no acquisition files nBaseline is zero
        if nBaseline == 0
            continue
        end

        % Otherwise, there is one more figure to plot
        iFigure = iFigure + 1;

        % Add current program and odor to current figure
        figures(iFigure).program = programNumber;
        figures(iFigure).odor = odor;

        % Computes Signal to Baseline Ratio of all outcomes
        % If there are no instances of the outcome, make it equal to NaN
        % This will make the code for the quantiles simpler later
        if nHit > 0
            figures(iFigure).hit = (avgHit-avgBaseline) ./ avgBaseline;
        else
            figures(iFigure).hit = NaN;
        end

        if nMiss > 0
            figures(iFigure).miss = (avgMiss-avgBaseline) ./ avgBaseline;
        else
            figures(iFigure).miss = NaN;
        end

        if nFalse > 0
            figures(iFigure).false = (avgFalse-avgBaseline) ./ avgBaseline;
        else
            figures(iFigure).false = NaN;
        end
    end
end

%% Plot the Figures

% Extreme values of the signal-to-baseline ratio can be so high that it
% is hard to see anything else in the figure. Thus, we use quantiles to
% exclude the outliers. The specific values (5% and 95%) are arbitrary.

% The MATLAB function quantile(

% Initiliaze limits
lowerLimit = NaN;
upperLimit = NaN;

% Compute limits figure by figure. The number lowerLimit will be the
% smallest of all the 5% quantiles and upperLimit will be the largest of
% all the 95% quantiles.

for i = 1:length(figures)    
    lowerLimit = min([ ...
        lowerLimit, ...
        quantile(figures(i).hit(:), 0.05), ...
        quantile(figures(i).miss(:), 0.05), ...
        quantile(figures(i).false(:), 0.05) ...
    ]);

    upperLimit = max([ ...
        upperLimit, ...
        quantile(figures(i).hit(:), 0.95), ...
        quantile(figures(i).miss(:), 0.95), ...
        quantile(figures(i).false(:), 0.95) ...
    ]);   
end

% Take the limit that is larger in absolute value
% It also turns it into a double because that is a requirement for plots
absoluteLimit = double(max(abs(upperLimit), abs(lowerLimit)));

% Make range symmetrical around zero (so white means zero in the plots).
% The same range will be used for all figures, for easy comparisons.
plotRange = [-absoluteLimit absoluteLimit];

for iFigure = 1:length(figures)
    fig = figure('Name', 'test_figure');

    % 3 possible outcomes in a row
    tl = tiledlayout("horizontal");

    nPlots = 0;

    if ~isnan(figures(iFigure).hit)
        nexttile
        imshow(figures(iFigure).hit, plotRange)
        title('Hits', 'FontSize', 16)
        nPlots = nPlots + 1;
    end

    if ~isnan(figures(iFigure).miss)
        nexttile
        imshow(figures(iFigure).miss, plotRange)
        title('Misses', 'FontSize', 16)
        nPlots = nPlots + 1;
    end

    if ~isnan(figures(iFigure).false)
        nexttile
        imshow(figures(iFigure).false, plotRange)
        title('False Choices', 'FontSize', 16)
        nPlots = nPlots + 1;
    end

    baseSize = 600;
    fig.Position = [200 100 baseSize * nPlots baseSize + 150];

    % Uses the colormap created at the start
    cb = colorbar;
    cb.Layout.Tile = "south";

    colormap(divergingGradient);
    clim(plotRange);
    
    ylabel(cb,'$dF/F$', 'interpreter', 'latex', 'FontSize', 16);

    tl.TileSpacing = "compact";
    tl.Padding = "compact";
    
    drawnow;
end