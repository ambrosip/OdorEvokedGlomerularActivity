%% Load MAT files (in case variables are not in the environment)
% Use the MAT file from timeSeriesFromFijiROIs
% load('/Users/priscilla/Documents/Local - Moss Lab/20250624/e1/processed/matlab/2025-08-26/20250624_m0055_00003_00122_mcor_timeSeriesFromFijiROIs.mat')

% TO DO
% take quantile of all outcomes together?

%% USER INPUT

% Define percentiles
% LOWER_QUANTILE = 0.0001;
% UPPER_QUANTILE = 0.9999;

% LOWER_QUANTILE = 0.05;
% UPPER_QUANTILE = 0.95;

LOWER_QUANTILE = 0.01;
UPPER_QUANTILE = 0.99;

% Define the colors
max_df_color = [103 0 31] / 255;
min_df_color = [5 48 97] / 255;

%% Extra inputs in case you run this before timeSeriesFromFijiROIs

% % set default firstFig and lastFig boundaries in case user does NOT want a
% % custom subset
% if plotSubset == 0
%     firstAcq = 1;
%     lastAcq = imgsToAnalyze_numberOf;
% end
% firstAcqName = imgsToAnalyzeDirs(1).name;
% lastAcqName = imgsToAnalyzeDirs(lastAcq).name;

%% Creates Diverging Colormap

% Creates 512 color steps between the two colors (512 is arbitrary)
% The number just need to be high enough for the gradient to be smooth
colorpoints = linspace(0.0, 1.0, 512);

% This function is from the following website:
% https://www.kennethmoreland.com/color-maps/
divergingGradient = divergingMap(colorpoints, min_df_color, max_df_color);

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

% Open the first frame of the first file to get image dimensions
filename = imgsToAnalyzeDirs(1).name;
fileDir = imgsToAnalyzeDirs(1).folder;
filepath = fullfile(fileDir, filename);
frameSize = size(read_file(filepath, 1, 1));

% There is one figure for each combination of program and odor
% iFigure is the index of the figure currently being computed
iFigure = 0;

% Iterate through programs
for program_name = string(fieldnames(s_olfactometer))'
    program = s_olfactometer.(program_name);

    % Get program number (number after 'program_')
    programSplit = split(program_name, '_');
    programNumber = str2double(programSplit(2));

    % Skip programs that don't have summaries
    if ~isfield(program, 'summary_by_trial')
        message = ['WARNING: Skipping program %d because ' ...
            'it was marked as ''ignore''.\n'];
        fprintf(message, programNumber)
        continue
    end

    % Get program type
    programType = s_olfactometer.(program_name).type;

    % Get odor and outcome info
    summaryByTrial = program.summary_by_trial;
    programOdors = unique(summaryByTrial.odor);

    % Iterating through odors in the program
    for odor = programOdors'
        % Get table with data only for the current odor
        odorRows = summaryByTrial.odor == odor;
        odorTable = summaryByTrial(odorRows, :);

        % Initialize variables that contain the current average images
        % ASSUMPTION: There is only one baseline for all outcomes!
        avgBaseline = zeros(frameSize);
        avgHit = zeros(frameSize);
        avgMiss = zeros(frameSize);
        avgFalse = zeros(frameSize);
        avgNa = zeros(frameSize);

        % Initialize counts of how many images where used in the average
        nBaseline = 0;
        nHit = 0;
        nMiss = 0;
        nFalse = 0;
        nNa = 0;

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
                fprintf('WARNING: File %s not found.\n', filename)
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
                case 'hit'
                    nHit = nHit + 1;
                    avgHit = ...
                        avgHit * (nHit - 1) / nHit + ...
                        mean(frames, ndims(frames)) / nHit;

                case 'miss'
                    nMiss = nMiss + 1;
                    avgMiss = ...
                        avgMiss * (nMiss - 1) / nMiss + ...
                        mean(frames, ndims(frames)) / nMiss;

                case 'false choice'
                    nFalse = nFalse + 1;
                    avgFalse = ...
                        avgFalse * (nFalse - 1) / nFalse + ...
                        mean(frames, ndims(frames)) / nFalse;

                case 'na'
                    nNa = nNa + 1;
                    avgNa = ...
                        avgNa * (nNa - 1) / nNa + ...
                        mean(frames, ndims(frames)) / nNa;

            end
        end

        % Don't compute the ratio if there are no acquisition files
        % If there were no acquisition files nBaseline is zero
        if nBaseline == 0
            message = ...
                'WARNING: No files for program %d (%s) and odor %d.\n';
            fprintf(message, programNumber, programType, odor);
            continue
        end

        message = 'Analyzed %d files for program %d (%s) and odor %d.\n';
        fprintf(message, nBaseline, programNumber, programType, odor);

        % Otherwise, there is one more figure to plot
        iFigure = iFigure + 1;

        % Add current program and odor to current figure
        % figures(iFigure).program = programNumber;
        figures(iFigure).type = programType;
        figures(iFigure).odor = odor;

        % Computes Signal to Baseline Ratio of all outcomes. If there are
        % no instances of that outcome the array will be filled with zeros.
        figures(iFigure).hit.image = ...
            (nHit > 0) * (avgHit - avgBaseline) ./ avgBaseline;
        figures(iFigure).miss.image = ...
            (nMiss > 0) * (avgMiss - avgBaseline) ./ avgBaseline;
        figures(iFigure).false.image = ...
            (nFalse > 0) * (avgFalse - avgBaseline) ./ avgBaseline;
        figures(iFigure).na.image = ...
            (nNa > 0) * (avgNa - avgBaseline) ./ avgBaseline;

        % We will keep track of how many instances of each outcome we got
        figures(iFigure).totalAcquisitions = nBaseline;
        figures(iFigure).hit.totalNumber = nHit;
        figures(iFigure).miss.totalNumber = nMiss;
        figures(iFigure).false.totalNumber = nFalse;
        figures(iFigure).na.totalNumber = nNa;
    end
end

%% Plot the Figures

% Extreme values of the signal-to-baseline ratio can be so high that it
% is hard to see anything else in the figure. Thus, we use quantiles to
% exclude the outliers. The specific values (5% and 95%) are arbitrary.

% If no files were analyzed stop computation
if isempty(figures)
    error('No figures to plot.')
end

% Initiliaze limits
lowerLimit = NaN;
upperLimit = NaN;

% Compute limits figure by figure. The number lowerLimit will be the
% smallest of all the 5% quantiles and upperLimit will be the largest of
% all the 95% quantiles (if LOWER_QUANTILE is set to 0.05 and
% UPPER_QUANTILE is set to 0.95)
for i = 1:length(figures)
    lowerLimit = min([ ...
        lowerLimit, ...
        quantile(figures(i).hit.image(:), LOWER_QUANTILE), ...
        quantile(figures(i).miss.image(:), LOWER_QUANTILE), ...
        quantile(figures(i).false.image(:), LOWER_QUANTILE), ...
        quantile(figures(i).na.image(:), LOWER_QUANTILE) ...
    ]);

    upperLimit = max([ ...
        upperLimit, ...
        quantile(figures(i).hit.image(:), UPPER_QUANTILE), ...
        quantile(figures(i).miss.image(:), UPPER_QUANTILE), ...
        quantile(figures(i).false.image(:), UPPER_QUANTILE), ...
        quantile(figures(i).na.image(:), UPPER_QUANTILE) ...
    ]);
end

% Take the limit that is larger in absolute value
% It also turns it into a double because that is a requirement for plots
absoluteLimit = double(max(abs(upperLimit), abs(lowerLimit)));

% Make range symmetrical around zero (so white means zero in the plots).
% The same range will be used for all figures, for easy comparisons.
plotRange = [-absoluteLimit absoluteLimit];

for iFigure = 1:length(figures)
    % Get figure name
    % Get the start of firstAcqName (before the third underline)
    figNameStart = split(string(firstAcqName), '_');
    figNameStart = join(figNameStart(1:end-1), '_');

    % Get the number of the last acquisition
    figNameMiddle = split(string(lastAcqName), '_');
    figNameMiddle = figNameMiddle(end-1);

    % Reformats program type string from "Fine 1" to "fine_1", for example
    figNameEnd = join(split(lower(figures(iFigure).type)), '_');

    % Join results in the correct format
    figName = sprintf("%s_to_%s_odor_%d_%s", figNameStart, ...
        figNameMiddle, figures(iFigure).odor, figNameEnd);

    fig = figure('Name', figName);

    % 3 possible outcomes in a row (or just 1 if there is NA)
    tl = tiledlayout('horizontal');
    title(tl, figName, 'Interpreter', 'none');

    if figures(iFigure).na.totalNumber == 0
        nPlots = 3;

        nexttile
        imshow(figures(iFigure).hit.image, plotRange)
        title('Hits', 'FontSize', 16)

        nexttile
        imshow(figures(iFigure).false.image, plotRange)
        title('False Choices', 'FontSize', 16)

        nexttile
        imshow(figures(iFigure).miss.image, plotRange)
        title('Misses', 'FontSize', 16)
    else
        nPlots = 1;

        nexttile
        imshow(figures(iFigure).na.image, plotRange)
        title('NA', 'FontSize', 16)
    end

    baseSize = 600;
    fig.Position = [200 100 baseSize * nPlots baseSize + 150];

    % Uses the colormap created at the start
    cb = colorbar;
    cb.Layout.Tile = 'south';

    colormap(divergingGradient);
    clim(plotRange);

    % ylabel(cb, '$dF/F$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel(cb, 'dF/F', 'FontSize', 16);

    tl.TileSpacing = 'compact';
    tl.Padding = 'compact';

    drawnow;
end