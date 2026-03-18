%{
DOCUMENTATION
written by VA on Mar/2026

GOAL:
    Calculate the z-scores of the dF/F for each ROI across all
    acquisitions.

DEPENDS on:
    mat file created by timeSeriesFromFijiROIs
    (it only uses timing information + the struct s + db_trials)
%}

%% USER INPUT

% Paths to the experiment folders to be included in the plots
% ASSUMPTION: folders are ordered chronologically
expFolders = [ ...
    "/Users/priscilla/Documents/Local - Moss Lab/20251007/sid260/e1", ...
    "/Users/priscilla/Documents/Local - Moss Lab/20251007/sid260/e2"
    ];

%% Load relevant .mat files

nFolders = length(expFolders);

% Create an struct array of experiment data/metadata
% Since there are only a few folders, we do it directly inside the loop

for iFolder = 1:nFolders
    % ALERT: loadData is defined at the end of the file.
    expData(iFolder) = loadData(expFolders(iFolder));
end


%% Compute Z-score

% Get baselines from the first acquisition on the first experiment
baselineAvgs = expData(1).baselineAvgs;

% ASSUMPTION: all experiments must have the same number of ROIs

% Add dF/F and z-scores to each expData
for iFolder = 1:nFolders
    % ALERT: computeZScores is defined at the end of the file.
    [dFF, zScore] = ...
        computeZScores(expData(iFolder), baselineAvgs);

    expData(iFolder).dFF = dFF;
    expData(iFolder).zScore = zScore;
end

%% Plot Z-Score

% Find z-scores range
zMin = Inf;
zMax = -Inf;

for data = expData
    zMin = min(zMin, min(data.zScore(:)));
    zMax = max(zMax, max(data.zScore(:)));
end

% Compute the time range for the plot

% We shift the time range for the first experiment to start at 0.
% ASSUMPTION: experiments are in chronological order

startTime = expData(1).startTime;
endTime = expData(end).endTime;
plotTimeRange = [0 (endTime - startTime)];

% ASSUMPTION: All experiments have the same number of ROIs
nROIs = expData(1).nROIs;

for iROI = 1:nROIs
    % Name figure
    figName = strcat(expData(1).name, ...
        '_to_', expData(end).name, ...
        '_ROI_', int2str(iROI));

    fig = figure('Name', figName);
    title(figName, 'Interpreter', 'none');

    % Plot acquisitions with different colors
    hold on;

    for data = expData
        for iFiles = 1:data.nFiles
            % Compute x range for each file (in minutes)
            signalStart = data.signalStartTimes(iFiles) - startTime;
            signalEnd = signalStart + data.signalDuration;

            timeRange = linspace(signalStart, signalEnd, data.nFrames);

            % Highlight odor presentation
            odorStart = data.odorStartTimes(iFiles) - startTime;
            odorEnd = odorStart + data.odorDuration;

            patch( ...
                [odorStart odorEnd odorEnd odorStart], ...
                [zMin zMin zMax zMax], ...
                [0.9 0.9 0.9], ...
                'EdgeColor', 'none');

            % Plot files
            zScoreROI = data.zScore(:, iFiles, iROI);
            plot(timeRange, zScoreROI(:));

        end
    end

    hold off;

    % Set plot ranges
    xlim(plotTimeRange);
    ylim([zMin zMax]);

    % Format the x-axis to show minutes and seconds
    xtickformat('mm:ss');
    xlabel('Time (min:sec)')

    ylabel('Z-Score for dF/F')

    % Distribute figures across the screen

    % ASSUMPTION: There are less than 16 ROIs (it still works with more,
    %             but they will be clumped in the right side of the screen)

    % ALERT: ROIs are arranged from bottom to top and from left to right.
    figWidth = 350;
    figHeight = 150;
    vPadding = 80;
    hPadding = 10;

    fig.Position = [ ...
        (50 + floor((iROI-1) / 4) * (figWidth + hPadding)) ...
        (20 + mod(iROI - 1, 4) * (figHeight + vPadding)) ...
        figWidth ...
        figHeight];

    drawnow;
end

%%  Auxiliary Functions

function expData = loadData(expFolder)
% LOADDATA Load and compute all data needed to compute the z-score

% Find the outputs of the timeSeriesFromFijiROIs script
matPaths = dir(fullfile( ...
    expFolder, '**', '*timeSeriesFromFijiROIs.mat'));

% There must be at least one output
if isempty(matPaths)
    error(['Did you run timeSeriesFromFijiROIs before running ' ...
        'this script?\nCould not find its .mat output inside %s'], ...
        expFolder);
end

% Pick the .mat file with the latest date
[~, maxDateIndex] = max([matPaths.datenum]);
matPath = matPaths(maxDateIndex);
matPath = fullfile(matPath.folder, matPath.name);

% Load relevant variables from .mat file
load(matPath, ...
    's', 'db_trials', 'firstAcqName', 'currentImgInfo', ...
    'photobleaching_window_frames', 'photobleaching_window_s', ...
    'frame_rate_hz', 'odor_onset_s', 'odor_dur_s');

% Create output struct

% Get experiment name (start of TIFF file names)
nameStart = split(string(firstAcqName), '_');
expData.name = join(nameStart(1:end-2), '_');

% Get first file signal
files = string(fieldnames(s));
firstFile = files(1);
firstSignal = s.(firstFile);

% Get baseline
baselineStart = ceil(photobleaching_window_frames);
baselineEnd = floor(odor_onset_s * frame_rate_hz);
expData.baselineAvgs = ...
    mean(firstSignal(baselineStart:baselineEnd, :), 1);

% Store sizes
expData.nFrames = size(firstSignal(baselineStart:end, :), 1);
expData.nFiles = length(files);
expData.nROIs = size(firstSignal, 2);

% ALERT: We remove the photobleaching window in all acquisitions, this
%        is why the nFrames above is not just size(firstSignal, 1).

% Initialize signals
expData.signals(expData.nFrames, expData.nFiles, expData.nROIs) = 0;

for iFile = 1:expData.nFiles
    % Get signal from the stored data
    signal = s.(files(iFile));

    % Crops signal to exclude photobleaching window
    signal = signal(baselineStart:end, :);

    % Group all signals together for easy processing
    expData.signals(:, iFile, :) = signal;
end

% Get timing data

% Get start time for the whole experiment
startTimeStr = extractBetween( ...
    currentImgInfo(1).ImageDescription,"epoch = [","]");
expData.startTime = ...
    datetime(startTimeStr, 'InputFormat', 'yyyy M d H m s.SSS');

% Get start/end times for each acquisition (in minutes)
expData.signalStartTimes = expData.startTime + minutes( ...
    db_trials.trial_locs_min + photobleaching_window_s / 60);
expData.signalDuration = minutes( ...
    expData.nFrames / (frame_rate_hz * 60));
expData.signalEndTimes = expData.signalStartTimes + ...
    expData.signalDuration;

% ALERT: since we are ignoring the photobleaching window, we have to
%        shift the start of the signals by the window duration

% Get start/end time for each odor presentation (in minutes)
expData.odorStartTimes = expData.startTime + ...
    minutes(db_trials.odor_locs_min);
expData.odorDuration = minutes(odor_dur_s / 60);
expData.odorEndTimes = expData.odorStartTimes + expData.odorDuration;

% Get end time for the whole experiment
expData.endTime = expData.signalEndTimes(end);

end

function [dFF, zScore] = computeZScores(expData, baselineAvgs)
% COMPUTEZSCORE Computes dF/F and its z-score using the provided
% baseline averages for each ROI.

% Get sizes
nFrames = expData.nFrames;
nFiles = expData.nFiles;
nROIs = expData.nROIs;

% Compute and store dF/F
dFF(nFrames, nFiles, nROIs) = 0;

for iROI = 1:nROIs
    dFF(:, :, iROI) = ...
        (expData.signals(:, :, iROI) - baselineAvgs(iROI)) ...
        / baselineAvgs(iROI);
end

% Compute and store z-score
zScore(nFrames, nFiles, nROIs) = 0;

for iROI = 1:nROIs
    meanROI = mean(dFF(:, :, iROI), "all");
    stdROI = std(dFF(:, :, iROI), 0, "all");

    zScore(:, :, iROI) = (dFF(:, :, iROI) - meanROI) / stdROI;
end

end
