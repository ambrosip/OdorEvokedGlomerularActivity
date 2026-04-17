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
% USE DOUBLE QUOTES, NOT SINGLE QUOTES
% PUT ",..." BETWEEN EXP
expFolders = [ ...
    "M:\ImagingData\20251002\M174\e1"...
    ];

plotRoiSubset = 1;
roiRange = [1:20];

% File extension for the plot images
fileExtension = ".svg";

% Ranges of y-axis for z-scores or dF/F plots
% Leave it as [Inf -Inf] for automatic range
dFFRange = [Inf -Inf];
zScoreRange = [-5 10];

% Save figures in saveDir or not?
% If false, it will save figure in:
%   EXPFOLDER/processed/matlab/SCRIPT_RUN_DATE
useSaveDir = false;


%% Load relevant .mat files

nFolders = length(expFolders);

% Create an struct array of experiment data/metadata
% Since there are only a few folders, we do it directly inside the loop

for iFolder = 1:nFolders
    % ALERT: loadData is defined at the end of the file.
    [expData(iFolder), saveDir] = loadData(expFolders(iFolder), iFolder);
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

% ALERT: 1) saves a mat file with data used in the figure
%        2) plotFigures is defined in the next section

if useSaveDir
    saveFolder = saveDir;
else
    saveFolder = '';
end

% Just some padding
fprintf('\n');

plotFigures(expData, fileExtension, 'dFF', dFFRange, saveFolder);
plotFigures(expData, fileExtension, 'zScore', zScoreRange, saveFolder);


%%  Auxiliary Functions

function [expData, saveDir] = loadData(expFolder, i)
% LOADDATA Load and compute all data needed to compute the z-score

% Find the outputs of the timeSeriesFromFijiROIs script
% matPaths = dir(fullfile( ...
%     expFolder, '**', '*timeSeriesFromFijiROIs.mat'));

% matPaths = dir(fullfile( ...
%     expFolder, '**', '*masksFromZScores.mat'));

matPaths = dir(fullfile( ...
    expFolder, '**', '*masksFromGUI.mat'));

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
% load(matPath, ...
%     's', 'db_trials', 'firstAcqName', 'currentImgInfo', ...
%     'photobleaching_window_frames', 'photobleaching_window_s', ...
%     'frame_rate_hz', 'odor_onset_s', 'odor_dur_s', 'saveDir'); % ALERT ALERT ALERT
load(matPath, ...
    'fileSignals', 'db_trials', 'firstAcqName', 'currentImgInfo', ...
    'photobleaching_window_s', ...
    'frame_rate_hz', 'odor_onset_s', 'odor_dur_s', 'saveDir'); % ALERT ALERT ALERT

% Create output struct

% Get experiment name (start of TIFF file names)
nameStart = split(string(firstAcqName), '_');
expData.name = join(nameStart(1:end-2), '_');

% Store folder input
expData.folder = expFolder;

% Store db_trials
expData.db_trials = db_trials;

% Get first file signal
% files = string(fieldnames(s)); % ALERT ALERT ALERT
files = string(fieldnames(fileSignals)); % ALERT ALERT ALERT
firstFile = files(1);
% firstSignal = s.(firstFile); % ALERT ALERT ALERT
firstSignal = fileSignals.(firstFile); % ALERT ALERT ALERT

% Get important "timestamps" (actually data points)
baselineStart = ceil(photobleaching_window_s * frame_rate_hz);
baselineEnd = floor(odor_onset_s * frame_rate_hz);
odorStart = ceil(odor_onset_s * frame_rate_hz);
odorEnd = ceil((odor_onset_s + odor_dur_s) * frame_rate_hz);

% get baseline
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
    % signal = s.(files(iFile)); % ALERT ALERT ALERT
    signal = fileSignals.(files(iFile)); % ALERT ALERT ALERT

    % Crops signal to exclude photobleaching window
    signal = signal(baselineStart:end, :);

    % Group all signals together for easy processing
    expData.signals(:, iFile, :) = signal;
end

% correct important "timestamps" in datapoints after cropping the
% photobleaching window
expData.baselineStart = baselineStart - baselineStart;
expData.baselineEnd = baselineEnd - baselineStart;
expData.odorStart = odorStart - baselineStart;
expData.odorEnd = odorEnd - baselineStart;

% Get timing data

% Get start time for the whole experiment
startTimeStr = extractBetween( ...
    currentImgInfo(1).ImageDescription,"epoch = [","]");
expData.startTime = ...
    datetime(startTimeStr, 'InputFormat', 'yyyy M d H m s.SSS');
    % datetime(startTimeStr, 'InputFormat', 'yyyy M d H m s.SSS') + minutes(i*40); % ALERT ALERT ALERT

% Get start/end times for each acquisition (in minutes)
expData.signalStartTimes = expData.startTime + minutes( ...
    db_trials.trial_locs_min + photobleaching_window_s / 60);
expData.signalDuration = minutes( ...
    expData.nFrames / (frame_rate_hz * 60));
expData.signalEndTimes = expData.signalStartTimes + ...
    expData.signalDuration;

% crop these in case you have missing acquisitions at the end, which
% happens in file 20260316_m357_e1
if length(expData.signalEndTimes) > expData.nFiles
    expData.signalEndTimes = expData.signalEndTimes(1:expData.nFiles);
    expData.signalStartTimes = expData.signalStartTimes(1:expData.nFiles);
end

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

function plotFigures(expData, fileExtension, ...
                                  plotType, yRange, saveFolder)
% PLOTFIGURES Plots and saves dF/F or z-scores to the saveFolder
% (default is EXPFOLDER/processed/matlab/DATE).
%
% Input Arguments
%   plotType - Choose whether to plot dF/F or z-scores
%       string 'dFF' | string 'zScore'

arguments
    expData
    fileExtension {mustBeTextScalar}
    plotType (1,1) string {mustBeMember(plotType, ["dFF", "zScore"])}
    yRange (2,1) double = [Inf -Inf]
    saveFolder {mustBeTextScalar} = ""
end

% Get plots save folder (if not supplied)
if strlength(saveFolder) == 0
    todayStr = string(datetime('Today'), 'yyyy-MM-dd');
    saveFolder = fullfile(expData(1).folder, ...
        'processed', 'matlab', todayStr);
end

fprintf(' Plotting %s figures\n', plotType);
fprintf('   Saving figures as %s\n', fileExtension);
fprintf('   Save folder:\n')
fprintf('     %s\n', saveFolder);

if ~isfolder(saveFolder)
    mkdir(saveFolder);
end

% Find range of y values, if not supplied
% Otherwise use the supplied values
if yRange(1) > yRange(2)
    yMin = Inf;
    yMax = -Inf;

    for data = expData
        yMin = min(yMin, min(data.(plotType)(:)));
        yMax = max(yMax, max(data.(plotType)(:)));
    end
else
    yMin = yRange(1);
    yMax = yRange(2);
end

% Compute the time range for the plot

% We shift the time range for the first experiment to start at 0.
% ASSUMPTION: experiments are in chronological order

startTime = expData(1).startTime;
endTime = expData(end).endTime;
plotTimeRange = [0 (endTime - startTime)];

% ASSUMPTION: All experiments have the same number of ROIs
nROIs = expData(1).nROIs;

if plotRoiSubset == 1
    roiRange = roiSubset;
else
    roiRange = 1:nROIs;
end

for iROI = roiRange
% for iROI = 1:nROIs
    % Name figure
    figName = strcat(expData(1).name, ...
        '_to_', expData(end).name, ...
        '_ROI_', int2str(iROI), ...
        '_SequentialZScore_', plotType);

    fig = figure('Name', figName, 'Visible', 'off');
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
                [yMin yMin yMax yMax], ...
                [0.9 0.9 0.9], ...
                'EdgeColor', 'none');

            % Plot files
            yROI = data.(plotType)(:, iFiles, iROI);
            plot(timeRange, yROI(:));

        end
    end

    hold off;

    % Set plot ranges
    xlim(plotTimeRange);
    ylim([yMin yMax]);

    % Format the x-axis to show minutes and seconds
    xtickformat('mm:ss');
    xlabel('Time (min:sec)')

    % Put the appropriate y-label
    if strcmp(plotType, 'zScore')
        ylabel('Z-Score for dF/F')
    elseif strcmp(plotType, 'dFF')
        ylabel('dF/F')
    end

    % Distribute figures across the screen

    % ASSUMPTION: There are less than 16 ROIs (it still works with more,
    %             but they will be clumped in the right side of the screen)

    % ALERT: ROIs are arranged from bottom to top and from left to right.
    figWidth = 1350;
    figHeight = 150;
    vPadding = 80;
    hPadding = 10;

    fig.Position = [ ...
        10 ...
        (20 + mod(iROI - 1, 4) * (figHeight + vPadding)) ...
        figWidth ...
        figHeight];

    saveas(fig, fullfile(saveFolder, strcat(figName, fileExtension)));
    close(fig);
end

% save workspace variables
matFileName = strcat(expData(1).name, "_to_", ...
    expData(end).name,'_sequentialZScore_', plotType);

fprintf(['   Saving mat file with plotFigure %s data' ...
                            ' (not the whole workspace).\n'], plotType);
fprintf(['   .mat file name:\n' ...
         '     %s.mat\n'], matFileName);
fprintf( '\n');

save(fullfile(saveFolder, matFileName));

end


%%

% save workspace variables
matFileName = strcat(expData(1).name, "_to_", ...
    expData(end).name,'_sequentialZScore');

save(fullfile(saveDir, matFileName));