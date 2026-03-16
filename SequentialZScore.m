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

%% Compute Z-score

% Get first file signal
files = string(fieldnames(s));
firstFile = files(1);
firstSignal = s.(firstFile);

% Get figure name (start of file name)
figNameStart = split(string(firstAcqName), '_');
figNameStart = join(figNameStart(1:end-2), '_');

% Get baseline
baselineStart = ceil(photobleaching_window_frames);
baselineEnd = floor(odor_onset_s * frame_rate_hz);
baselineAvgs = mean(firstSignal(baselineStart:baselineEnd, :), 1);

% Initialize dFF
nFiles = length(files);
nFrames = size(firstSignal(baselineStart:end, :), 1);
nROIs = size(firstSignal, 2);

% ALERT: We remove the photobleaching window in all acquisitions, this is
%        why the nFrames above is not just size(firstSignal, 1).

dFF(nFrames, nFiles, nROIs) = 0;

% Compute dFF for all files
for iFile = 1:nFiles
    signal = s.(files(iFile));

    % Crops signal to exclude photobleaching window
    signal = signal(baselineStart:end, :);
    
    for iROI = 1:nROIs
        dFF(:, iFile, iROI) = ...
            (signal(:, iROI) - baselineAvgs(iROI)) / baselineAvgs(iROI);
    end
end

% Compute Z-score
zScore(nFrames, nFiles, nROIs) = 0;

for iROI = 1:nROIs
    meanROI = mean(dFF(:, :, iROI), "all");
    stdROI = std(dFF(:, :, iROI), 0, "all");
    
    zScore(:, :, iROI) = (dFF(:, :, iROI) - meanROI) / stdROI;
end

%% Plot Z-Score

% Find y range by getting max and min across ROIs
zMin = min(zScore(:));
zMax = max(zScore(:));

% Get start times for each acquisition (in minutes)
% ALERT: since we are ignoring the photobleaching window, we have to shift
%        the start of the plots by the window duration

acqStartTimes = db_trials.trial_locs_min + photobleaching_window_s / 60;
acqDuration = nFrames / (frame_rate_hz * 60);

% Get start time for each odor presentation
odorStartTimes = db_trials.odor_locs_min;

% Compute the time range for the whole experiment
expEnd = acqStartTimes(end) + acqDuration;
expTimeRange = [0 expEnd];

for iROI = 1:nROIs
    figName = strcat(figNameStart, '_ROI_', int2str(iROI));

    fig = figure('Name', figName);
    title(figName, 'Interpreter', 'none');
    
    % Plot acquisitions with different colors
    hold on;
    for iFiles = 1:nFiles
        % Compute x range for each file (in minutes)
        acqStart = acqStartTimes(iFiles);
        acqEnd = acqStart + acqDuration;
        
        timeRange = linspace(acqStart, acqEnd, nFrames);
        timeRange = minutes(timeRange);

        % Highlight odor presentation
        odorStart = minutes(odorStartTimes(iFiles));
        odorEnd = odorStart + minutes(odor_dur_s / 60);

        patch( ...
            [odorStart odorEnd odorEnd odorStart], ...
            [zMin zMin zMax zMax], ...
            [0.9 0.9 0.9], ...
            'EdgeColor', 'none');
        
        % Plot files
        signalROI = zScore(:, iFiles, iROI);
        plot(timeRange, signalROI(:));

    end
    hold off;

    % Set plot ranges
    xlim(minutes(expTimeRange));
    ylim([zMin zMax]);

    % Format the x-axis to show minutes and seconds
    xtickformat('mm:ss');
    xlabel('Time (min:sec)')

    ylabel('Z-Score for dF/F')

    % Distribute figures across the screen
    % ASSUMPTION: There are less than 16 ROIs (it still works with more,
    %             but they will be clumped in the right side of the screen)
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

