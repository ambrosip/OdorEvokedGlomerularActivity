%{ 
DOCUMENTATION
written by VA on Mar/2026

GOAL:
    Calculate the z-scores of the dF/F for each ROI across all
    acquisitions.

DEPENDS on:
    mat file created by timeSeriesFromFijiROIs
    (it only uses timing information + the struct s)

TO DO:
    - Highlight odor presentation intervals
    - Change x-axis to time (using h5 file)
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
yRange = [min(zScore(:)) max(zScore(:))];

xLimits = [0 nFiles*nFrames];

for iROI = 1:nROIs
    fig = figure('Name', strcat(figNameStart, '_ROI_', int2str(iROI)));
    
    signalROI = zScore(:,:,iROI);
    plot(signalROI(:));

    xlim(xLimits);
    ylim(yRange);

    drawnow;
end

