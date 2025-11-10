%% Get dir of fiji files

% Get fiji ROI file (roi = regions of interest)
% ALERT: only analyzing the 1st .zip file found
ROIFile = dir(fullfile(expDir,'processed','fiji','*.zip'));
ROIFile = remove_stupid_mac_hidden_files(ROIFile);
ROIs = ReadImageJROI(fullfile(ROIFile.folder, ROIFile.name));

% Get fiji Avg Intensity Projection file (zProj = z Projection)
zProjFiles = dir(fullfile(expDir, 'processed', 'fiji', '*.tif'));
zProjFiles = remove_stupid_mac_hidden_files(zProjFiles);

% ALERT: only using the 1st .tif file found
zProjFileToAnalyzeDir = fullfile(zProjFiles(1).folder, zProjFiles(1).name);
zProjFileToAnalyze = imread(zProjFileToAnalyzeDir);

% vnImageSize is one of the inputs for ROIs2Regions
% for whatever reason, it is a translated version of the image size
vnImageSize = flip(size(zProjFileToAnalyze));
regions = ROIs2Regions(ROIs,vnImageSize);
nROIs = length(ROIs);

disp('Loaded Fiji data')

%% Get fluorescence in fiji ROIs for each file/acquisition

meanIntPerRoi = [];
for file = 1:imgsToAnalyze_numberOf
    
    % get OS-appropriate file dir
    imgToAnalyzeFileDir = fullfile(imgsToAnalyzeFolder, ...
        imgsToAnalyzeDirs(file).name);

    % get img file name without extension (stored in "f")
    [p,f,e] = fileparts(imgsToAnalyzeNames(file));

    % add "a" in front of the filename in f to build a structure later
    % why? matlab freaks out if field names of a structure start with numbers
    f = {strcat('a', cell2mat(f))};

    % get img info
    imgInfo = imfinfo(imgToAnalyzeFileDir);
    nFrames = length(imgInfo);

    % Load whole movie and shift to be positive
    movieToAnalyze = single(tiffreadVolume(imgToAnalyzeFileDir));
    movieToAnalyze = movieToAnalyze - min(movieToAnalyze, [], 'all');

    % Initialize a matrix of the correct size
    % The value 0.0 will be overwritten later
    meanPerROI(nFrames, nROIs) = 0.0;

    for iROI = 1:nROIs
        ROIMask = single(labelmatrix(regions) == iROI);
        ROIMask = ROIMask';
        nPixelInROI = sum(ROIMask, "all");
        
        meanPerROI(:,iROI) = sum(ROIMask .* movieToAnalyze, [1, 2]);
        meanPerROI(:,iROI) = meanPerROI(:, iROI) / nPixelInROI;
    end


    fileSignals.(f{1}) = meanPerROI;
    disp(strcat("Finished processing ROIs in file ", f))
end

%% Calculate Z-Scores of dF/F for each ROI across files

if plotSubset
    firstAcq = 1;
    lastAcq = length(imgToAnalyzeFileDir);
end

aFilenames = fieldnames(fileSignals);
firstAcquisitionName = aFilenames{firstAcq};
lastAcquisitionName = aFilenames{lastAcq};

photobleachingWindowEnd = round(photobleaching_window_s * frame_rate_hz);
adjustedBaselineEnd = round( ...
    (baseline_dur_s - photobleaching_window_s) * frame_rate_hz);


for file = firstAcq:lastAcq
    croppedSignal = fileSignals.(aFilenames{file});
    croppedSignal = croppedSignal(photobleachingWindowEnd:end, :);

    baseline = croppedSignal(1:adjustedBaselineEnd, :);
    meanBaseline = mean(baseline, 1, 'omitnan');
    
    dFF = (croppedSignal - meanBaseline) ./ meanBaseline;
    baseline = dFF(1:adjustedBaselineEnd, :);
    
    zdFF = (dFF - mean(dFF, 1)) ./ std(dFF, 0, 1);
    fileZdFF.(aFilenames{file}) = zdFF;
end

% create x axis in seconds and adjust it based on the photobleaching window
croppedSignalFrames = size(croppedSignal, 1);
adjustedImageDuration = croppedSignalFrames / frame_rate_hz;
adjusted_baseline_dur_s = baseline_dur_s - photobleaching_window_s;

odorOnset = adjusted_baseline_dur_s;
odorOffset = adjusted_baseline_dur_s + odor_dur_s;
xs = linspace(-odorOnset, adjustedImageDuration - odorOffset, ...
    croppedSignalFrames);

disp("Calculated Z-scores of dF/F")

%% PLOT data in ROIs organized by odor (rows) and program type (columns)

% Gets the first integer below min and first above max
ymax = ceil(max(structfun(@(x) max(x,[],'all'), ...
    fileZdFF,'UniformOutput',true)));
ymin = floor(min(structfun(@(x) min(x,[],'all'), ...
    fileZdFF,'UniformOutput',true)));

% Get the min and max value for x
xmin = round(xs(1));
xmax = round(xs(end));

% Get the max number of odors used in this experiment
maxOdorNumber = 0;
for iProgram = 1:length(programFieldNames)
    programFieldName = programFieldNames(iProgram);

    if s_olfactometer.(programFieldName).type ~= "ignore"
        odorNumber = length(s_olfactometer.(programFieldName).odorList);
        if maxOdorNumber < odorNumber
            maxOdorNumber = odorNumber;
        end
    end
end

% Number of rows and columns on the tiledlayout
rows = maxOdorNumber;
columns = programsToAnalyze + 1;

for iROI = 1:nROIs
    figName = strcat( ...
        firstAcquisitionName, '_to_', lastAcquisitionName, ...
        '_roi_', num2str(iROI), '_zdFF');

    fig = figure('Name', figName);

    % Plot settings
    set(gca, 'FontName', 'Arial');
    set(gcf, 'OuterPosition', [100 100 500 900]);
    set(gca, 'LineWidth', 0.75);

    t = tiledlayout(rows, columns);
    title(t, figName);

    iProgramRelative = 0;
    for iProgram = 1:length(programFieldNames)
        programFieldName = programFieldNames(iProgram);
        
        if s_olfactometer.(programFieldName).type == "ignore"
            continue
        end

        iProgramRelative = iProgramRelative + 1;
        odorList = s_olfactometer.(programFieldName).odorList;
        programData = s_olfactometer.(programFieldName);
        
        for iOdor = 1:length(odorList)
            nexttile(iProgramRelative + (iOdor-1)*columns)
            odorID = extractBetween(odorList(iOdor), "I", " -");
            
            color = odor_color.colorID( ...
                odor_color.odorID == str2double(odorID), :);

            hold on;
            rectangle('Position',[0 ymin odor_dur_s ymax-ymin], ...
                'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
            yline(0,'k--')
            
            axis([xmin xmax ymin ymax])
            title(strcat(odorFieldName, '_', programData.type));
            xlabel('Time from odor onset (s)')
            ylabel('dF/F Z-score')

            idx = programData.summary_by_trial.odor==str2double(odorID);
            for acqIdx = programData.summary_by_trial.acqIdx(idx)'
                % Plot Z-score
                if isnan(acqIdx) | acqIdx > lastAcq
                    continue
                end

                plot(xs',fileZdFF.(aFilenames{acqIdx})(:,iROI), ...
                    'Color',[color 0.5],'LineWidth',0.5);
            end

            hold off;
            disp(strcat("plot odor ", odorID, " done"))
        end
    end

    nexttile(columns, [maxOdorNumber, 1])
    imshow(imadjust(zProjFileToAnalyze,[0.5 0.65]))
    
    hold on
    thetas = linspace(0,2*pi,200);
    ellipseR1 = (ROIs{iROI}.vnRectBounds(4) - ROIs{iROI}.vnRectBounds(2))/2;
    ellipseR2 = (ROIs{iROI}.vnRectBounds(3) - ROIs{iROI}.vnRectBounds(1))/2;
    ellipseA = (ROIs{iROI}.vnRectBounds(4) + ROIs{iROI}.vnRectBounds(2))/2;
    ellipseB = (ROIs{iROI}.vnRectBounds(3) + ROIs{iROI}.vnRectBounds(1))/2;
    ellipseX = ellipseR1*cos(thetas)+ellipseA;
    ellipseY = ellipseR2*sin(thetas)+ellipseB; 
    plot(ellipseX,ellipseY,'Color','y','LineWidth',1);
    hold off  
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    disp(strcat("plot roi ", num2str(iROI), " done"))
end