%% USER INPUT
% REMEMBER TO CLOSE GUI BEFORE ADVANCING PAST SAVE MASK SESSION

expDir = '/Volumes/MossLab/ImagingData/20260311/m357';

useRoisFromOtherExp = 0;
useRoisFromPython = 0;
matDirWithMaskInfo = "/Users/priscilla/Documents/Local - Moss Lab/20251007/sid260/e1/processed/matlab/2026-04-15/a20251007_sid260_e1_00001_mcor_to_a20251007_sid260_e1_00090_mcor_2026-04-15_masksFromGUI.mat";
should_I_plot_dFF = 0;
visibility = 'off';

plotRoiSubset = 1;
roiSubset = [1:20];

% if using mcor files from odyn (caiman python)
convertFrom32bit = true;

plotSubset = 0;
firstAcq = 1;
lastAcq = 2; 

saveFigs = 1;
saveWorkspace = 1;

manual_y_limits = 1;
ymax = 25;
ymin = -5;

manual_x_limits = 1;
xmin = -2;
xmax = 4;


%% Saves final mask from GUI

if useRoisFromOtherExp == 1
    load(matDirWithMaskInfo,'finalMask','nROIs','labelRGB')
elseif useRoisFromPython == 1
    load(matDirWithMaskInfo,'all_mask')
    nROIs = max(all_mask(:));
    finalMask = all_mask;
    labelRGB = label2rgb(finalMask, 'jet', 'k', 'shuffle'); 
else
    [finalMask, nROIs] = bwlabel(segmentationGUI.UserData.maskAfterExclusion > 0, 4);
    labelRGB = label2rgb(finalMask, 'jet', 'k', 'shuffle'); 
    segmentationGUI_UserData = segmentationGUI.UserData;
    % CLOSE GUI BEFORE CONTINUING
    close(segmentationGUI)
end


%% Update dirs

if exist('patchwarp','var')
    if patchwarp == 1
        getFileDirs_patchwarp
    else
        getFileDirs
    end
else
    getFileDirs
end

getImgDirs

% Create naming variables
todayStr = string(datetime('Today'), 'yyyy-MM-dd');
saveFolder = fullfile(expDir, 'processed', 'matlab', todayStr);
firstAcquisitionName = firstAcqName(1:end-9);
lastAcquisitionName = lastAcqName(1:end-9);

% create saveFolder if needed
% check if saveFolder exists
if not(isfolder(saveFolder))
    % create saveFolder
    mkdir(saveFolder);
end

% store mat file name
matFileName = strcat(firstAcquisitionName, '_to_', ...
    lastAcquisitionName, '_', todayStr, '_masksFromGUI');


%% Save workspace so far

close all
save(fullfile(saveFolder, matFileName));
disp('I saved the mat file');


%% Make z-proj of 1st mcor file for ROI labeling

filename = imgsToAnalyzeDirs(1).name;
fileDir = imgsToAnalyzeDirs(1).folder;
filepath = fullfile(fileDir, filename);
imgToProject = imread(filepath);
avgProjection = mean(imgToProject, 3);


%% Show ROIS

figName = strcat(firstAcqName(1:end), '_to_', lastAcqName(1:end), '_rois');
fig = figure('Name',figName);

if convertFrom32bit
    avgProjection = cast(avgProjection,'uint16');
    imshow(imadjust(avgProjection))
else
    imshow(imadjust(avgProjection,[0.5 0.65])) 
end
hold on;
finalMaskRGB = label2rgb(finalMask, 'jet', 'k', 'shuffle'); 
hOverlay = imshow(finalMaskRGB);
alphaMap = double(finalMask > 0) * 0.5; 
set(hOverlay, 'AlphaData', alphaMap);
hold off;

FigName = fig.Name;
set(0, 'CurrentFigure', fig);
% forces matlab to save fig as a vector
fig.Renderer = 'painters';  
% actually saves a vector file
saveas(fig,fullfile(saveFolder, [FigName '.pdf']));
close(fig);


%% Get fluorescence in fiji ROIs for each file/acquisition
% ADAPTED FROM TimeSeriesFromFijiROIsZscores (TODO: CLEAN UP)

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

    % Load whole movie and shift to be positive by adding min int16 value
    movieToAnalyze = single(tiffreadVolume(imgToAnalyzeFileDir));
    movieToAnalyze = movieToAnalyze + 32768;
    
    % Initialize a matrix of the correct size
    % The value 0.0 will be overwritten later
    meanPerROI(nFrames, nROIs) = 0.0;

    for iROI = 1:nROIs
        ROIMask = single(finalMask == iROI);
        nPixelInROI = sum(ROIMask, "all");
        
        meanPerROI(:,iROI) = sum(ROIMask .* movieToAnalyze, [1, 2]);
        meanPerROI(:,iROI) = meanPerROI(:, iROI) / nPixelInROI;
    end

    fileSignals.(f{1}) = meanPerROI;
    disp(strcat("Finished processing ROIs in file ", f))
end


%% Calculate Z-Scores of dF/F for each ROI across files

lastAcq = imgsToAnalyze_numberOf;

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
    baseline = dFF(1:adjustedBaselineEnd, :); % comment out if you want to calculate the zscore from F instead of dF
    file_dFF.(aFilenames{file}) = dFF;
    
    zdFF = (dFF - mean(baseline, 1)) ./ std(baseline, 0, 1); % comment if you want to calculate the zscore from F instead of dF
    % zdFF = (croppedSignal - meanBaseline) ./ std(baseline, 0, 1); % uncomment out if you want to calculate the zscore from F instead of dF
    fileZdFF.(aFilenames{file}) = zdFF;
end

% create x axis in seconds and adjust it based on the photobleaching window
croppedSignalFrames = size(croppedSignal, 1);
adjustedImageDuration = croppedSignalFrames / frame_rate_hz;
adjusted_baseline_dur_s = baseline_dur_s - photobleaching_window_s;

odorOnset = adjusted_baseline_dur_s;
odorOffset = adjusted_baseline_dur_s + odor_dur_s;
xs = linspace(-odorOnset, adjustedImageDuration - odorOnset, ...
    croppedSignalFrames);

disp("Calculated Z-scores of dF/F")


%% Save mat file so far

close all
save(fullfile(saveFolder, matFileName));
disp('I saved the mat file')


%% Prep for z plots

if manual_y_limits == 0
    % Gets the first integer below min and first above max
    ymax = ceil(max(structfun(@(x) max(x,[],'all'), ...
        fileZdFF,'UniformOutput',true)));
    ymin = floor(min(structfun(@(x) min(x,[],'all'), ...
        fileZdFF,'UniformOutput',true)));
end

if manual_x_limits == 0
    % Get the min and max value for x
    xmin = round(xs(1));
    xmax = round(xs(end));
end

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


%% Get average z-score per odor per block (i.e. program) per roi

% Get size from first zdFF matrix (frames per ROIs)
fileZdFF_fieldnames = fieldnames(fileZdFF);
[nFrames, nROIs] = size(fileZdFF.(fileZdFF_fieldnames{1}));

for programNum = unique(db_trials.programNum)'
    rowsToKeep_thisProgram = ismember(db_trials.programNum, programNum);

    for odorID = unique(db_trials.odorID)'
        rowsToKeep_thisOdor = ismember(db_trials.odorID, odorID);

        for outcome = unique(db_trials.outcome)'
            rowsToKeep_thisOutcome = ismember(db_trials.outcome, outcome);

            rowsToKeep =...
                logical(rowsToKeep_thisProgram .*...
                    rowsToKeep_thisOdor .*...
                    rowsToKeep_thisOutcome);

            acqsIdxToKeep = db_trials.acqIdx(rowsToKeep);
            acqsIdxToKeep = rmmissing(acqsIdxToKeep); % remove nans
            acqsIdxToKeep = acqsIdxToKeep(acqsIdxToKeep ~= 0); % remove zeros
            acqsTotal = numel(acqsIdxToKeep);

            fieldnamesToKeep = fileZdFF_fieldnames(acqsIdxToKeep);
            nAcq = length(fieldnamesToKeep);

            % Skips this combination of program, odor, outcome if there are
            % no acquisition with those parameters
            if nAcq == 0
                continue
            end
            
            % Preallocate the matrix by adding a 0 to the end.
            % It might have been easier to take the mean first, to not
            % store a big matrix, but this is easier to understand
            clear allFrames_allROIs_allAcq;
            allFrames_allROIs_allAcq(nFrames, nROIs, nAcq) = 0;

            % Concatenate all zdFF into a nFrames x nROIs x nAcq matrix
            for iAcquisition = 1:length(fieldnamesToKeep)
                fieldname = fieldnamesToKeep{iAcquisition};
                allFrames_allROIs_thisAcq = fileZdFF.(fieldname);
                allFrames_allROIs_allAcq(:, :, iAcquisition) = ...
                    allFrames_allROIs_thisAcq;
            end

            % Average over acquisitions (3rd dimension)
            allFrames_allROIs_meanAcq = ...
                mean(allFrames_allROIs_allAcq, 3, 'omitnan');
            
            % Fix field names for the structure
            programFieldName = strcat("program_", num2str(programNum));
            odorFieldName = strcat("odor_", num2str(odorID));

            if outcome == "false choice"
                outcome = "false_choice";
            end
    
            outcomeFieldName = outcome;

            % Add new matrix to the structure
            s_mean_zscore.(programFieldName).(odorFieldName).(outcomeFieldName) = ...
                allFrames_allROIs_meanAcq;

            % Also, storing the matrix before averaging
            s_zscore.(programFieldName).(odorFieldName).(outcomeFieldName) = ...
                allFrames_allROIs_allAcq;

            % TO DO: s_mean_zscore is redundant, remove it by changing
            % ------ the code that comes later (adding the relevant means)
        end
    end
end


%% PLOT AVGS - all outcomes (with SEM), color-coded by outcome (or odor, if outcome is nan)

allProgramFieldName = fieldnames(s_zscore);

% (number of odors) x (number of programs + extra columns)
% The extra columns are to put the plot showing the ROI
extraColumns = 3;

rows = length(allProgramFieldName);
columns = maxOdorNumber + extraColumns;

% xsx is the palindromic version of xs (needed for the fill plot)
xsx = [xs flip(xs)];

if plotRoiSubset == 1
    roiRange = roiSubset;
else
    roiRange = 1:nROI;
end

for iROI = roiRange
% for iROI = 1:nROIs
    % Create one figure per ROI
    figName = strcat( ...
        firstAcquisitionName, '_to_', lastAcquisitionName, ...
        '_roi_', num2str(iROI), '_zdFF_SEM');

    fig = figure('Name', figName, 'Visible',visibility);
    % fig = figure('Name', figName, 'Visible',"on");

    % Create tiledlayout of shape odors by programs
    t = tiledlayout(rows, columns);
    title(t, figName, 'Interpreter','none');

    % Iterate through programs and odors, getting the fieldnames as we go
    for iProgram = 1:length(allProgramFieldName)
        % Get the struct to shorten the names in the next for loops
        programFieldName = allProgramFieldName{iProgram};        
        programStruct = s_zscore.(programFieldName);
        allOdorFieldName = fieldnames(programStruct);

        % Get program type to make tile title
        programType = s_olfactometer.(programFieldName).type;
        odorList = s_olfactometer.(programFieldName).odorList;

        for iOdor = 1:length(allOdorFieldName)
            odorFieldName = allOdorFieldName{iOdor};
            odorStruct = programStruct.(odorFieldName);
            allOutcomeFieldName = fieldnames(odorStruct);
            
            % Create the plot
            nexttile(iOdor + (iProgram - 1) * columns)
            hold on;
            
            % Getting the odorID like Priscilla did
            odorID = extractBetween(odorList(iOdor), "I ", " -");

            % More plot settings + highlighting odor presentation window
            rectangle('Position',[0 ymin odor_dur_s ymax-ymin], ...
                'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
            yline(0,'k--')
            
            % Iterate through outcomes, plotting both on the same graph
            for iOutcome = 1:length(allOutcomeFieldName)
                outcomeFieldName = allOutcomeFieldName{iOutcome};

                if olfactory_task == "passive_odor_presentations"
                    % Get the odor color from user input in preProcessing_v2
                    color = odor_color.colorID(odor_color.odorID==str2double(odorID),:);
                else
                    % Get the outcome color from struct in user input
                    color = outcomeColors.(outcomeFieldName);
                end

                % Gets data for this ROI and all acquisitions
                allFrames_allROIs_allAcq = odorStruct.(outcomeFieldName);
                allFrames_thisROI_allAcq = ...
                    allFrames_allROIs_allAcq(:,iROI,:);

                % Get the mean and stderr (SEM) along acquisitions
                allFrames_thisROI_meanAcq = ...
                    mean(allFrames_thisROI_allAcq, 3);
                allFrames_thisROI_stderrAcq = ...
                    std(allFrames_thisROI_allAcq, 0, 3) / ...
                    sqrt(size(allFrames_thisROI_allAcq, 3));

                % Compute fill limits (xs were done at the start)
                upperBound = allFrames_thisROI_meanAcq + ...
                    allFrames_thisROI_stderrAcq;
                lowerBound = allFrames_thisROI_meanAcq - ...
                    allFrames_thisROI_stderrAcq;

                fill(xsx', [lowerBound; flip(upperBound)], color, ...
                    'FaceAlpha', 0.3, 'EdgeColor', 'none');

                % Lines are a little less transparent and thicker than
                % in other graphs (0.5 to 0.7 for both of them).
                plot(xs', allFrames_thisROI_meanAcq, ...
                    'Color', [color 1], 'LineWidth', 1);
            end

            hold off;
            disp(strcat("plot odor ", odorID, " done"))

            % Plot aesthetics            
            axis([xmin xmax ymin ymax])
            if olfactory_task == "passive_odor_presentations"
                set(gcf, 'OuterPosition', [100 100 1600 400]);
            else 
                set(gcf, 'OuterPosition', [100 100 1600 900]);
            end
            set(gca, 'LineWidth', 0.75);
            set(gca, 'FontName', 'Arial');
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
            title(strcat(odorFieldName, '_', programType), ...
                'Interpreter', 'none');
            xlabel('Time from odor onset (s)')
            ylabel('dF/F Z-score')
            xticks([xmin,0,1,xmax]);
            yticks([ymin,0,ymax]);           
        end
    end

    % Last tile showing the ROI
    nexttile(columns - extraColumns + 1, [rows, extraColumns])
    
    hold on;
    if convertFrom32bit
        avgProjection = cast(avgProjection,'uint16');
        imshow(imadjust(avgProjection));
    else
        imshow(imadjust(avgProjection,[0.5 0.65]));
    end
    
    hOverlay = imshow(labelRGB);
    alphaMap = double(finalMask == iROI) * 0.5; 
    set(hOverlay, 'AlphaData', alphaMap);

    hold off;  

    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    disp(strcat("plot roi ", num2str(iROI), " done"))

    if saveFigs
        FigName = fig.Name;
        set(0, 'CurrentFigure', fig);
        % forces matlab to save fig as a vector
        fig.Renderer = 'painters';  
        % actually saves a vector file
        saveas(fig,fullfile(saveFolder, [FigName '.svg']));
        close(fig);
    end

end


%% Plot dFF (if user wants it)

if should_I_plot_dFF == 1
    plot_dFF
end


%% Save figs

% if saveFigs
%     FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% 
%     % save all open figs
%     for iFig = 1:length(FigList)
%       FigHandle = FigList(iFig);
%       FigName = FigList(iFig).Name;
%       set(0, 'CurrentFigure', FigHandle);
%       % forces matlab to save fig as a vector
%       FigHandle.Renderer = 'painters';  
%       % actually saves a vector file
%       saveas(FigHandle,fullfile(saveFolder, [FigName '.svg']));
%     end 
%     disp('saved all figs')
%     close all
% end


%% Save workspace

close all
save(fullfile(saveFolder, matFileName));
disp('I saved the mat file');



