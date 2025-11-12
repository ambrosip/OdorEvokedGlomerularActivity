%% USER INPUT

plotROIsubset = 0;
ROIsubset = [1, 2, 4, 7, 12, 14, 18, 19, 21, 24, 25, 26, 33, 34, 35];

% plotSubset = 0;
saveFigs = 1;
saveWorkspace = 1;

% colorblind-safe colors (source:
% https://www.nature.com/articles/nmeth.1618)
black_color = [0 0 0]/255;
orange_color = [230 159 0]/255;
sky_blue_color = [86 180 233]/255;
bluish_green_color = [0 158 115]/255;
yellow_color = [240 228 66]/255;
blue_color = [0 114 178]/255;
vermillion_color = [213 94 0]/255;
reddish_purple_color = [204 121 167]/255;

% colors for poster - source:
% https://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3
soft_green_color = [102 194 165]/255;
soft_orange_color = [252 141 98]/255;
soft_purple_color = [141 160 203]/255;

% colors for poster - source:
% https://colorbrewer2.org/#type=diverging&scheme=BrBG&n=4
dark_brown_color = [166 97 26]/255;
light_brown_color = [223 194 125]/255;
light_teal_color = [128 205 193]/255;
dark_teal_color = [1 133 113]/255;

% odor color-coding
odor_ids = [1; 2;...
    17; 18; 19; 20;...
    21; 22; 23; 24;...
    25; 26; 27; 28];
color_ids = [black_color; reddish_purple_color;...
    dark_brown_color; light_brown_color; light_teal_color; dark_teal_color;...
    blue_color; sky_blue_color; vermillion_color; orange_color;...
    blue_color; sky_blue_color; vermillion_color; orange_color];
odor_color = table(odor_ids,color_ids,'VariableNames',{'odorID','colorID'});


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

% % set default firstFig and lastFig boundaries in case user does NOT want a
% % custom subset
% if plotSubset == 0
%     firstAcq = 1;
%     lastAcq = imgsToAnalyze_numberOf;
% end

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
xs = linspace(-odorOnset, adjustedImageDuration - odorOnset, ...
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
    set(gcf, 'OuterPosition', [100 100 1000 900]);
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
            title(strcat(odorFieldName, '_', programData.type),'Interpreter','none');
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
    imshow(imadjust(cast(zProjFileToAnalyze,'uint16')))
    
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


%% Get average z-score per odor per block (program) per roi

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
                boolean(rowsToKeep_thisProgram .*...
                    rowsToKeep_thisOdor .*...
                    rowsToKeep_thisOutcome);

            acqsIdxToKeep = db_trials.acqIdx(rowsToKeep);
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


%% PLOT AVGS - hits and false choices, color-coded by outcome

allProgramFieldName = fieldnames(s_mean_zscore);

% (number of odors) x (number of programs + 1)
% The extra column is to put the plot showing the ROI
rows = maxOdorNumber;
columns = length(allProgramFieldName) + 1;

% Struct with the outcome colors, for easy reference
outcomeColors = struct( ...
    'hit', [56 83 163]/255, ...
    'false_choice', [236 30 36]/255, ...
    'miss',  [127 127 127]/255 ...
    );

for iROI = 1:nROIs
    % Create one figure per ROI
    figName = strcat( ...
        firstAcquisitionName, '_to_', lastAcquisitionName, ...
        '_roi_', num2str(iROI), '_zdFF_AVG');

    fig = figure('Name', figName);

    % Plot settings
    set(gca, 'FontName', 'Arial');
    set(gcf, 'OuterPosition', [100 100 1000 900]);
    set(gca, 'LineWidth', 0.75);

    % Create tiledlayout of shape odors by programs
    t = tiledlayout(rows, columns);
    title(t, figName, 'Interpreter','none');

    % Iterate through programs and odors, getting the fieldnames as we go
    for iProgram = 1:length(allProgramFieldName)
        % Get the struct to shorten the names in the next for loops
        programFieldName = allProgramFieldName{iProgram};        
        programStruct = s_mean_zscore.(programFieldName);
        allOdorFieldName = fieldnames(programStruct);

        % Get program type to make tile title
        programType = s_olfactometer.(programFieldName).type;
        odorList = s_olfactometer.(programFieldName).odorList;

        for iOdor = 1:length(allOdorFieldName)
            odorFieldName = allOdorFieldName{iOdor};
            odorStruct = programStruct.(odorFieldName);
            allOutcomeFieldName = fieldnames(odorStruct);
            
            % Create the plot
            nexttile(iProgram + (iOdor - 1) * columns)
            hold on;
            
            % Getting the odorID like Priscilla did
            odorID = extractBetween(odorList(iOdor), "I ", " -");

            % More plot settings + highlighting odor presentation window
            rectangle('Position',[0 ymin odor_dur_s ymax-ymin], ...
                'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
            yline(0,'k--')

            axis([xmin xmax ymin ymax])
            title(strcat(odorFieldName, '_', programType), ...
                'Interpreter', 'none');
            xlabel('Time from odor onset (s)')
            ylabel('dF/F Z-score')
            
            % Iterate through outcomes, plotting both on the same graph
            for iOutcome = 1:length(allOutcomeFieldName)
                outcomeFieldName = allOutcomeFieldName{iOutcome};

                % This has all the ROIs so we get a slice below
                allFrames_allROIs = odorStruct.(outcomeFieldName);

                % Get the color from struct
                color = outcomeColors.(outcomeFieldName);

                % Lines are a little less transparent and thicker than
                % in other graphs (0.5 to 0.7 for both of them).
                plot(xs', allFrames_allROIs(:, iROI), ...
                    'Color', [color 0.7], 'LineWidth', 0.7);
            end

            hold off;
            disp(strcat("plot odor ", odorID, " done"))
        end
    end

    % Last tile showing the ROI
    nexttile(columns, [maxOdorNumber, 1])
    
    hold on
    imshow(imadjust(cast(zProjFileToAnalyze,'uint16')))

    thetas = linspace(0,2*pi,200);
    ellipseR1 = (ROIs{iROI}.vnRectBounds(4) - ROIs{iROI}.vnRectBounds(2))/2;
    ellipseR2 = (ROIs{iROI}.vnRectBounds(3) - ROIs{iROI}.vnRectBounds(1))/2;
    ellipseA = (ROIs{iROI}.vnRectBounds(4) + ROIs{iROI}.vnRectBounds(2))/2;
    ellipseB = (ROIs{iROI}.vnRectBounds(3) + ROIs{iROI}.vnRectBounds(1))/2;
    ellipseX = ellipseR1 * cos(thetas) + ellipseA;
    ellipseY = ellipseR2 * sin(thetas) + ellipseB; 
    plot(ellipseX, ellipseY, 'Color', 'y', 'LineWidth', 1);

    hold off  

    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    disp(strcat("plot roi ", num2str(iROI), " done"))
end


%% PLOT AVGS - hits only, color-coded by odor

allProgramFieldName = fieldnames(s_mean_zscore);

% (number of programs) x (columns spanned by ROI plot + 1)
% One column for all odors and the  rest for ROI plot
rows = length(allProgramFieldName);
columns = 4;

if plotROIsubset == 1   
    for iROI = ROIsubset
        % Create one figure per ROI
        figName = strcat( ...
            firstAcquisitionName, '_to_', lastAcquisitionName, ...
            '_roi_', num2str(iROI), '_zdFF_AVG_only_hits');
    
        fig = figure('Name', figName);
    
        % Create tiledlayout of shape odors by programs
        t = tiledlayout(rows, columns);
        title(t, figName, 'Interpreter', 'none');
    
        % Iterate through programs and odors, getting the fieldnames as we go
        for iProgram = 1:length(allProgramFieldName)
            % Get the struct to shorten the names in the next for loops
            programFieldName = allProgramFieldName{iProgram};        
            programStruct = s_mean_zscore.(programFieldName);
            allOdorFieldName = fieldnames(programStruct);
    
            % Get program type to make tile title
            programType = s_olfactometer.(programFieldName).type;
            odorList = s_olfactometer.(programFieldName).odorList;
    
            % Create the plot
            nexttile(1 + (iProgram - 1) * columns)
            hold on;
    
            % More plot settings + highlighting odor presentation window
            rectangle('Position',[0 ymin odor_dur_s ymax-ymin], ...
                'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
            yline(0,'k--')
    
            axis([xmin xmax ymin ymax])
            title(strcat(odorFieldName, '_', programType), ...
                'Interpreter', 'none');
            xlabel('Time from odor onset (s)')
            ylabel('dF/F Z-score')
    
            for iOdor = 1:length(allOdorFieldName)
                odorFieldName = allOdorFieldName{iOdor};
                odorStruct = programStruct.(odorFieldName);
                allOutcomeFieldName = fieldnames(odorStruct);
                
                % Getting the odorID like Priscilla did
                odorID = extractBetween(odorList(iOdor), "I ", " -");
    
                color = color_ids(find(odor_ids == str2double(odorID), 1), :);
    
                % This has all the ROIs so we get a slice below
                allFrames_allROIs = odorStruct.hit;
    
                % Lines are a little less transparent and thicker than
                % in other graphs (0.5 to 0.7 for both of them).
                plot(xs', allFrames_allROIs(:, iROI), ...
                        'Color', [color 0.7], 'LineWidth', 2);
    
                disp(strcat("plot odor ", odorID, " done"))
            end
    
            hold off;
    
            % Plot settings
            set(gca, 'FontName', 'Arial');
            set(gcf, 'OuterPosition', [100 100 1000 900]);
            set(gca, 'LineWidth', 0.75);
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
    
        end
    
        % Last tile showing the ROI
        nexttile(2, [rows, columns - 1])
        
        hold on
        imshow(imadjust(cast(zProjFileToAnalyze,'uint16')))
    
        thetas = linspace(0,2*pi,200);
        ellipseR1 = (ROIs{iROI}.vnRectBounds(4) - ROIs{iROI}.vnRectBounds(2))/2;
        ellipseR2 = (ROIs{iROI}.vnRectBounds(3) - ROIs{iROI}.vnRectBounds(1))/2;
        ellipseA = (ROIs{iROI}.vnRectBounds(4) + ROIs{iROI}.vnRectBounds(2))/2;
        ellipseB = (ROIs{iROI}.vnRectBounds(3) + ROIs{iROI}.vnRectBounds(1))/2;
        ellipseX = ellipseR1 * cos(thetas) + ellipseA;
        ellipseY = ellipseR2 * sin(thetas) + ellipseB; 
        plot(ellipseX, ellipseY, 'Color', 'y', 'LineWidth', 1);
    
        hold off  
    
        t.TileSpacing = 'compact';
        t.Padding = 'compact';
        disp(strcat("plot roi ", num2str(iROI), " done"))
    end
else
    for iROI = 1:nROIs
        % Create one figure per ROI
        figName = strcat( ...
            firstAcquisitionName, '_to_', lastAcquisitionName, ...
            '_roi_', num2str(iROI), '_zdFF_AVG_only_hits');
    
        fig = figure('Name', figName);
    
        % Create tiledlayout of shape odors by programs
        t = tiledlayout(rows, columns);
        title(t, figName, 'Interpreter', 'none');
    
        % Iterate through programs and odors, getting the fieldnames as we go
        for iProgram = 1:length(allProgramFieldName)
            % Get the struct to shorten the names in the next for loops
            programFieldName = allProgramFieldName{iProgram};        
            programStruct = s_mean_zscore.(programFieldName);
            allOdorFieldName = fieldnames(programStruct);
    
            % Get program type to make tile title
            programType = s_olfactometer.(programFieldName).type;
            odorList = s_olfactometer.(programFieldName).odorList;
    
            % Create the plot
            nexttile(1 + (iProgram - 1) * columns)
            hold on;
    
            % More plot settings + highlighting odor presentation window
            rectangle('Position',[0 ymin odor_dur_s ymax-ymin], ...
                'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
            yline(0,'k--')
    
            axis([xmin xmax ymin ymax])
            title(strcat(odorFieldName, '_', programType), ...
                'Interpreter', 'none');
            xlabel('Time from odor onset (s)')
            ylabel('dF/F Z-score')
    
            for iOdor = 1:length(allOdorFieldName)
                odorFieldName = allOdorFieldName{iOdor};
                odorStruct = programStruct.(odorFieldName);
                allOutcomeFieldName = fieldnames(odorStruct);
                
                % Getting the odorID like Priscilla did
                odorID = extractBetween(odorList(iOdor), "I ", " -");
    
                color = color_ids(find(odor_ids == str2double(odorID), 1), :);
    
                % This has all the ROIs so we get a slice below
                allFrames_allROIs = odorStruct.hit;
    
                % Lines are a little less transparent and thicker than
                % in other graphs (0.5 to 0.7 for both of them).
                plot(xs', allFrames_allROIs(:, iROI), ...
                        'Color', [color 0.7], 'LineWidth', 2);
    
                disp(strcat("plot odor ", odorID, " done"))
            end
    
            hold off;
    
            % Plot settings
            set(gca, 'FontName', 'Arial');
            set(gcf, 'OuterPosition', [100 100 1000 900]);
            set(gca, 'LineWidth', 0.75);
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
    
        end
    
        % Last tile showing the ROI
        nexttile(2, [rows, columns - 1])
        
        hold on
        imshow(imadjust(cast(zProjFileToAnalyze,'uint16')))
    
        thetas = linspace(0,2*pi,200);
        ellipseR1 = (ROIs{iROI}.vnRectBounds(4) - ROIs{iROI}.vnRectBounds(2))/2;
        ellipseR2 = (ROIs{iROI}.vnRectBounds(3) - ROIs{iROI}.vnRectBounds(1))/2;
        ellipseA = (ROIs{iROI}.vnRectBounds(4) + ROIs{iROI}.vnRectBounds(2))/2;
        ellipseB = (ROIs{iROI}.vnRectBounds(3) + ROIs{iROI}.vnRectBounds(1))/2;
        ellipseX = ellipseR1 * cos(thetas) + ellipseA;
        ellipseY = ellipseR2 * sin(thetas) + ellipseB; 
        plot(ellipseX, ellipseY, 'Color', 'y', 'LineWidth', 1);
    
        hold off  
    
        t.TileSpacing = 'compact';
        t.Padding = 'compact';
        disp(strcat("plot roi ", num2str(iROI), " done"))
    end
end

%% PLOT AVGS - hits and false choices (with ranges), color-coded by outcome

allProgramFieldName = fieldnames(s_zscore);

% (number of odors) x (number of programs + extra columns)
% The extra columns are to put the plot showing the ROI
extraColumns = 3;

rows = length(allProgramFieldName);
columns = maxOdorNumber + extraColumns;

% Struct with the outcome colors, for easy reference
outcomeColors = struct( ...
    'hit', [56 83 163]/255, ...
    'false_choice', [236 30 36]/255, ...
    'miss',  [127 127 127]/255 ...
    );

% xsx is the palindromic version of xs (needed for the fill plot)
xsx = [xs flip(xs)];

for iROI = 1:nROIs
    % Create one figure per ROI
    figName = strcat( ...
        firstAcquisitionName, '_to_', lastAcquisitionName, ...
        '_roi_', num2str(iROI), '_zdFF_stderr');

    fig = figure('Name', figName);

    % Plot settings
    set(gca, 'FontName', 'Arial');
    set(gcf, 'OuterPosition', [100 100 1600 900]);
    set(gca, 'LineWidth', 0.75);

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

            axis([xmin xmax ymin ymax])
            title(strcat(odorFieldName, '_', programType), ...
                'Interpreter', 'none');
            xlabel('Time from odor onset (s)')
            ylabel('dF/F Z-score')
            
            % Iterate through outcomes, plotting both on the same graph
            for iOutcome = 1:length(allOutcomeFieldName)
                outcomeFieldName = allOutcomeFieldName{iOutcome};

                % Get the color from struct
                color = outcomeColors.(outcomeFieldName);

                % Gets data for this ROI and all acquisitions
                allFrames_allROIs_allAcq = odorStruct.(outcomeFieldName);
                allFrames_thisROI_allAcq = ...
                    allFrames_allROIs_allAcq(:,iROI,:);

                % Get the mean and stderr along acquisitions
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
                    'Color', [color 0.7], 'LineWidth', 0.7);
            end

            hold off;
            disp(strcat("plot odor ", odorID, " done"))
        end
    end

    % Last tile showing the ROI
    nexttile(columns - extraColumns + 1, [rows, extraColumns])
    
    hold on
    imshow(imadjust(cast(zProjFileToAnalyze,'uint16')))

    thetas = linspace(0,2*pi,200);
    ellipseR1 = (ROIs{iROI}.vnRectBounds(4) - ROIs{iROI}.vnRectBounds(2))/2;
    ellipseR2 = (ROIs{iROI}.vnRectBounds(3) - ROIs{iROI}.vnRectBounds(1))/2;
    ellipseA = (ROIs{iROI}.vnRectBounds(4) + ROIs{iROI}.vnRectBounds(2))/2;
    ellipseB = (ROIs{iROI}.vnRectBounds(3) + ROIs{iROI}.vnRectBounds(1))/2;
    ellipseX = ellipseR1 * cos(thetas) + ellipseA;
    ellipseY = ellipseR2 * sin(thetas) + ellipseB; 
    plot(ellipseX, ellipseY, 'Color', 'y', 'LineWidth', 1);

    hold off  

    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    disp(strcat("plot roi ", num2str(iROI), " done"))
end

%% Save figs

if saveFigs
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    
    % save all open figs
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName = FigList(iFig).Name;
      set(0, 'CurrentFigure', FigHandle);
      % forces matlab to save fig as a vector
      FigHandle.Renderer = 'painters';  
      % actually saves a vector file
      saveas(FigHandle,fullfile(saveDir, [FigName '.svg']));
    end 
    disp('saved all figs')
    close all
end


%% Save workspace

if saveWorkspace
    % save workspace variables
    matFileName = strcat(imgsToAnalyzeNames{1}(1:end-9),'_',imgsToAnalyzeNames{end}(end-13:end-4),'_preProcessing');
    save(fullfile(saveDir,matFileName));     
    disp('saved mat file')
end