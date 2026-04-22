%{

created to compare dFF to zscore images and fix zscore calculation!
image_dF_v4 calculates z = (x-mean)/std where
    x = mean dF/F across time for one pixel
    mean = mean dF/F across time and pixels
    std = mean dF/F across time and pixels
this calculation compares each pixel to all other pixels instead of
comparing each pixel to its past self

1) Run pre-processing OR load preProcessing MAT files (in case variables are not in the environment)

order of things
    1) get F over time - collect frames from key time windows (baseline, 
    odor presentation, or both)
    2) get F over window - calculate avg across frames for key time windows
    (baseline, odor, or both)
    3) calculate dFF, zScore1 and zScore2 for odor presentation window
        dFF: use baseline mean
        zScore1: use baseline mean
        zScore2: use baseline + odor mean
    4) avg across acquisitions

%}


%% USER INPUT

% set expDir to adjust dirs in case you're loading a MAT file
previousFile = 1;
expDir = "M:\ImagingData\20260311\m316";

% Put NaN for automatic limits
absoluteLimit = 1;

% Define percentiles
LOWER_QUANTILE = 0.01;
UPPER_QUANTILE = 0.99;

% Define the colors
max_df_color = [103 0 31] / 255;
min_df_color = [5 48 97] / 255;

% Plot one figure for each program and odor combination or plot only hits
% for all programs and odors in a single figure?
plotOnlyHits = false;

compareFigTypes = 1;

% if odor presentation is 1s and you choose divideWindowsByFactorOf = 2,
% you will compare 0.5 s of odor presentation to 0.5 s of baseline
% divideWindowsByFactorOf = 1;


%% Extra inputs in case you run this before timeSeriesFromFijiROIs

% set default firstFig and lastFig boundaries in case user does NOT want a
% custom subset
if plotSubset == 0
    firstAcq = 1;
    lastAcq = imgsToAnalyze_numberOf;
end
firstAcqName = imgsToAnalyzeDirs(1).name;
lastAcqName = imgsToAnalyzeDirs(lastAcq).name;


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

% if exist('divideWindowsByFactorOf','var')
%     frameSignalDuration = round(frameOdorDuration / divideWindowsByFactorOf);
% end


%% Update dirs

if previousFile == 1
    if exist('patchwarp','var')
        if patchwarp == 1
            getFileDirs_patchwarp
        else
            getFileDirs
        end
    else
        getFileDirs
    end
end
    
% Open the first frame of the first file to get image dimensions
filename = imgsToAnalyzeDirs(1).name;
fileDir = imgsToAnalyzeDirs(1).folder;
filepath = fullfile(fileDir, filename);
frameSize = size(read_file(filepath, 1, 1));


%% Iterate through programs

% Blank image for outcomes with no corresponding acquisition
blankImage = zeros(frameSize);

% There is one figure for each combination of program and odor
% iFigure is the index of the figure currently being computed
iFigure = 0;

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
        % Add one figure to the figure list and record its program type,
        % odor number, and number of acquisitions
        iFigure = iFigure + 1;

        figures(iFigure).type = programType;
        figures(iFigure).odor = odor;
        figures(iFigure).acquisitions = 0;
        figures(iFigure).programNumber = programNumber;

        % Get table with data only for the current odor
        odorRows = summaryByTrial.odor == odor;
        odorTable = summaryByTrial(odorRows, :);

        % Initialize counts of how many images where used in the average
        figures(iFigure).hit.total = 0;
        figures(iFigure).miss.total = 0;
        figures(iFigure).false.total = 0;
        figures(iFigure).na.total = 0;

        % Initialize images
        for outcomeType = ["hit", "miss", "false", "na"]
            for figType = ["image_dFF", "image_zScore_v1", "image_zScore_v2"]
                figures(iFigure).(outcomeType).(figType) = zeros(frameSize);
            end
        end

        for row = 1:height(odorTable)
            acqIdx = odorTable{row, 'acqIdx'};

            % Skip trials without acquisition
            if isnan(acqIdx) || acqIdx == 0
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

            % Load and computes the mean of the relevant frames

            % Note that since read_file returns a matrix of integers from
            % [-32768, 32767] we need to add 32768 to that matrix to shift its
            % values to the [0, 65535] range.

            % get 
            % Load frames for the baseline window
            % this will collect the value of each pixel (x y) for all frames (time) during the baseline window into a (x y time) matrix
            baselineFrames = single(read_file( ...
                filepath, frameBaselineStart, frameOdorDuration)) + 32768;
            % calculate the average accross time (aka frames) for each pixel into a (x y) matrix
            baselineImageMean = mean(baselineFrames, ndims(baselineFrames)); 
            % calculate the std accross time (aka frames) for each pixel into a (x y) matrix
            baselineImageStd = std(baselineFrames, 0, ndims(baselineFrames));

            % Load frames for the signal window 
            % this will collect the value of each pixel (x y) for all
            % frames (time) during the odor presentation into a (x y time)
            % matrix
            signalFrames = single(read_file( ...
                filepath, frameOdorOnset, frameOdorDuration)) + 32768;
            % calculate the average across time for each pixel
            signalImageMean = mean(signalFrames, ndims(signalFrames));

            % Load frames for baseline window + odor presentation window
            allFrames = single(read_file( ...
                filepath, frameBaselineStart, 2*frameOdorDuration)) + 32768;
            % take the average across time for each pixel
            allFramesMean = mean(allFrames, ndims(allFrames));
            % calculate the std across time for each pixel
            allFramesStd = std(allFrames, 0, ndims(allFrames));
            
            % calculate the dF/F for each pixel using the mean during the
            % baseline window
            dFF = ...
                (signalImageMean - baselineImageMean) ./ baselineImageMean;

            % calculate the zScore for each pixel using the mean and std of
            % the baseline window
            zScore_v1 = ...
                (signalImageMean - baselineImageMean) ./ baselineImageStd;

            % take the zScore for each pixel using the mean and std of the
            % baseline & odor presentation window
            zScore_v2 = ...
                (signalImageMean - allFramesMean) ./ allFramesStd;

            % outcome is the first word of the outcome name
            % Possibilities are "hit", "miss", "false", and "na".
            outcome = split(odorTable{row, 'outcome'});
            outcome = outcome(1);

            % Increments counters and create a copy (with a shorter name)
            figures(iFigure).acquisitions = ...
                figures(iFigure).acquisitions + 1;
            figures(iFigure).(outcome).total = ...
                figures(iFigure).(outcome).total + 1;

            n = figures(iFigure).(outcome).total;

            % Average across aquisitions
            % We use a simple trick to compute the average on the fly:
            % https://stackoverflow.com/a/23493727
            figures(iFigure).(outcome).image_dFF = ...
                figures(iFigure).(outcome).image_dFF * (n - 1) / n + ...
                dFF / n;

            figures(iFigure).(outcome).image_zScore_v1 = ...
                figures(iFigure).(outcome).image_zScore_v1 * (n - 1) / n + ...
                zScore_v1 / n;

            figures(iFigure).(outcome).image_zScore_v2 = ...
                figures(iFigure).(outcome).image_zScore_v2 * (n - 1) / n + ...
                zScore_v2 / n;
        end

        % for outcome = ["hit" "miss" "false" "na"]
        %     % REMOVE THIS CALCULATION (this code uses average accross pixels) 
        %     % figures(iFigure).(outcome).image = ...
        %     %     (figures(iFigure).(outcome).image - ...
        %     %     mean(figures(iFigure).(outcome).image(:))) / ...
        %     %     std(figures(iFigure).(outcome).image(:));
        % 
        %     figures(iFigure).(outcome).image_dFF = figures(iFigure).(outcome).image_dFF;
        %     figures(iFigure).(outcome).image_zScore_v1 = figures(iFigure).(outcome).image_zScore_v1; 
        %     figures(iFigure).(outcome).image_zScore_v2 = figures(iFigure).(outcome).image_zScore_v2;
        % end

        % Don't compute if there are no acquisition files
        % If there were no acquisition files nBaseline is zero
        if figures(iFigure).acquisitions == 0
            message = ...
                'WARNING: No files for program %d (%s) and odor %d.\n';
            fprintf(message, programNumber, programType, odor);
            continue
        end

        message = 'Analyzed %d files for program %d (%s) and odor %d.\n';
        fprintf(message, figures(iFigure).acquisitions, ...
            programNumber, programType, odor);
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

% set the color scale limits if user did not specify an absolute limit 
if isnan(absoluteLimit)
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
            quantile(figures(i).hit.image_zScore_v1(:), LOWER_QUANTILE), ...
            quantile(figures(i).miss.image_dFF(:), LOWER_QUANTILE), ...
            quantile(figures(i).false.image_dFF(:), LOWER_QUANTILE), ...
            quantile(figures(i).na.image_dFF(:), LOWER_QUANTILE) ...
            ]);

        upperLimit = max([ ...
            upperLimit, ...
            quantile(figures(i).hit.image_zScore_v1(:), UPPER_QUANTILE), ...
            quantile(figures(i).miss.image_dFF(:), UPPER_QUANTILE), ...
            quantile(figures(i).false.image_dFF(:), UPPER_QUANTILE), ...
            quantile(figures(i).na.image_dFF(:), UPPER_QUANTILE) ...
            ]);
    end

    % Take the limit that is larger in absolute value. It also turns it
    % into a double because that is a requirement for plots.
    absoluteLimit = double(max(abs(upperLimit), abs(lowerLimit)));
end

% Make range symmetrical around zero (so white means zero in the plots).
% The same range will be used for all figures, for easy comparisons.
plotRange = [-absoluteLimit absoluteLimit];

if plotOnlyHits % ignoring this for now!
    % % Get figure name
    % % Get the start of firstAcqName (before the third underline)
    % figNameStart = split(string(firstAcqName), '_');
    % figNameStart = join(figNameStart(1:end-1), '_');
    % 
    % % Get the number of the last acquisition
    % figNameEnd = split(string(lastAcqName), '_');
    % figNameEnd = figNameEnd(end-1);
    % 
    % % Join results in the correct format
    % figName = sprintf("%s_to_%s", figNameStart, figNameEnd);
    % 
    % fig = figure('Name', figName);
    % 
    % % Get the number of unique programs and odors
    % programTypes = unique([figures.type], 'stable');
    % odors = unique([figures.odor]);
    % 
    % tl = tiledlayout(length(odors), length(programTypes));
    % title(tl, figName, 'Interpreter', 'none', 'FontSize', 16);
    % 
    % % Add blank images for pairs of program and odor that don't have a
    % % corresponding figure. This is needed to make sure column and row labels
    % % appear in the final figure.
    % for row = 1:length(odors)
    %     for col = 1:length(programTypes)
    %         % Check if there is a figure for this program and odor
    %         rows = [figures.odor] == odors(row);
    %         cols = strcmp([figures.type], programTypes(col));
    % 
    %         if sum(rows & cols) == 0
    %             % If there is no figure, add a blank image
    %             nexttile((row - 1) * length(programTypes) + col)
    %             imshow(blankImage, plotRange);
    %         end
    %     end
    % end
    % 
    % % Plot all figures
    % for iFigure = 1:length(figures)
    %     % Get the row and column of the current figure
    %     row = find(odors == figures(iFigure).odor);
    %     col = find(strcmp(programTypes, figures(iFigure).type));
    % 
    %     nexttile((row - 1) * length(programTypes) + col)
    %     imshow(figures(iFigure).hit.image, plotRange);
    % end
    % 
    % % Create column labels
    % for column = 1:length(programTypes)
    %     ax = nexttile(column);
    %     ax.XAxisLocation = 'top';
    %     xlabel(programTypes(column), 'FontSize', 12);
    % end
    % 
    % % Create row labels
    % for row = 1:length(odors)
    %     ax = nexttile((row - 1) * length(programTypes) + 1);
    %     ax.YAxisLocation = 'left';
    %     label = sprintf('Odor %d', odors(row));
    %     ylabel(label, 'FontSize', 12);
    % end
    % 
    % fig.Position = [200 100 frameSize(2)/2.5 * length(programTypes) frameSize(1)/2.5 * length(odors) + 50];
    % 
    % % Uses the colormap created at the start
    % cb = colorbar;
    % cb.Layout.Tile = 'south';
    % 
    % colormap(divergingGradient);
    % clim(plotRange);
    % 
    % % ylabel(cb, 'dF/F', 'FontSize', 12);
    % ylabel(cb, 'Z-Score', 'FontSize', 12);
    % 
    % tl.TileSpacing = 'tight';
    % tl.Padding = 'tight';
    % 
    % drawnow;

elseif compareFigTypes == 1
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
        figNameBase = sprintf("%s_to_%s_program_%d_odor_%d_%s", figNameStart, ...
            figNameMiddle, figures(iFigure).programNumber, ...
            figures(iFigure).odor, figNameEnd);

        figName = strcat(figNameBase, "_", "comparison");
        
        fig = figure('Name', figName);

        % 3 fig types in a row
        tl = tiledlayout('horizontal');
        title(tl, figName, 'Interpreter', 'none');
        nPlots = 3;

        for figType = ["image_dFF", "image_zScore_v1", "image_zScore_v2"]
            legendName = split(figType, "_");
            legendName = legendName(2);
            nexttile
            imshow(figures(iFigure).na.(figType), plotRange)
            title(legendName, 'FontSize', 16)
        end

        baseSize = 600;
        fig.Position = [200 100 baseSize * nPlots baseSize + 150];

        % Uses the colormap created at the start
        cb = colorbar;
        cb.Layout.Tile = 'south';

        colormap(divergingGradient);
        clim(plotRange);

        % ylabel(cb, '$dF/F$', 'interpreter', 'latex', 'FontSize', 16);
        % ylabel(cb, 'dF/F', 'FontSize', 16);
        ylabel(cb, legendName, 'FontSize', 16);

        tl.TileSpacing = 'compact';
        tl.Padding = 'compact';

        drawnow;
    end
else % also ignoring for now
    % for iFigure = 1:length(figures)
    %     % Get figure name
    %     % Get the start of firstAcqName (before the third underline)
    %     figNameStart = split(string(firstAcqName), '_');
    %     figNameStart = join(figNameStart(1:end-1), '_');
    % 
    %     % Get the number of the last acquisition
    %     figNameMiddle = split(string(lastAcqName), '_');
    %     figNameMiddle = figNameMiddle(end-1);
    % 
    %     % Reformats program type string from "Fine 1" to "fine_1", for example
    %     figNameEnd = join(split(lower(figures(iFigure).type)), '_');
    % 
    %     % Join results in the correct format
    %     figNameBase = sprintf("%s_to_%s_program_%d_odor_%d_%s", figNameStart, ...
    %         figNameMiddle, figures(iFigure).programNumber, ...
    %         figures(iFigure).odor, figNameEnd);
    % 
    %     for figType = ["image_dFF", "image_zScore_v1", "image_zScore_v2"]
    % 
    %         figName = strcat(figNameBase, "_", figType);
    %         legendName = split(figType, "_");
    %         legendName = legendName(2);
    % 
    %         fig = figure('Name', figName);
    % 
    %         % 3 possible outcomes in a row (or just 1 if there is NA)
    %         tl = tiledlayout('horizontal');
    %         title(tl, figName, 'Interpreter', 'none');
    % 
    %         if figures(iFigure).na.total == 0
    %             nPlots = 3;
    % 
    %             nexttile
    %             imshow(figures(iFigure).hit.(figType), plotRange)
    %             title('Hits', 'FontSize', 16)
    % 
    %             nexttile
    %             imshow(figures(iFigure).false.(figType), plotRange)
    %             title('False Choices', 'FontSize', 16)
    % 
    %             nexttile
    %             imshow(figures(iFigure).miss.(figType), plotRange)
    %             title('Misses', 'FontSize', 16)
    %         else
    %             nPlots = 1;
    % 
    %             nexttile
    %             imshow(figures(iFigure).na.(figType), plotRange)
    %             title('NA', 'FontSize', 16)
    %         end
    % 
    %         baseSize = 600;
    %         fig.Position = [200 100 baseSize * nPlots baseSize + 150];
    % 
    %         % Uses the colormap created at the start
    %         cb = colorbar;
    %         cb.Layout.Tile = 'south';
    % 
    %         colormap(divergingGradient);
    %         clim(plotRange);
    % 
    %         % ylabel(cb, '$dF/F$', 'interpreter', 'latex', 'FontSize', 16);
    %         % ylabel(cb, 'dF/F', 'FontSize', 16);
    %         ylabel(cb, legendName, 'FontSize', 16);
    % 
    %         tl.TileSpacing = 'compact';
    %         tl.Padding = 'compact';
    % 
    %         drawnow;
    %     end
    % end
end


%% Save figs

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

% save all open figs
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = FigList(iFig).Name;
    set(0, 'CurrentFigure', FigHandle);
    % forces matlab to save fig as a vector
    % FigHandle.Renderer = 'painters';
    % actually saves a vector file
    saveas(FigHandle, fullfile(saveDir, [FigName '.png']));
end
disp('saved all figs')
close all


%% Save workspace

% save workspace variables
matFileName = strcat(imgsToAnalyzeNames{1}(1:end-9),'_',imgsToAnalyzeNames{end}(end-13:end-4),'_preProcessing');
save(fullfile(saveDir,matFileName));     
disp('saved mat file')