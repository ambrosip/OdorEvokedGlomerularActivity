% DEPENDS ON OUTPUTS FROM:
% 1) image_dF_post_odor (*_mcor_postOdor.mat)

% TODO: 1) Split blobs
%       2) Plot z-scores instead of dF/F
%       3) Make it work with hits too.

%% USER INPUT

% Range for the z-score plots
plotRange = [-5, 5];

% Lower threshold to exclude from the z-score mask
lowerThreshold = 1.0;

% Disks for erosion and dilation
erosionDisk = strel("disk", 5);
dilationDisk = strel("disk", 5);

%% Generating Mask

fig = figure('Name', 'Z-Score Mask Steps');
tl = tiledlayout(2,3);

minEnvelope = zeros(size(figures(1).na.image));
maxEnvelope = zeros(size(figures(1).na.image));

maxBlur = zeros(size(figures(1).na.image));

for iFigure = 1:length(figures)
    minEnvelope = min(minEnvelope, figures(iFigure).na.image);
    maxEnvelope = max(maxEnvelope, figures(iFigure).na.image);
    maxBlur = max(maxBlur, abs(imgaussfilt(figures(iFigure).na.image, 3.5)));
end

% Plotting the z-score envelope
zScoreEnvelope = (abs(minEnvelope) >= abs(maxEnvelope)) .* minEnvelope + ...
        (abs(minEnvelope) < abs(maxEnvelope)) .* maxEnvelope;

nexttile;
imshow(zScoreEnvelope, plotRange);
colormap(divergingGradient);

% Lower threshold for z-score
ROI = abs(maxBlur) > lowerThreshold;

nexttile;
imshow(maxBlur, plotRange);
colormap(divergingGradient);

% Distance from boundary
distanceMask = - bwdist(~ROI);

% Make regions more smooth to split less
smoothedMask = medfilt2(distanceMask);

nexttile;
imshow(smoothedMask, []);

% Split regions
splitMask = watershed(smoothedMask, 4);
splitMask(~ROI) = 0.0;

% Color regions independently
splitMaskRGB = label2rgb(splitMask,'jet', [.5, .5, .5], 'shuffle');

nexttile;
imshow(splitMaskRGB);

% Make regions rounder
erodedMask = imerode(splitMask, erosionDisk);
dilatedMask = imdilate(erodedMask, dilationDisk);

% Plot eroded
nexttile;
imshow(erodedMask, []);

% Plot final result with initial blurred max projection
nexttile;
imshow(maxBlur, plotRange);
colormap(divergingGradient);

hold on;
finalMaskRGB = label2rgb(dilatedMask, 'jet', 'k', 'shuffle'); 

hOverlay = imshow(finalMaskRGB);
alphaMap = double(dilatedMask > 0) * 0.5; 
set(hOverlay, 'AlphaData', alphaMap);
hold off;

% Fix spacing
tl.TileSpacing = 'tight';
tl.Padding = 'tight';

fig.Position = [30 30 1450 800];

% Rearrange regions numbers to make sequential
% TODO: do this in a more efficient way
[finalMask, nROIs] = bwlabel(dilatedMask > 0);

% ADAPTED FROM TimeSeriesFromFijiROIsZscores (TODO: CLEAN UP)
%% Get fluorescence in fiji ROIs for each file/acquisition

% meanInt = mean intensity
meanIntPerRoi = [];
for file = 1:imgsToAnalyze_numberOf

    % get OS-appropriate file dir
    imgToAnalyzeFileDir = fullfile(imgsToAnalyzeFolder, imgsToAnalyzeDirs(file).name);

    % get img file name without extension (stored in "f")
    [p,f,e] = fileparts(imgsToAnalyzeNames(file));

    % add "a" in front of the filename in f to build a structure later
    % why? matlab freaks out if field names of a structure start with numbers
    f = {strcat('a', cell2mat(f))};
    
    % get img info
    imgInfo = imfinfo(imgToAnalyzeFileDir);
    frames_per_img = length(imgInfo);
    
    % iterate frame by frame (each frame is a time point)
    for frame = 1:frames_per_img

        % Read img frame by frame
        % Raw files from scanimage and mcor files from the matlab motion
        % correction script are stored as int16 (signed 16-bit integer,
        % with values ranging from -32768 to 32767). To avoid errors with
        % dF/F calculation due to negative values, I used to convert these
        % files to uint16 (unsigned 16-bit integer, with values ranging
        % from 0 to 65535). To allow further processing, I then converted
        % the numbers to single (similar to double, but with less precision). 
        % Mcor files from the python script are stored as int32 (signed
        % 32-bit integer, with values ranging from -2,147,483,648 to
        % 2,147,483,647), but the actual range of the data is just a bit
        % wider than the int16 range. To deal with all of these files
        % without improperly compressing the data (e.g. the function cast
        % does NOT work properly to convert int32 to uint16), I decided to
        % simply shift the dataset to the right by subtracting the min
        % value if it's below zero. If the min value is above zero, we do
        % nothing, because we don't have negative values to worry about. In
        % sum, regardless of bit depth, we will make sure that the lower 
        % range of values is >= zero. FYI, imread saves data as double.
        imgToAnalyze = imread(imgToAnalyzeFileDir,frame);

        % ------------- COMMENTED OUT FOR TESTING --------------------- %
        % Make sure that all entries of imgToAnalize are >= zero
        % if min(imgToAnalyze, [], 'all') < 0
        %     imgToAnalyze = imgToAnalyze - min(imgToAnalyze, [], 'all');
        % end
        % ------------------------------------------------------------- %

        % % old code, kept for archiving purposes
        % % deal with 32-bit files - this does NOT WORK - DO NOT USE
        % if convertFrom32bit
        %     imgToAnalyze = cast(imgToAnalyze, 'uint16');
        % end
        % 
        % % convert img to uint16 (range: 0 to 65535)
        % imgToAnalyze = im2uint16(imgToAnalyze);
        
        % convert img to single
        imgToAnalyze = single(imgToAnalyze);
    
        % iterate ROI by ROI
        for roiNumber = 1:nROIs
            labeledRoi = finalMask == roiNumber;
            % figure; imshow(labeledRoi) % use this to troubleshoot 
            nPixelsInRoi = sum(labeledRoi,'all');
            % labeledRoiAsInt16 = int16(labeledRoi);
            labeledRoiAsSingle = single(labeledRoi);
            % maskedImg = labeledRoiAsInt16.*imgToAnalyze;
            maskedImg = labeledRoiAsSingle.*imgToAnalyze;
            % figure; imshow(imadjust(maskedImg,[0.5 0.65])) % use this to troubleshoot 
            % it is safe to sum uint16 variables: https://www.mathworks.com/matlabcentral/answers/5401-matlab-function-mean-returns-the-exact-same-value-for-uint16-and-double-values-not-for-single
            meanIntInRoi = sum(maskedImg,'all')/nPixelsInRoi;
            % store mean fluorescence per frame and roi
            meanIntPerRoi(frame,roiNumber) = meanIntInRoi;
        end
    end

    % store info for all files in a structure
    s.(f{1})=meanIntPerRoi;

    disp(strcat("processing rois in file ", f, " done"))
end


%% Calculate dF/F for each ROI across files

% set default firstFig and lastFig boundaries in case user does NOT want a
% custom subset
if plotSubset == 0
    firstAcq = 1;
    lastAcq = imgsToAnalyze_numberOf;
end

fns = fieldnames(s);
firstAcqName = fns{firstAcq};
lastAcqName = fns{lastAcq};

% delete the first few seconds of the data because of photobleaching
photobleaching_window_frames = photobleaching_window_s * frame_rate_hz;
adjusted_baseline_dur_s = baseline_dur_s - photobleaching_window_s;
adjusted_baseline_frames = adjusted_baseline_dur_s * frame_rate_hz;

% calculate dF/F
% dF/F = (F - mean F in baseline) / mean F in baseline
dF_per_file=[];
for file=firstAcq:lastAcq
    f_per_file = s.(fns{file});
    f_per_file = f_per_file(photobleaching_window_frames:end,:);
    for roi=1:nROIs
        mean_baseline_f = mean(f_per_file(1:adjusted_baseline_frames,roi),'omitnan');
        dF_per_file(:,roi) = (f_per_file(:,roi) - mean_baseline_f) / mean_baseline_f;
    end
    s_dF.(fns{file})=dF_per_file(1:end,:);
end

% create x axis in seconds and adjust it based on the photobleaching window
adjusted_img_dur_dataPts = size(dF_per_file,1);
adjusted_img_dur_s = adjusted_img_dur_dataPts / frame_rate_hz;
adjusted_odor_onset_s = adjusted_baseline_dur_s;
adjusted_odor_offset_s = adjusted_baseline_dur_s + odor_dur_s;
xAxisInSec = linspace(0 - adjusted_odor_onset_s, adjusted_img_dur_s - adjusted_odor_onset_s, adjusted_img_dur_dataPts);

disp("calculated dF/F")

%% PLOT data in ROIs organized by odor (rows) and program type (columns)

% get max and min value cof dF/F to set y axis limits
% ymax = round(max(structfun(@(x) max(x,[],'all'),s_dF,'UniformOutput',true)), TieBreaker='plusinf');
% ymin = round(min(structfun(@(x) min(x,[],'all'),s_dF,'UniformOutput',true)), TieBreaker='minusinf');

% Gets the first integer below min and first above max
ymax = ceil(max(structfun(@(x) max(x,[],'all'),s_dF,'UniformOutput',true)));
ymin = floor(min(structfun(@(x) min(x,[],'all'),s_dF,'UniformOutput',true)));

% get max and min value of xAxis to set x acis limits
xmin = round(min(xAxisInSec),TieBreaker='minusinf');
xmax = round(max(xAxisInSec),TieBreaker='plusinf');

% adjust ymin to -0.1 in case it's zero
if ymin == 0
   ymin = -0.1;
end

% adjust max to +0.1 in case it's zero
if ymax == 0
   ymax = 0.1;
end

% get max number of odors used in this experiment
max_odor_num = 0;
for programNum = 1:size(programFieldNames,1) % ALERT added ",1" on Mar 19 2026
    programFieldName = programFieldNames(programNum);
    if s_olfactometer.(programFieldName).type ~= "ignore"
        odor_num = size(s_olfactometer.(programFieldName).odorList,1);
        if max_odor_num < odor_num
            max_odor_num = odor_num;
        end
    end
end

% iterate through rois
ignoredPrograms = size(programFieldNames,1) - programsToAnalyze;
labelRGB = label2rgb(finalMask, 'jet', 'k', 'shuffle'); 

for roi=1:nROIs
    figName = strcat(firstAcqName(2:end), '_to_', lastAcqName(2:end), '_roi_', num2str(roi), '_dF');
    fig = figure('Name',figName);
    set(gca,'FontName','Arial');
    set(gca,'LineWidth', 0.75);
    % make a layout with "odor" rows and "program" + 1 columns
    rows = max_odor_num;
    columns = programsToAnalyze+1;
    t = tiledlayout(rows,columns);
    title(t,figName,'Interpreter','none');
    relativeProgramNum = 0;
    for programNum = 1:size(programFieldNames,1) % ALERT added ,1 
        programFieldName = programFieldNames(programNum);
        if s_olfactometer.(programFieldName).type ~= "ignore"
            relativeProgramNum = relativeProgramNum + 1;
            for odorNum = 1:length(s_olfactometer.(programFieldName).odorList)
                nexttile(relativeProgramNum + (odorNum-1)*columns)
                odorID = extractBetween(s_olfactometer.(programFieldName).odorList(odorNum),"I "," -");
                odorFieldName = s_olfactometer.(programFieldName).odorFieldNames(odorNum);
                color = odor_color.colorID(odor_color.odorID==str2double(odorID),:);
                hold on;
                rectangle('Position',[0 ymin odor_dur_s ymax-ymin],'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
                yline(0,'k--')
                axis([xmin xmax ymin ymax])
                title(strcat(odorFieldName, '_', s_olfactometer.(programFieldName).type), 'Interpreter','none');
                xlabel('Time from odor onset (s)')
                ylabel('dF/F')

                for acqIdx = s_olfactometer.(programFieldName).summary_by_trial.acqIdx(s_olfactometer.(programFieldName).summary_by_trial.odor==str2double(odorID))'
                    % plot dF/F   
                    if ~isnan(acqIdx)
                        plot(xAxisInSec',s_dF.(fns{acqIdx})(:,roi),'Color',[color 0.5],'LineWidth',0.5);  
                    end
                end

                % plot mean dF/F
                % plot(xAxisInSec',s_mean_dF.(odorFieldNames(odor))(:,roi),'Color',color,'LineWidth',1);
                hold off;

                disp(strcat("plot odor ", odorID, " done"))
            end
        end
    end
    % show ROI location
    nexttile(columns,[max_odor_num,1])
    if convertFrom32bit
        zProjFileToAnalyze = cast(zProjFileToAnalyze,'uint16');
        imshow(imadjust(zProjFileToAnalyze))
    else
        imshow(imadjust(zProjFileToAnalyze,[0.5 0.65])) 
    end


    hold on

    hOverlay = imshow(labelRGB);
    alphaMap = double(finalMask == roi) * 0.5; 
    set(hOverlay, 'AlphaData', alphaMap);

    hold off  
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    fig.Position = [0 0 800 1200];

    drawnow;

    disp(strcat("plot roi ", num2str(roi), " done"))
end

disp("plot done")

%% Save Figures

todayStr = string(datetime('Today'), 'yyyy-MM-dd');
saveFolder = fullfile(expDir, 'processed', 'matlab', todayStr);

if ~isfolder(saveFolder)
    mkdir(saveFolder);
end

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
      saveas(FigHandle,fullfile(saveFolder, [FigName '.svg']));
    end 

    disp('saved all figs')
    close all
end
