%{ 
DOCUMENTATION
written by PA on Aug/2025

GOAL:
    Calculate and show dF/F in each ROI (manually drawn as ellipses in
    fiji) organized by odor (rows and color) and program type (columns).

IMPORTANT re fiji ROIs:
    Do not create overlapping ROIs because the fiji ROI code I'm using will
    not properly analyze them

ASSUMPTIONS:
    !! search for ALERT and ASSUMPTION to read important info

DEPENDS on:
    mat file created by preProcessing.m
    From others > ReadImageJROI (https://www.mathworks.com/matlabcentral/fileexchange/32479-readimagejroi)
    From others > ScanImageTiffReader

TO DO: 
    check if code still works for passive odor presentations (probably not cuz fml)
    re-calculate mean_dF
%}


%% Get dir of fiji files 

roi_file_name = dir(fullfile(expDir,'processed','fiji','*.zip')).name;
roi_file_folder = dir(fullfile(expDir,'processed','fiji','*.zip')).folder;
roi_file_dir = fullfile(roi_file_folder,roi_file_name);

% Get fiji Avg Intensity Projection file (zProj = z Projection)
zProjFileDirs = dir(fullfile(expDir, 'processed', 'fiji', '*.tif'));
zProjFileNames = {zProjFileDirs.name}';
zProjFileFolders = {zProjFileDirs.folder}';
zProj_numberOf = length(zProjFileNames);
% ALERT: only using the 1st .tif file found
zProjFileToAnalyzeDir = fullfile(zProjFileFolders{1}, zProjFileNames{1});
zProjFileToAnalyze = imread(zProjFileToAnalyzeDir);
% % quality control of user input
% if size(zProjFileToAnalyze) == size(mcorImgToAnalyze)
%     disp('zProj and mcor imgs match in size')
% else
%     disp('WARNING: zProj and mcor imgs DO NOT match in size')
% end

% Get fiji ROI file (roi = regions of interest)
roisFileDirs = dir(fullfile(expDir, 'processed', 'fiji', '*.zip'));
roisFileNames = {roisFileDirs.name}';
roisFileFolders = {roisFileDirs.folder}';
% ALERT: only analyzing the 1st .zip file found
roisFileToAnalyze = fullfile(roisFileFolders{1}, roisFileNames{1});
rois = ReadImageJROI(roisFileToAnalyze);
% vnImageSize is one of the inputs for ROIs2Regions
% for whatever reason, it is a translated version of the image size
vnImageSize=[size(zProjFileToAnalyze,2),size(zProjFileToAnalyze,1)];
regions=ROIs2Regions(rois,vnImageSize);
rois_numberOf = length(rois);

disp('loaded fiji data')


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
        imgToAnalyze = imread(imgToAnalyzeFileDir,frame);
    
        % iterate ROI by ROI
        for roiNumber = 1:rois_numberOf
            labeledRoi = labelmatrix(regions) == roiNumber;
            labeledRoi = labeledRoi';
            % figure; imshow(labeledRoi) % use this to troubleshoot 
            nPixelsInRoi = sum(labeledRoi,'all');
            labeledRoiAsInt16 = int16(labeledRoi);
            maskedImg = labeledRoiAsInt16.*imgToAnalyze;
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
    for roi=1:rois_numberOf
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


%% Calculate mean dF/F for each odor and ROI

% for programNum = 1:size(programFieldNames)
%         programFieldName = programFieldNames(programNum);
%         if s_olfactometer.(programFieldName).type ~= "ignore"
%             for odorNum = 1: length(s_olfactometer.(programFieldName).odorList)
%                 odorID = extractBetween(s_olfactometer.(programFieldName).odorList(odorNum),"I "," -");
%                 odorFieldName = s_olfactometer.(programFieldName).odorFieldNames(odorNum);


% % mean dF/F for each odor and ROI across files
% for odor = 1:length(odorFieldNames)
%     mean_dF_per_ROI = [];
%     for roi=1:rois_numberOf  
%         dF_per_file=[];
%         for acq = struct_acqNum_by_odor.(odorFieldNames(odor))' 
%             dF_per_file = [dF_per_file s_dF.(fns{acq})(:,roi)];
%         end
%         mean_dF_per_ROI(:,roi) = mean(dF_per_file, 2 ,'omitnan');        
%     end
%     s_mean_dF.(odorFieldNames(odor)) = mean_dF_per_ROI;
% end


%% PLOT data in ROIs organized by odor (rows) and program type (columns)

% get max and min value of dF/F to set y axis limits
ymax = round(max(structfun(@(x) max(x,[],'all'),s_dF,'UniformOutput',true)), TieBreaker='plusinf');
ymin = round(min(structfun(@(x) min(x,[],'all'),s_dF,'UniformOutput',true)), TieBreaker='minusinf');

% get max and min value of xAxis to set x acis limits
xmin = round(min(xAxisInSec),TieBreaker='minusinf');
xmax = round(max(xAxisInSec),TieBreaker='plusinf');

% adjust ymin to some negative number in case it's zero
if ymin == 0
   ymin = -1;
end

% get max number of odors used in this experiment
max_odor_num = 0;
for programNum = 1:size(programFieldNames)
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
for roi=1:rois_numberOf
    figName = strcat(firstAcqName(2:end), '_to_', lastAcqName(2:end), '_roi_', num2str(roi), '_dF');
    fig = figure('Name',figName);
    set(gca,'FontName','Arial');
    set(gcf,'OuterPosition',[100 100 1200 900]); % [left bottom width height]
    set(gca,'LineWidth', 0.75);
    % make a layout with "odor" rows and "program" + 1 columns
    rows = max_odor_num;
    columns = programsToAnalyze+1;
    t = tiledlayout(rows,columns);
    title(t,figName,'Interpreter','none');
    for programNum = 1:size(programFieldNames)
        programFieldName = programFieldNames(programNum);
        if s_olfactometer.(programFieldName).type ~= "ignore"
            for odorNum = 1:length(s_olfactometer.(programFieldName).odorList)
                nexttile(programNum - ignoredPrograms + (odorNum-1)*columns)
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
    imshow(imadjust(zProjFileToAnalyze,[0.5 0.65])) 
    hold on
    thetas = linspace(0,2*pi,200);
    ellipseR1 = (rois{roi}.vnRectBounds(4) - rois{roi}.vnRectBounds(2))/2;
    ellipseR2 = (rois{roi}.vnRectBounds(3) - rois{roi}.vnRectBounds(1))/2;
    ellipseA = (rois{roi}.vnRectBounds(4) + rois{roi}.vnRectBounds(2))/2;
    ellipseB = (rois{roi}.vnRectBounds(3) + rois{roi}.vnRectBounds(1))/2;
    ellipseX = ellipseR1*cos(thetas)+ellipseA;
    ellipseY = ellipseR2*sin(thetas)+ellipseB; 
    plot(ellipseX,ellipseY,'Color','y','LineWidth',1);
    hold off  
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    disp(strcat("plot roi ", num2str(roi), " done"))
end

disp("plot done")


%% Save figs

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
disp('I saved the figs')
close all


%% Save workspace

% save workspace variables
matFileName = strcat(imgsToAnalyzeNames{1}(1:end-9),'_',imgsToAnalyzeNames{end}(end-13:end-4),'_timeSeriesFromFijiROIs');
save(fullfile(saveDir,matFileName));     
disp('I saved the mat file')


%% ARCHIVE

% %% Calculate z score for each ROI across files

% % z-score = (dF/F - mean(dF/F) in baseline) / sd(dF/F) in baseline
% fns = fieldnames(s_dF);
% zScorePerFile=[];
% for file=1:mcor_numberOf
%     dFPerFile = s_dF.(fns{file});
%     for roi=1:rois_numberOf
%         meanBaseline_dF = mean(dFPerFile(1:ceil(adjustedBaselineInFrames),roi),'omitnan');
%         sdBaseline_dF = std(dFPerFile(1:ceil(adjustedBaselineInFrames),roi),'omitnan');
%         zScorePerFile(:,roi) = (dFPerFile(:,roi) - meanBaseline_dF) / sdBaseline_dF;
%     end
%     s_zS.(fns{file})=zScorePerFile(1:end,:);
% end
% 
% % mean z-score in ROI across files
% fns = fieldnames(s_zS);
% zScorePerFile=[];
% for roi=1:rois_numberOf
%     zSPerROI = [];
%     for file=1:mcor_numberOf 
%         zScorePerFile(:,file) = s_zS.(fns{file})(:,roi);
%     end
%     mean_zS_PerROI(:,roi) = mean(zScorePerFile,2,'omitnan');
%     if filter == 1
%         mean_zS_PerROI_filtered(:,roi) = smooth(mean_zS_PerROI(:,roi),span);
%     end
% end


% for file=firstAcq:lastAcq      
%     % plot dF/F            
%     plot(xAxisInSec',s_dF.(fns{file})(:,roi));
%     if addAnnotations == 1
%         annotatedFile = ismember(filteredLog.acqNum,file);
%         if sum(annotatedFile) >= 1
%             description = filteredLog.description(annotatedFile);
%             acqStartTime = s_xAxis.(fns{file})(1);
%             xline(acqStartTime,'b--',description, 'LabelVerticalAlignment','top')
%         end
%     end
% end


% %% code hoarding

% for roi=1:rois_numberOf
%     figName = strcat(firstAcqName(2:end), '_to_', lastAcqName(2:end), '_roi_', num2str(roi), '_dF');
%     fig = figure('Name',figName);
%     set(gca,'FontName','Arial');
%     set(gcf,'OuterPosition',[100 100 1200 900]); % [left bottom width height]
%     set(gca,'LineWidth', 0.75);
%     % make a layout with "odor" rows and 2 columns
%     t = tiledlayout(length(odorList),2);
%     for odor = 1:length(odorList)
%         nexttile(2*odor - 1)
%         odorNum = str2double(extractBetween(odorList(odor),"I "," -"));
%         color = odor_color.colorID(odor_color.odorID==odorNum, :);        
%         hold on;
%         rectangle('Position',[0 ymin odor_dur_s ymax-ymin],'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
%         yline(0,'k--')
%         axis([xmin xmax ymin ymax])
%         title(strcat(odorFieldNames(odor), '_', figName), 'Interpreter','none');
%         xlabel('Time from odor onset (s)')
%         ylabel('dF/F')
%         for acq = struct_acqNum_by_odor.(odorFieldNames(odor))'
%             % plot dF/F            
%             plot(xAxisInSec',s_dF.(fns{acq})(:,roi),'Color',[color 0.5],'LineWidth',0.5);            
%         end
%         % plot mean dF/F
%         plot(xAxisInSec',s_mean_dF.(odorFieldNames(odor))(:,roi),'Color',color,'LineWidth',1);
%         hold off;
%         disp(strcat("plot odor ", num2str(odor_ids(odor)), " done"))
%     end
%     % show ROI location
%     nexttile(2,[odor,1])
%     imshow(imadjust(zProjFileToAnalyze,[0.5 0.65])) 
%     hold on
%     thetas = linspace(0,2*pi,200);
%     ellipseR1 = (rois{roi}.vnRectBounds(4) - rois{roi}.vnRectBounds(2))/2;
%     ellipseR2 = (rois{roi}.vnRectBounds(3) - rois{roi}.vnRectBounds(1))/2;
%     ellipseA = (rois{roi}.vnRectBounds(4) + rois{roi}.vnRectBounds(2))/2;
%     ellipseB = (rois{roi}.vnRectBounds(3) + rois{roi}.vnRectBounds(1))/2;
%     ellipseX = ellipseR1*cos(thetas)+ellipseA;
%     ellipseY = ellipseR2*sin(thetas)+ellipseB; 
%     plot(ellipseX,ellipseY,'Color','y');
%     hold off  
%     t.TileSpacing = 'compact';
%     t.Padding = 'compact';
% 
%     disp(strcat("plot roi ", num2str(roi), " done"))
% end