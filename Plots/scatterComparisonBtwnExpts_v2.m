%% USER INPUT

% LOAD sequential mat file

ylimit_min.dFF = -0.1;
ylimit_max.dFF = 0.1;
ylimit_min.zScore = -5;
ylimit_max.zScore = 5;
label.dFF = "dF/F";
label.zScore = "Z-score";

olfactory_task = "paint me gray";
% olfactory_task = "passive_odor_presentations";

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

% odor color-coding
odor_ids = [1; 2; 17; 18; 19; 20;...
    21; 22; 23; 24;...
    25; 26; 27; 28;...
    3];
color_ids = [black_color; reddish_purple_color; blue_color; sky_blue_color; vermillion_color; orange_color;...
    blue_color; sky_blue_color; vermillion_color; orange_color;...
    blue_color; sky_blue_color; vermillion_color; orange_color;...
    black_color];
odor_color = table(odor_ids,color_ids,'VariableNames',{'odorID','colorID'});

grayColor = [.7 .7 .7];




%% CODE

% check if this is a comparison between multiple experiments or between
% multiple programs within an experiment (multiple expDirs vs single
% expDir)
if length(expData) > 1
    comparisonType = "between experiments";
elseif length(unique(expData.db_trials.programNum)) > 1
    comparisonType = "between programs";
else
    disp("bruh wtf u trying to compare. there is only 1 exp and 1 program here")
    return;
end

namePrefix = strcat(...
    expData(1).name, "_to_", ...
    expData(end).name,"_acq_",...
    num2str(expData(end).nFiles), '_sequential');

if comparisonType == "between experiments"

    % check for mismatching number of acqs between expts
    if min(expData.nFiles) ~= max(expData.nFiles)
        disp("discarding extra acqs in one expt")
    end
    % crop number of acqs that will be analyzed
    nFilesToCompare = min(expData.nFiles);

    if expData(1).nROIs ~= expData(2).nROIs
        disp("CRITICAL FLAG: different number of ROIs in first 2 expts")
    end

    for dataType = ["dFF", "zScore"]

        % create nan array that will be filled
        % "acq" dim 1, "ROI" dim 2 and "exp" dim 3
        % ASSUMPTION: all expts have the same number of ROIs
        baselinesToCompare = nan(nFilesToCompare, expData(1).nROIs, length(expData));

        for expDataNum = 1:length(expData)
    
            % adjust first data point to 1 instead of 0
            expData(expDataNum).baselineStart = expData(expDataNum).baselineStart + 1;     

            dataToAnalyze = expData(expDataNum).(dataType);
        
            % == BASELINE =================================================
            % average across frames
            baselinesAll = mean(dataToAnalyze(expData(expDataNum).baselineStart:expData(expDataNum).baselineEnd,:,:));
            % creates array with "acq" rows and "ROI" columns
            baselinesAll = squeeze(baselinesAll); 
            % crop out extra acqs
            baselinesAll = baselinesAll(1:nFilesToCompare,:);
            % save data for this exp
            baselinesToCompare(:,:,expDataNum) = baselinesAll;
            % mean across acqs, "roi" rows, "exp" columns
            baselineAvgsToCompare = mean(baselinesToCompare);
            baselineAvgsToCompare = squeeze(baselineAvgsToCompare);

            % check if we have multiple programs inside an experiment
            if length(unique(expData(expDataNum).db_trials.programNum)) > 1
                disp('we have multiple programs inside multiple experiments, thanks for the mess!')
                disp('Im not gonna deal w this rn')

                % % started editing
                % % create array with "ROI" rows and "program" columns
                % baselinesToCompare = nan(expData(expDataNum).nROIs,length(unique(expData(expDataNum).db_trials.programNum)));
                % % create arrays with 1 row and "program" columns
                % lastAcqIdx = nan(1,length(unique(expData(expDataNum).db_trials.programNum)));
                % firstAcqIdx = nan(1,length(unique(expData(expDataNum).db_trials.programNum)));
                % totalAcqs = nan(1,length(unique(expData(expDataNum).db_trials.programNum)));
                % for programNum = 1:length(unique(expData(expDataNum).db_trials.programNum))
                %     allAcqIdx = expData(expDataNum).db_trials.acqIdx(expData(expDataNum).db_trials.programNum == programNum,:);
                %     lastAcqIdx(1,expDataNum,programNum) = max(allAcqIdx);
                %     firstAcqIdx(1,expDataNum,programNum) = min(allAcqIdx(allAcqIdx > 0));
                %     totalAcqs(1,expDataNum,programNum) = lastAcqIdx(1,expDataNum,programNum) - firstAcqIdx(1,expDataNum,programNum) + 1;
                %     baselinesToCompare(:,expDataNum,programNum) = mean(baselinesAll(firstAcqIdx(1,expDataNum,programNum):lastAcqIdx(1,expDataNum,programNum),:))';
                % end               

                % % did not even start editing
                % if totalAcqs(1) < totalAcqs(2)
                %     disp("we have more acqs in exp2 bro")
                %     lastAcqIdx(2) = lastAcqIdx(2) - (totalAcqs(2) - totalAcqs(1)); % CHECK MATH ALERT
                % elseif totalAcqs(1) > totalAcqs(2)
                %     disp("we have more acqs in exp1 bro")
                %     lastAcqIdx(1) = lastAcqIdx(1) - (totalAcqs(1) - totalAcqs(2));
                % else
                %     disp("good job, we have the same number of acqs in both expts")
                % end
            end
        end

        dataToAnalyzeExp1 = expData(1).(dataType);
        dataToAnalyzeExp2 = expData(2).(dataType);

        % == ODOR-EVOKED ==============================================
        % hard coded comparison between exp 1 and exp 2
        % average across frames
        allOdorDataExp1 = mean(dataToAnalyzeExp1(expData(1).odorStart:expData(1).odorEnd,:,:));
        allOdorDataExp2 = mean(dataToAnalyzeExp2(expData(2).odorStart:expData(2).odorEnd,:,:));
        % create array with "acq" rows and "ROI" columns
        allOdorDataExp1 = squeeze(allOdorDataExp1);
        allOdorDataExp2 = squeeze(allOdorDataExp2);

        % get data from odors that exist in both exp 1 and exp 2
        for odorID = unique(expData(1).db_trials.odorID)'
            if ismember(odorID, unique(expData(2).db_trials.odorID)')
                thisOdorExp1AcqIdx = expData(1).db_trials.acqIdx(expData(1).db_trials.odorID == odorID,:);
                thisOdorExp2AcqIdx = expData(2).db_trials.acqIdx(expData(2).db_trials.odorID == odorID,:);
                if length(thisOdorExp1AcqIdx) ~= length(thisOdorExp2AcqIdx)
                    nAcqsToCompare = min(length(thisOdorExp1AcqIdx), length(thisOdorExp2AcqIdx));
                    thisOdorExp1AcqIdx = thisOdorExp1AcqIdx(1:nAcqsToCompare,:);
                    thisOdorExp2AcqIdx = thisOdorExp2AcqIdx(1:nAcqsToCompare,:);
                end

                if olfactory_task == "passive_odor_presentations"
                    % Get the odor color from user input in preProcessing_v2
                    color = odor_color.colorID(odor_color.odorID==str2double(odorID),:);
                else
                    color = grayColor;
                end

                scatterNicePlot(namePrefix,...
                    strcat("odor ", num2str(odorID)),...
                    expData(expDataNum).name,...
                    dataType,...
                    label,ylimit_min.(dataType),ylimit_max.(dataType),...
                    mean(allOdorDataExp1(thisOdorExp1AcqIdx',:)),...
                    mean(allOdorDataExp2(thisOdorExp2AcqIdx',:)),...
                    color)
            end
        end

        % == BASELINE PLOT ================================================
        if olfactory_task == "passive_odor_presentations"
            % Get the odor color from user input in preProcessing_v2
            color = 'k';
        else
            color = grayColor;
        end

        scatterNicePlot(namePrefix,...
                    "baseline",...
                    expData(expDataNum).name,...
                    dataType,...
                    label,ylimit_min.(dataType),ylimit_max.(dataType),...
                    baselineAvgsToCompare(:,1),...
                    baselineAvgsToCompare(:,2),...
                    color)  
    end

elseif comparisonType == "between programs"

    for dataType = ["dFF", "zScore"]

        % adjust first data point to 1 instead of 0
        expData.baselineStart = expData.baselineStart + 1;

        dataToAnalyze = expData.(dataType);
        
        % == BASELINE =========================================================
        % average accorss frames
        baselinesAll = mean(dataToAnalyze(expData.baselineStart:expData.baselineEnd,:,:));
        % creates array with as many rows as acquisitions and as many columns
        % as ROIs
        baselinesAll = squeeze(baselinesAll);
        % create array with "ROI" rows and "program" columns
        baselinesToCompare = nan(expData.nROIs,length(unique(expData.db_trials.programNum)));
        % create arrays with 1 row and "program" columns
        lastAcqIdx = nan(1,length(unique(expData.db_trials.programNum)));
        firstAcqIdx = nan(1,length(unique(expData.db_trials.programNum)));
        totalAcqs = nan(1,length(unique(expData.db_trials.programNum)));
        for programNum = 1:length(unique(expData.db_trials.programNum))
            allAcqIdx = expData.db_trials.acqIdx(expData.db_trials.programNum == programNum,:);
            lastAcqIdx(1,programNum) = max(allAcqIdx);
            firstAcqIdx(1,programNum) = min(allAcqIdx(allAcqIdx > 0));
            totalAcqs(1,programNum) = lastAcqIdx(1,programNum) - firstAcqIdx(1,programNum) + 1;
            baselinesToCompare(:,programNum) = mean(baselinesAll(firstAcqIdx(1,programNum):lastAcqIdx(1,programNum),:))';
        end
        if totalAcqs(1) < totalAcqs(2)
            disp("we have more acqs in exp2 bro")
            lastAcqIdx(2) = lastAcqIdx(2) - (totalAcqs(2) - totalAcqs(1)); % CHECK MATH ALERT
        elseif totalAcqs(1) > totalAcqs(2)
            disp("we have more acqs in exp1 bro")
            lastAcqIdx(1) = lastAcqIdx(1) - (totalAcqs(1) - totalAcqs(2));
        else
            disp("good job, we have the same number of acqs in both expts")
        end

        if olfactory_task == "passive_odor_presentations"
            % Get the odor color from user input in preProcessing_v2
            color = 'k';
        else
            color = grayColor;
        end

        scatterNicePlot(namePrefix,...
                    "baseline",...
                    expData(1).name,...
                    dataType,...
                    label,ylimit_min.(dataType),ylimit_max.(dataType),...
                    baselinesToCompare(:,1),...
                    baselinesToCompare(:,2),...
                    color)  
    
        %== ODOR-EVOKED =======================================================
        % average across frames
        odorAll = mean(dataToAnalyze(expData.odorStart:expData.odorEnd,:,:));
        % creates array with as many rows as acquisitions and as many columns
        % as ROIs
        odorAll = squeeze(odorAll);
        for odorID = unique(expData.db_trials.odorID)'
            thisOdorAllAcqIdx = expData.db_trials.acqIdx(expData.db_trials.odorID == odorID,:);
            thisOdorProg1AcqIdx = thisOdorAllAcqIdx(thisOdorAllAcqIdx > 0 & thisOdorAllAcqIdx <= lastAcqIdx(1,1));
            thisOdorProg2AcqIdx = thisOdorAllAcqIdx(thisOdorAllAcqIdx > 0 & thisOdorAllAcqIdx > lastAcqIdx(1,1));
            if length(thisOdorProg1AcqIdx) ~= length(thisOdorProg2AcqIdx)
                nAcqsToCompare = min(length(thisOdorProg1AcqIdx), length(thisOdorProg2AcqIdx));
                thisOdorProg1AcqIdx = thisOdorProg1AcqIdx(1:nAcqsToCompare,:);
                thisOdorProg2AcqIdx = thisOdorProg2AcqIdx(1:nAcqsToCompare,:);
            end  

            if olfactory_task == "passive_odor_presentations"
                % deal with mismatch in odorID type in code written by me
                % vs Vinicius
                if isnumeric(odorID)
                    color = odor_color.colorID(odor_color.odorID==odorID,:);
                else
                    color = odor_color.colorID(odor_color.odorID==str2double(odorID),:);
                end
            else
                color = grayColor;
            end

            scatterNicePlot(namePrefix,...
                    strcat("odor ", num2str(odorID)),...
                    expData(1).name,...
                    dataType,...
                    label,ylimit_min.(dataType),ylimit_max.(dataType),...
                    mean(odorAll(thisOdorProg1AcqIdx',:)),...
                    mean(odorAll(thisOdorProg2AcqIdx',:)),...
                    color)
        end
    end
end

        
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
disp('saved all figs')
close all


% %% Save workspace
% 
% % save workspace variables
% matFileName = strcat(expData(1).name, "_to_", ...
%      expData(end).name,'_sequentialZScore');
% save(fullfile(saveDir,matFileName));     
% disp('saved mat file')
