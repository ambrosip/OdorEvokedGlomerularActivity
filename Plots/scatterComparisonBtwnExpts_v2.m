%% USER INPUT

ylimit.dFF = 0.1;
ylimit.zScore = 5;
label.dFF = "dF/F";
label.zScore = "Z-score";

% ALERT: hard coded comparison between exp 1 and exp 2 or program 1 and
% program 2


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

                scatterNicePlot(namePrefix,...
                    strcat("odor ", num2str(odorID)),...
                    expData(expDataNum).name,...
                    dataType,...
                    label,-ylimit.(dataType),ylimit.(dataType),...
                    mean(allOdorDataExp1(thisOdorExp1AcqIdx',:)),...
                    mean(allOdorDataExp2(thisOdorExp2AcqIdx',:)))
            end
        end

        % == BASELINE PLOT ================================================
        scatterNicePlot(namePrefix,...
                    "baseline",...
                    expData(expDataNum).name,...
                    dataType,...
                    label,-ylimit.(dataType),ylimit.(dataType),...
                    baselineAvgsToCompare(:,1),...
                    baselineAvgsToCompare(:,2))  
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

        scatterNicePlot(namePrefix,...
                    "baseline",...
                    expData(1).name,...
                    dataType,...
                    label,-ylimit.(dataType),ylimit.(dataType),...
                    baselinesToCompare(:,1),...
                    baselinesToCompare(:,2))  
    
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

            scatterNicePlot(namePrefix,...
                    strcat("odor ", num2str(odorID)),...
                    expData(1).name,...
                    dataType,...
                    label,-ylimit.(dataType),ylimit.(dataType),...
                    mean(odorAll(thisOdorProg1AcqIdx',:)),...
                    mean(odorAll(thisOdorProg2AcqIdx',:)))
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


%% Save workspace

% save workspace variables
matFileName = strcat(expData(1).name, "_to_", ...
     expData(end).name,'_sequentialZScore');
save(fullfile(saveDir,matFileName));     
disp('saved mat file')
