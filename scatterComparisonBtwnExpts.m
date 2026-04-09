%% USER INPUT

ymax = 0.8;
ymin = -ymax;



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

% if comparisonType == "between experiments"
%     % check for mismatching number of acqs between expts
%     if min(expData.nFiles) ~= max(expData.nFiles)
%         disp("discarding extra acqs in one expt")
%     end
%     nFilesToCompare = min(expData.nFiles);
%     baselinesToCompare = nan(expData.nROIs, length(expData));
%     for expDataNum = 1:length(expData)
%         baselinesToCompare(:,expDataNum) = 

if comparisonType == "between programs"

    % == BASELINE =========================================================
    % average accorss frames
    baselinesAll = mean(dFF(expData.baselineStart:expData.baselineEnd,:,:));
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
    figure('Name',strcat(namePrefix, "_baselines"))
    scatter(baselinesToCompare(:,1),baselinesToCompare(:,2))
    hold on
    plot([ymin,ymax],[ymin,ymax],'k')
    hold off
    axis([ymin,ymax,ymin ymax])
    ylabel('dF/F post-treatment')
    xlabel('dF/F pre-treatment')
    title('Pre-odor signal per ROI',expData(1).name,'Interpreter','none')
    axis square

    %== ODOR-EVOKED =======================================================
    % average across frames
    odorAll = mean(dFF(expData.odorStart:expData.odorEnd,:,:));
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
        figure('Name',strcat(namePrefix, "_odor_", num2str(odorID)))
        nexttile
            for roi = 1:expData.nROIs
                scatter(odorAll(thisOdorProg1AcqIdx',roi),odorAll(thisOdorProg2AcqIdx',roi))
                hold on;
            end   
            plot([ymin,ymax],[ymin,ymax], 'k')
            axis([ymin,ymax,ymin ymax])
            ylabel('dF/F post-treatment')
            xlabel('dF/F pre-treatment')
            title('Odor-evoked signal per ROI',strcat(expData(1).name, "_odor_", num2str(odorID)),'Interpreter','none')
            axis square
        nexttile
            scatter(mean(odorAll(thisOdorProg1AcqIdx',:)),mean(odorAll(thisOdorProg2AcqIdx',:)),"*")  
            hold on;
            plot([ymin,ymax],[ymin,ymax], 'k')
            hold off;
            axis([ymin,ymax,ymin ymax])
            ylabel('dF/F post-treatment')
            xlabel('dF/F pre-treatment')
            title('Mean odor-evoked signal per ROI',strcat(expData(1).name, "_odor_", num2str(odorID)),'Interpreter','none')
            axis square
            set(gcf, 'Position', [50 50 700 350])    % x y width height
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
save(fullfile(saveDir,matFileName));     
disp('saved mat file')
