%{ 
scatterCompareBtwnPrograms

ALERT: hard coded comparison between program 1 and program 2!

2 ways of using this code:
    1) LOAD sequential mat file containing data for prog 1 and prog 2
    2) LOAD masksFromGUI mat file containing data for prog 2 and specify path
         to masksFromGUI mat file containing data for prog 1. 

GOAL: compare odor-evoked responses for all odors in program 1 vs program 2
    - 1 plot per odor, each dot is an ROI
    - average accross acquisitions
    - NOT sequential: zscores and dFF are relative to pre-odor baseline,
        NOT baseline of first aquisition
%}


%% USER INPUT

saveFigs = 1;
saveDir = "/Volumes/T7 Shield/PA/Poster/scatter_NOT_sequential";

% deal with instances where data was collected in separate loops
loadExp1Data = 0;
exp1MatFileDir = '/Volumes/T7 Shield/PA/Poster/20260306m357/20260306_m357_e1_00001_to_20260306_m357_e1_00090_2026-04-16_masksFromGUI.mat';

% use color coding or gray
useOdorColorCoding = 0;

% define what data to plot 
dataType = "zScore"; % dFF or zScore
outcome = "na"; % use "na" for passive odor presentations

% plot aesthetics
ylimit_min.dFF = -0.1;
ylimit_max.dFF = 0.1;
ylimit_min.zScore = -5;
ylimit_max.zScore = 5;
label.dFF = "dF/F";
label.zScore = "Z-score";

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


%% LOAD DATA AND PLOT

namePrefix = strcat(figNameStart, "_to_", figNameMiddle, "_", figNameEnd,...
    "_NOT-sequential");

if loadExp1Data == 0
    figureTitle = expData(1).name;
    % get odorIDs as numbers
    allOdorIDs = unique(expData(1).db_trials.odorID)';
else
    % remove acq number for fig name
    split_figNameStart = split(figNameStart,"_");
    figureTitle = join(split_figNameStart(1:end-1),'_');
    % get odorIDs as numbers
    allOdorIDs = unique(db_trials.odorID)';
end

for odorID = allOdorIDs

    odorFieldName = strcat('odor_',num2str(odorID));

    if useOdorColorCoding == 1
        color = odor_color.colorID(odor_color.odorID==odorID,:);
    else
        color = grayColor;
    end

    if loadExp1Data == 0
        % hard coding program 1 and program 2
        if dataType == "zScore"
            xAxisData = s_mean_zscore.program_1.(odorFieldName).(outcome);
            yAxisData = s_mean_zscore.program_2.(odorFieldName).(outcome);
        else
            xAxisData = s_mean_dFF.program_1.(odorFieldName).(outcome);
            yAxisData = s_mean_dFF.program_2.(odorFieldName).(outcome);
        end
    else
        % hard coding program 1 for both experiments
        if dataType == "zScore"
            % exp 1, prog 1 is xAxis
            exp1Data = load(exp1MatFileDir,'s_mean_zscore');
            xAxisData = exp1Data.s_mean_zscore.program_1.(odorFieldName).(outcome);
            % exp 2, prog 1 is yAxis
            yAxisData = s_mean_zscore.program_1.(odorFieldName).(outcome);
            % saniy check
            disp(strcat("prog 1 has ",...
                num2str(size(exp1Data.s_mean_zscore.program_1.(odorFieldName).(outcome),1)),...
                " frames and ",...
                num2str(size(exp1Data.s_mean_zscore.program_1.(odorFieldName).(outcome),2)),...
                " rois"));
            disp(strcat("prog 2 has ",...
                num2str(size(s_mean_zscore.program_1.(odorFieldName).(outcome),1)),...
                " frames and ",...
                num2str(size(s_mean_zscore.program_1.(odorFieldName).(outcome),2)),...
                " rois"));
        else
            % exp 1, prog 1 is xAxis
            exp1Data = load(exp1MatFileDir,'s_mean_dFF');
            xAxisData = exp1Data.s_mean_dFF.program_1.(odorFieldName).(outcome);
            % exp 2, prog 1 is yAxis
            yAxisData = s_mean_dFF.program_1.(odorFieldName).(outcome);
            % sanity check
            disp(strcat("prog 1 has ",...
                num2str(size(exp1Data.s_mean_dFF.program_1.(odorFieldName).(outcome),1)),...
                " frames and ",...
                num2str(size(exp1Data.s_mean_dFF.program_1.(odorFieldName).(outcome),2)),...
                " rois"));
            disp(strcat("prog 2 has ",...
                num2str(size(s_mean_dFF.program_1.(odorFieldName).(outcome),1)),...
                " frames and ",...
                num2str(size(s_mean_dFF.program_1.(odorFieldName).(outcome),2)),...
                " rois"));
        end
    end

    scatterNicePlot(namePrefix,...
                    strcat("odor ", num2str(odorID)),...
                    figureTitle,...
                    dataType,...
                    label,ylimit_min.(dataType),ylimit_max.(dataType),...
                    mean(xAxisData),...
                    mean(yAxisData),...
                    color)
end


%% PRE-ODOR BASELINE



%% SAVE

if saveFigs == 1
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