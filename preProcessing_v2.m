%% USER INPUT - experiment directory and others - EDIT ME

% experiment dir to be analyzed
expDir = '/Users/priscilla/Documents/Local - Moss Lab/20250902/e3';

% set img-specific inputs
photobleaching_window_s = 2; % duration of data in senconds that will be removed from baseline to account for photobleaching
analyzeMcorImgs = 1; % 0 (no) or 1 (yes)
plotSubset = 0; % 0 (no) or 1 (yes)
    firstAcq = 1; % ignore if plotSubset == 0; edit if plotSubset == 1
    lastAcq = 3; % ignore if plotSubset == 0; edit if plotSubset == 1

% label relevant events from olfactometer program that are not
% automatically found
trial_start_label = "Output 1";

% label relevant events from h5 file
% ALERT: ImagingWindow was used as trial_start in some files and as
% imaging_window in other files!
trial_start_label_h5 = '/ImagingWindow';
odor_dur_h5 = '/OdorDelivery';

% set program type x odor x action x outcome relationships
olfactory_task = "2afc_fine_coarse_fine";
olfactory_task = "passive_odor_presentations";
minLicksToTriggerReward = 3;

% set max time interval (in min) between file acq time and trial onset 
% trigger to assign acq numbers to specific trials when you accidentally 
% miss the acq of some trials
tolerance = 0.1;

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
odor_ids = [1; 2; 17; 18; 19; 20];
color_ids = [black_color; reddish_purple_color; blue_color; sky_blue_color; vermillion_color; orange_color];
odor_color = table(odor_ids,color_ids,'VariableNames',{'odorID','colorID'});

disp('loaded inputs')


%% Pre-processing

%%
getFileDirs
%%
getImgMetadata
%%
getOlfactometerData
%%
getScopeH5Timestamps
%%
matchAcqsToTrials


%% Calculate baseline duration 

% test if you missed the first trial peak, causing the size of odor_locs to
% be larger than the size of trial_locs
if size(trial_locs,1) < size(odor_locs,1)
    % drop off the first odor_locs and the first odor_end_locs
    odor_locs = odor_locs(2:end);
    odor_end_locs = odor_end_locs(2:end);
end
baseline_dur_s = mean((odor_locs - trial_locs)*60);
img_dur_s = frames_per_img / frame_rate_hz;
baseline_dur_frames = baseline_dur_s * frame_rate_hz;
odor_dur_s = mean((odor_end_locs - odor_locs)*60);
odor_dur_frames = odor_dur_s * frame_rate_hz;
odor_onset_s = baseline_dur_s;
odor_offset_s = odor_onset_s + odor_dur_s;
disp('calculated baseline')


%% FIG 2 - OLFACTOMETER ODOR PRESENTATIONS RASTER

fig2 = figure('name', strcat(s_olfactometer.program_1.shortName, '_', analysisDate, ' - raster'));
t=tiledlayout('horizontal');
for programNum = 1:size(programFieldNames,1)
    programFieldName = programFieldNames(programNum);
    if s_olfactometer.(programFieldName).type ~= "ignore"  
        nexttile
        yticks_nums_and_labels = ["1", "Trials"];
        plot(s_olfactometer.(programFieldName).startMin_by_trial,1,'|','Color','k','LineWidth',1)
        hold on;
        for odorNum = 1:length(s_olfactometer.(programFieldName).odorList)
            odorID = extractBetween(s_olfactometer.(programFieldName).odorList(odorNum),"I "," -");
            odorFieldName = s_olfactometer.(programFieldName).odorFieldNames(odorNum);
            color = odor_color.colorID(odor_color.odorID==str2double(odorID),:);
            plot(s_olfactometer.(programFieldName).(odorFieldName).startMin_by_odor,odorNum+1,'|','LineWidth',1,'Color',color)
            yticks_nums_and_labels = [yticks_nums_and_labels; [odorNum, odorID]];
        end
        hold off;
        yMinForRaster = 0;
        yMaxForRaster = length(s_olfactometer.(programFieldName).odorList)+2;
        % yMaxForRaster = 5;
        xMinForRaster = 0;
        xMaxForRaster = ceil(max(s_olfactometer.(programFieldName).summary_by_trial.min)/5)*5; % round up to nearest multiple of 5
        axis([xMinForRaster xMaxForRaster yMinForRaster yMaxForRaster])
        yticks(1:odorNum+1);
        yticklabels(yticks_nums_and_labels(:,2)')
        ylabel('Odors')
        xticks([0,xMaxForRaster]);
        xlabel('Time (min)');
        title(s_olfactometer.(programFieldName).type, 'Interpreter','none');
    end
end 

set(fig2, 'Position', [0 0 500 250])    % x y width height
disp('plot fig2 complete')


%% FIG 3 - plot task performance

if olfactory_task == "2afc_fine_coarse_fine"

    hits_per_block = [];
    misses_per_block = [];
    false_choices_per_block = [];
    trials_per_block = [];
    program_types = "";

    for programNum = 1:size(programFieldNames,1)
        programFieldName = programFieldNames(programNum);
        if s_olfactometer.(programFieldName).type ~= "ignore" 
            hits_per_block = [hits_per_block; sum(strcmp(s_olfactometer.(programFieldName).outcome_by_trial,"hit"))];
            misses_per_block = [misses_per_block; sum(strcmp(s_olfactometer.(programFieldName).outcome_by_trial,"miss"))];
            false_choices_per_block = [false_choices_per_block; sum(strcmp(s_olfactometer.(programFieldName).outcome_by_trial,"false choice"))];
            trials_per_block = [trials_per_block; size(s_olfactometer.(programFieldName).outcome_by_trial,1)];
            program_types = strcat(program_types, '_', s_olfactometer.(programFieldName).type);
        end
    end

    fig3 = figure('name', strcat(s_olfactometer.program_1.shortName, '_', analysisDate, ' - performance'));
    hold on;
    plot(hits_per_block./trials_per_block,'Color',bluish_green_color,'LineWidth',1,'DisplayName','Hits');
    plot(false_choices_per_block./trials_per_block,'Color',vermillion_color,'LineWidth',1,'DisplayName','False choices');
    plot(misses_per_block./trials_per_block,'Color',black_color,'LineWidth',1,'DisplayName','Misses');
    yline(0.5,'--')
    axis([1 programsToAnalyze 0 1])
    yticks([0,1]);
    xticks([0,programsToAnalyze]);
    xlabel('Block');
    ylabel('Events/Trials');
    title([s_olfactometer.program_1.shortName, strcat("analyzed_on_", analysisDate), strip(program_types,'_'), ""], 'Interpreter','none');
    hold off;
    legend('Hits', 'False choices', 'Misses','Location','northwest')
    legend('boxoff')
    set(fig3, 'Position', [100 100 200 300])    % x y width height
    disp('plot fig3 complete')
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
matFileName = strcat(imgsToAnalyzeNames{1}(1:end-9),'_',imgsToAnalyzeNames{end}(end-13:end-4),'_preProcessing');
save(fullfile(saveDir,matFileName));     
disp('saved mat file')
