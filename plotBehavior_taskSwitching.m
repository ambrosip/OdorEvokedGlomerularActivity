%{ 
DOCUMENTATION

run analyzeBehavior before using this script

%}

function plotBehavior_taskSwitching(date, db_trials, saveDir, mouseNum)

close all;

%% USER INPUT

% % set data to analyze
% dateToPlot = "2025_10_24";
dateToPlot = date;

% fine-tune data display
exclude_First_N_Trials = 0;
trials_per_bin = 5; % trials per session

% ALERT HARD CODED FOR EASE, change to SOFT CODE later
odorDur_s = 1;
responseDur_s = 3;
totalDur_s = odorDur_s + responseDur_s;

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

% outcome color-coding
hit_color = [56 83 163]/255;
false_choice_color = [236 30 36]/255;
miss_color = [127 127 127]/255;

% odor colo-coding
odor17_color = blue_color;
odor18_color = sky_blue_color;
odor19_color = vermillion_color;


%% FILTER DATA BY DATE & SANITY CHECK

% filter data by date based on user input
filtered_rows_by_date = matches(db_trials.date, dateToPlot);

% figure out if multiple programs were started that day
programNumList = unique(db_trials(filtered_rows_by_date,:).programNum);
totalNumOfPrograms = size(programNumList,1);
first_programNum = programNumList(1);
last_programNum = programNumList(end);

% sanity checks that this is a day of task-switching in the order expected
% by the code
if totalNumOfPrograms == 4    
    % let's brute code for sake of time 
    first_program_rows = (db_trials.programNum == programNumList(1));
    second_program_rows = (db_trials.programNum == programNumList(2));
    third_program_rows = (db_trials.programNum == programNumList(3));
    fourth_program_rows = (db_trials.programNum == programNumList(4));

    first_program_name = db_trials(first_program_rows,:).programName(1);
    second_program_name = db_trials(second_program_rows,:).programName(1);
    third_program_name = db_trials(third_program_rows,:).programName(1);
    fourth_program_name = db_trials(fourth_program_rows,:).programName(1);

    if contains(first_program_name, "fine") &&...
        contains(second_program_name, "coarse") &&...
        contains(third_program_name, "fine") &&...
        contains(fourth_program_name, "coarse")
        disp('sanity check passed');
    else
        error('fine and coarse programs not in expected order');
    end
else
    error(strcat("I found ", num2str(totalNumOfPrograms), " programs, but this code can only handle 4 programs (fine, coarse, fine, coarse)"));
end


%% FILTER DATA BY PROGRAM & SANITY CHECK

% ALERT hard-coding fine and coarse label (which is why I added the sanity
% checks above)
fine1_rows = (db_trials.programNum == programNumList(1) &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials);
filtered_trials_fine1 = db_trials(fine1_rows,:);
fine1_total_trials = size(filtered_trials_fine1,1);
if sum(contains(filtered_trials_fine1.programName,"fine"),1) ~= fine1_total_trials
    error('your fine1 label is fishy')
end

coarse1_rows = (db_trials.programNum == programNumList(2) &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials);
filtered_trials_coarse1 = db_trials(coarse1_rows,:);
coarse1_total_trials = size(filtered_trials_coarse1,1);
if sum(contains(filtered_trials_coarse1.programName,"coarse"),1) ~= coarse1_total_trials
    error('your coarse1 label is fishy')
end

fine2_rows = (db_trials.programNum == programNumList(3) &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials);
filtered_trials_fine2 = db_trials(fine2_rows,:);
fine2_total_trials = size(filtered_trials_fine2,1);
if sum(contains(filtered_trials_fine2.programName,"fine"),1) ~= fine2_total_trials
    error('your fine2 label is fishy')
end

coarse2_rows = (db_trials.programNum == programNumList(4) &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials);
filtered_trials_coarse2 = db_trials(coarse2_rows,:);
coarse2_total_trials = size(filtered_trials_coarse2,1);
if sum(contains(filtered_trials_coarse2.programName,"coarse"),1) ~= coarse2_total_trials
    error('your coarse2 label is fishy')
end


%% FIG (lick raster (R vs L) and performance vs program, color-coded by outcome)

% create plot of licks & outcomes in chronological order
figName_prefix = strcat(mouseNum,'_', dateToPlot,...
    '_program_',num2str(first_programNum),'_to_',num2str(last_programNum));
figName_suffix = "_lick rasters_joint";
figure('name',strcat(figName_prefix, figName_suffix));

% 4 rows: fine, coarse, fine, coarse
% 9 columns: L lick raster (2), R lick raster (2), odor ID (1), performance
% (2), outcome by odor (2)
tl = tiledlayout(4,9, 'Padding', 'compact', 'TileSpacing', 'tight');
title(tl,figName_prefix,'Interpreter','none');

% general variables for plot beautification
hide_title = false;
hide_x_labels = true;
hide_y_labels = false;

% plot L licks (spans 1 row and 2 columns)
nexttile(tl,[1,2])
plot_licks("left", filtered_trials_fine1, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);

% plot R licks (spans 1 row and 2 columns)
nexttile(tl,[1,2])
plot_licks("right", filtered_trials_fine1, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, true);

% plot odor ID (spans 1 row and 1 column)
nexttile(tl,[1,1])
plot_odorID_color(filtered_trials_fine1, hit_color, false_choice_color, miss_color, true, hide_x_labels, hide_y_labels)
% plot_odorID(filtered_trials_fine1, odor17_color, odor18_color, odor19_color, true, hide_x_labels, hide_y_labels);

% plot performance (spans 1 row and 2 columns)
nexttile(tl,[1,2])
plot_performance_rotated(filtered_trials_fine1, trials_per_bin, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);

% plot outcome by odor (spans 1 row and 2 columns)
nexttile(tl,[1,2])
plot_performance_bars(filtered_trials_fine1, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels)


% general variables for plot beautification
hide_title = true;
hide_x_labels = true;
hide_y_labels = false;


% plot L licks, R licks, odor ID, and performance for other programs
nexttile(tl,[1,2])
plot_licks("left", filtered_trials_coarse1, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);
nexttile(tl,[1,2])
plot_licks("right", filtered_trials_coarse1, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, true);
nexttile(tl,[1,1])
plot_odorID_color(filtered_trials_coarse1, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);
nexttile(tl,[1,2])
plot_performance_rotated(filtered_trials_coarse1, trials_per_bin, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);
nexttile(tl,[1,2])
plot_performance_bars(filtered_trials_coarse1, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels)

nexttile(tl,[1,2])
plot_licks("left", filtered_trials_fine2, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);
nexttile(tl,[1,2])
plot_licks("right", filtered_trials_fine2, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, true);
nexttile(tl,[1,1])
plot_odorID_color(filtered_trials_fine2, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);
nexttile(tl,[1,2])
plot_performance_rotated(filtered_trials_fine2, trials_per_bin, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);
nexttile(tl,[1,2])
plot_performance_bars(filtered_trials_fine2, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels)

nexttile(tl,[1,2])
plot_licks("left", filtered_trials_coarse2, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, false, hide_y_labels);
nexttile(tl,[1,2])
plot_licks("right", filtered_trials_coarse2, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, false, true);
nexttile(tl,[1,1])
plot_odorID_color(filtered_trials_coarse2, hit_color, false_choice_color, miss_color, hide_title, false, hide_y_labels);
nexttile(tl,[1,2])
plot_performance_rotated(filtered_trials_coarse2, trials_per_bin, hit_color, false_choice_color, miss_color, hide_title, false, hide_y_labels);
nexttile(tl,[1,2])
plot_performance_bars(filtered_trials_coarse2, hit_color, false_choice_color, miss_color, hide_title, false, hide_y_labels)


% adjust final figure position and size
set(gcf,'OuterPosition',[0 100 1000 800]);


%% SAVE PLOTS

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