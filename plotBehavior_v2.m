%{ 
DOCUMENTATION

run analyzeBehavior before using this script

%}

function plotBehavior_v2(date, db_trials, saveDir, mouseNum)

close all;

%% USER INPUT

% dateToPlot = "2025_10_24";
dateToPlot = date;
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


%% FILTER DATA & PLOT

% filter data by date based on user input
filtered_rows_by_date = matches(db_trials.date, dateToPlot);

% figure out if multiple programs were started that day
totalNumOfPrograms = size(unique(db_trials(filtered_rows_by_date,:).programNum),1);

% iterate through programs
for program = unique(db_trials(filtered_rows_by_date,:).programNum)'

    % filter trials for main fig
    rows_to_analyze = (db_trials.programNum == program &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials);
    trials_to_analyze = db_trials(rows_to_analyze,:);
    trials_to_analyze_total = size(trials_to_analyze,1);

    % create filters for plotting outcome by odor
    odor17_rows = (db_trials.programNum == program &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials &...
        db_trials.odorID == 17);
    odor18_rows = (db_trials.programNum == program &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials &...
        db_trials.odorID == 18);
    trials_odor17 = db_trials(odor17_rows,:);
    trials_odor18 = db_trials(odor18_rows,:);
    trials_odor17_total = size(trials_odor17,1);
    trials_odor18_total = size(trials_odor18,1);

    % calculate outcomes by odor (if you have data)
    if trials_odor17_total >= 2 * trials_per_bin &&...
            trials_odor18_total >= 2 * trials_per_bin

        % create plot of licks & outcomes in chronological order
        figName_prefix = strcat(mouseNum,'_', dateToPlot,...
            '_program_',num2str(program),'_',...
            trials_to_analyze.programName(1));
        figName_suffix = "_all";

        % create performance plot
        figure('name',strcat(figName_prefix, "_vertical"));
        plot_performance_vertical(trials_to_analyze, trials_per_bin, hit_color, false_choice_color, miss_color, false, false, false)
      
        % create main figure
        figure('name',strcat(figName_prefix, figName_suffix));
        
        % 4 rows: chronological (2), odor 17 (1), odor 18 (1)
        % 9 columns: L lick raster (2), R lick raster (2), odor ID (1), performance
        % (2), outcome by odor (2)
        tl = tiledlayout(4,9, 'Padding', 'compact', 'TileSpacing', 'tight');
        title(tl,figName_prefix,'Interpreter','none');
        
        % general variables for plot beautification
        hide_title = false;
        hide_x_labels = false;
        hide_y_labels = false;

        % plot L licks (spans 2 rows and 2 columns)
        nexttile([2,2])
        plot_licks("left", trials_to_analyze, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);
        
        % plot R licks (spans 2 rows and 2 columns)
        nexttile([2,2])
        plot_licks("right", trials_to_analyze, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, true);
        
        % plot odor ID (spans 2 rows and 1 column)
        nexttile([2,1])
        plot_odorID_color(trials_to_analyze, hit_color, false_choice_color, miss_color, true, hide_x_labels, hide_y_labels)
        
        % plot performance (spans 2 rows and 2 columns)
        nexttile([2,2])
        plot_performance_rotated(trials_to_analyze, trials_per_bin, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);
        
        % plot outcome by odor (spans 2 rows and 2 columns)
        nexttile([2,2])
        plot_performance_bars(trials_to_analyze, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels)
        
        % % empty slot
        % nexttile(tl,[1,2])

        % general variables for plot beautification
        hide_title = false;
        hide_x_labels = true;
        hide_y_labels = false;

        % plot lick raster, odor id, and performance for odor 17
        nexttile(19,[1,2])
        plot_licks("left", trials_odor17, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);
        nexttile(21,[1,2])
        plot_licks("right", trials_odor17, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, true);
        nexttile(23,[1,1])
        plot_odorID_color(trials_odor17, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);
        nexttile(24,[1,2])
        plot_performance_rotated(trials_odor17, trials_per_bin, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels);
        
        % % empty
        % nexttile(tl,[1,2])
        hide_title = true;

        % plot lick raster, odor id, and performance for odor 18
        nexttile(28,[1,2])
        plot_licks("left", trials_odor18, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, false, hide_y_labels);
        nexttile(30,[1,2])
        plot_licks("right", trials_odor18, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, false, true);
        nexttile(32,[1,1])
        plot_odorID_color(trials_odor18, hit_color, false_choice_color, miss_color, hide_title, false, hide_y_labels);
        nexttile(33,[1,2])
        plot_performance_rotated(trials_odor18, trials_per_bin, hit_color, false_choice_color, miss_color, hide_title, false, hide_y_labels);
        
        % % empty
        % nexttile(tl,[1,2])

        % adjust final figure position and size
        set(gcf,'OuterPosition',[0 100 1000 800]);
    end
end


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
