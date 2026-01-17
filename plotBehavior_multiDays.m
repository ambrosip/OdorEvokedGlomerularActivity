



function plotBehavior_multiDays(datesToPlot, db_trials, saveDir, mouseNum)


%% USER INPUT

% dateToPlot = ["2025_10_24", "2025_10_25"];
exclude_First_N_Trials = 20;

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
filtered_rows_by_date = matches(db_trials.date, datesToPlot);

% figure out if multiple programs were started that day
totalNumOfPrograms = size(unique(db_trials(filtered_rows_by_date,:).programNum),1);

% create empty matrices
hit_per_trials_per_day = [];
false_choices_per_trials_per_day = [];
misses_per_trials_per_day = [];

% iterate through dates
for date = datesToPlot
    
    % iterate through programs and append "2AFC fine" programs
    programNums_to_append = [];
    for program = unique(db_trials(filtered_rows_by_date,:).programNum)'
        program_to_check = db_trials(db_trials.programNum == program,:).programName(1);
        if contains(program_to_check,"2AFC fine")
            programNums_to_append = [programNums_to_append, program];
        end
    end
    
    if size(programNums_to_append,2) > 1
        programNums = num2str(programNums_to_append);
    
        % filter trials for main fig    
        rows_to_analyze_appended = (ismember(db_trials.programNum,programNums_to_append) &...
            matches(db_trials.date, dateToPlot) &...
            db_trials.trialNum_rel > exclude_First_N_Trials);
        trials_to_analyze = db_trials(rows_to_analyze_appended,:);
        trials_to_analyze_total = size(trials_to_analyze,1);
    end
    
        
    
    title_prefix = strcat(trials_to_analyze(1,:).mouse, '_from_',...
        trials_to_analyze(1,:).date(1), '_to_', trials_to_analyze(1,:).date(end));   
    
    % set number of trials to skip from beginning (to discount the first 10
    % pavlovian trials and the following 10 trials)
    trials_to_skip = exclude_First_N_Trials;
    
    % get total number of trials and cap if beyond min number
    total_trials = size(trials_to_analyze,1);
    if total_trials >= 160
        total_trials = 140;
    end
    
    % adjust total_trials to account for skipped trials from beginning
    total_trials = total_trials - trials_to_skip;
    
    avg_trial_dur_s = trials_to_analyze(1,:).avg_trial_dur_s;
    total_dur_plotted_min = round((total_trials * avg_trial_dur_s)/60);
    
    hold on;
    trialBins = floor(total_trials/trials_per_bin);
    if trialBins >= 2
        for trialBin = 1:trialBins            
            last_trial = trials_to_skip + trials_per_bin * trialBin;
            first_trial = last_trial - trials_per_bin + 1;
            binned_hit(trialBin) = 100 * sum(matches(trials_to_analyze(first_trial:last_trial,:).outcome, "hit")) / trials_per_bin;
            binned_fc(trialBin) = 100 * sum(matches(trials_to_analyze(first_trial:last_trial,:).outcome, "false choice")) / trials_per_bin;
            binned_miss(trialBin) = 100 * sum(matches(trials_to_analyze(first_trial:last_trial,:).outcome, "miss")) / trials_per_bin;
        end       
        
        plot(binned_hit,'Color',hit_color,'LineWidth',1,'DisplayName','Hit');
        plot(binned_fc,'Color',false_choice_color,'LineWidth',1,'DisplayName','False choice');
        plot(binned_miss,'Color',miss_color,'LineWidth',1,'DisplayName','Miss');
    
        % beautifying
        % note the rotation of x and y axis
        yline(50,'--')
        axis([1 trialBins 0 100])
        yticks([0,100]);
        xticks([1,trialBins]);
    
        if ~hide_x_labels
            xlabel(strcat("Session (", num2str(trials_per_bin), " trials)"));
        end
    
        if ~hide_x_labels
            ylabel('Outcome (%)');
        end
    
        full_title = [title_prefix,...     
            trials_to_analyze(1,:).programName,...
            strcat("Total time ~", num2str(total_dur_plotted_min), " min"),...
            ""];
    
        if ~hide_title
            title(full_title, 'Interpreter','none');       
        end
    
        legend('Hit', 'False choice', 'Miss','Location','eastoutside')
        legend('boxoff')
        % set(gca, 'YDir', 'reverse'); 
        set(gca,'FontName','Arial');
        set(gca,'TickLength', [0.025, 0.025]);
        set(gca,'LineWidth', 0.75);
        set(findall(gcf,'-property','FontSize'),'FontSize',12)
        set(gcf,'OuterPosition',[0 100 350 350]);
        hold off; 
    
    end
    end