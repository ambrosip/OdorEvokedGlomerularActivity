%% USER INPUT

dateToPlot = "2025_10_22";
exclude_First_N_Trials = 0;


%% FILTER DATA & PLOT

filtered_rows_by_date = matches(db_trials.date, dateToPlot);
totalNumOfPrograms = size(unique(db_trials(filtered_rows_by_date,:).programNum),1);
for program = unique(db_trials(filtered_rows_by_date,:).programNum)'

    % create filters for plotting licks per trial
    filtered_fc_trials = (db_trials.programNum == program &...
        matches(db_trials.outcome,"false choice") &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials);
    filtered_hit_trials = (db_trials.programNum == program &...
        matches(db_trials.outcome,"hit") &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials);

    % create figure (licks per trial) if you have data to plot
    if ~isempty(db_trials(filtered_hit_trials,:))

        figure('name', strcat(db_trials(filtered_hit_trials,:).mouse(1),...
            '_', dateToPlot,'_program_',num2str(program),'_',...
            db_trials(filtered_hit_trials,:).programName(1),...
            '_licks per trial'));
        hold on;
        
        % plot hits
        plot(db_trials(filtered_hit_trials,:).R_lick_latency_s_from_odor_onset,...
            db_trials(filtered_hit_trials,:).R_lick_total_during_odor_and_response,...
            "o",'Color','b','DisplayName','R lick (Hit)')
        plot(db_trials(filtered_hit_trials,:).L_lick_latency_s_from_odor_onset,...
            db_trials(filtered_hit_trials,:).L_lick_total_during_odor_and_response,...
            "x",'Color','b','DisplayName','L lick (Hit)')
        
        % plot false choices
        plot(db_trials(filtered_fc_trials,:).R_lick_latency_s_from_odor_onset,...
            db_trials(filtered_fc_trials,:).R_lick_total_during_odor_and_response,...
            "o",'Color','r','DisplayName','R lick (FC)')
        plot(db_trials(filtered_fc_trials,:).L_lick_latency_s_from_odor_onset,...
            db_trials(filtered_fc_trials,:).L_lick_total_during_odor_and_response,...
            "x",'Color','r','DisplayName','L lick (FC)')
        
        % beautify plot
        axis([0,4,0,30])
        xlabel('Latency from odor onset (s)');
        ylabel('Total licks')
        title([strcat(db_trials(filtered_hit_trials,:).mouse(1),...
            '_', dateToPlot, '_licks per trial'),...
            strcat('program_',num2str(program), '_',...
            db_trials(filtered_hit_trials,:).programName(1))],...
            'Interpreter','none');
        legend('R lick (Hit)', 'L lick (Hit)', 'R lick (FC)', 'L lick (FC)', 'Location','northeast')
        legend('boxoff')
        hold off;
    end

    % create filters for plotting outcome by odor
    filtered_odor17_trials = (db_trials.programNum == program &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials &...
        db_trials.odorID == 17);
    filtered_odor18_trials = (db_trials.programNum == program &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials &...
        db_trials.odorID == 18);

    % calculate performance and create bar plot
    if ~isempty(db_trials(filtered_odor17_trials,:))

        total_trials_odor17 = size(db_trials(filtered_odor17_trials,:),1);
        total_hit_odor17 = sum(matches(db_trials(filtered_odor17_trials,:).outcome, "hit")) / total_trials_odor17;
        total_fc_odor17 = sum(matches(db_trials(filtered_odor17_trials,:).outcome, "false choice")) / total_trials_odor17;
        total_miss_odor17 = sum(matches(db_trials(filtered_odor17_trials,:).outcome, "miss")) / total_trials_odor17;
    
        total_trials_odor18 = size(db_trials(filtered_odor18_trials,:),1);
        total_hit_odor18 = sum(matches(db_trials(filtered_odor18_trials,:).outcome, "hit")) / total_trials_odor18;
        total_fc_odor18 = sum(matches(db_trials(filtered_odor18_trials,:).outcome, "false choice")) / total_trials_odor18;
        total_miss_odor18 = sum(matches(db_trials(filtered_odor18_trials,:).outcome, "miss")) / total_trials_odor18;
    
        figure('name', strcat(db_trials(filtered_odor17_trials,:).mouse(1),...
            '_', dateToPlot,'_program_',num2str(program),'_',...
            db_trials(filtered_odor17_trials,:).programName(1),...
            '_outcome by odor'));
        hold on;
    
        x = [17 18];
        y = [total_hit_odor17 total_fc_odor17 total_miss_odor17; total_hit_odor18 total_fc_odor18 total_miss_odor18];
    
        bar(x,y,'stacked')
        legend('Hit', 'False choice', 'Miss', 'Location','northeast')
        xlabel('Odors');
        ylabel('Outcomes/Trials')
        title([strcat(db_trials(filtered_odor17_trials,:).mouse(1),...
            '_', dateToPlot, '_outcomes by odor'),...
            strcat('program_',num2str(program), '_',...
            db_trials(filtered_odor17_trials,:).programName(1))],...
            'Interpreter','none');
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