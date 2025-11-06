%{ 
DOCUMENTATION

run analyzeBehavior before using this script

%}

% function plotBehavior(date, db_trials, saveDir)

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

    % clearing things that seem to linger from one run to the next
    clear x_bins;
    clear binned_hit_odor17;
    clear binned_fc_odor17;
    clear binned_miss_odor17;
    clear binned_hit_odor18;
    clear binned_fc_odor18;
    clear binned_miss_odor18;
    clear binned_hit;
    clear binned_fc;
    clear binned_miss;

    % clearing other things just to be safe
    clear total_trials_odor17; 
    clear total_hit_odor17; 
    clear total_fc_odor17; 
    clear total_miss_odor17; 
    clear total_trials_odor18; 
    clear total_hit_odor18; 
    clear total_fc_odor18; 
    clear total_miss_odor18; 
    clear x;
    clear y;
    clear trialBins;
    clear filtered_data;
    clear first_trial;
    clear last_trial;
    
%%  FIG (lick latency vs lick number, color-coded by outcome)

    % create filters for plotting licks per trial
    filtered_fc_trials = (db_trials.programNum == program &...
        matches(db_trials.outcome,"false choice") &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials);
    filtered_hit_trials = (db_trials.programNum == program &...
        matches(db_trials.outcome,"hit") &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials);

    % create licks per trial plot if you have data to plot
    if ~isempty(db_trials(filtered_hit_trials,:))

        figure('name', strcat(db_trials(filtered_hit_trials,:).mouse(1),...
            '_', dateToPlot,'_program_',num2str(program),'_',...
            db_trials(filtered_hit_trials,:).programName(1),...
            '_licks per trial'));
        hold on;
        
        % plot hits
        plot(db_trials(filtered_hit_trials,:).first_R_lick_latency_s_from_odor_onset,...
            db_trials(filtered_hit_trials,:).R_lick_total_during_odor_and_response / totalDur_s,...
            "o",'Color',hit_color,'DisplayName','R (Hit)')
        plot(db_trials(filtered_hit_trials,:).first_L_lick_latency_s_from_odor_onset,...
            db_trials(filtered_hit_trials,:).L_lick_total_during_odor_and_response / totalDur_s,...
            "x",'Color',hit_color,'DisplayName','L (Hit)')
               
        % plot false choices
        plot(db_trials(filtered_fc_trials,:).first_R_lick_latency_s_from_odor_onset,...
            db_trials(filtered_fc_trials,:).R_lick_total_during_odor_and_response / totalDur_s,...
            "o",'Color',false_choice_color,'DisplayName','R (False choice)')
        plot(db_trials(filtered_fc_trials,:).first_L_lick_latency_s_from_odor_onset,...
            db_trials(filtered_fc_trials,:).L_lick_total_during_odor_and_response / totalDur_s,...
            "x",'Color',false_choice_color,'DisplayName','L (False choice)')
        
        % beautify plot
        axis([0,totalDur_s,0,7])        
        xlabel('First lick latency from odor onset (s)');
        ylabel('Lick rate (Hz)')
        yticks([0,7]);
        xticks([0,odorDur_s,totalDur_s]);
        title([strcat(db_trials(filtered_hit_trials,:).mouse(1),...
            '_', dateToPlot),...
            'licks per trial',...
            strcat('program_',num2str(program)),...
            db_trials(filtered_hit_trials,:).programName(1),...
            ''],...
            'Interpreter','none');
        legend('R (Hit)', 'L (Hit)', 'R (FC)', 'L (FC)', 'Location','eastoutside') % Alt location: northeast
        legend('boxoff')
        set(gca,'FontName','Arial');
        set(gca,'TickLength', [0.025, 0.025]);
        set(gca,'LineWidth', 0.75);
        set(findall(gcf,'-property','FontSize'),'FontSize',12)
        set(gcf,'OuterPosition',[0 100 300 350]);
        hold off;
    else
        disp('bruh you dont got no licks to plot');
    end

    
%%  FIG (outcomes vs odors, stacked bars)

    % TO DO: add if statement for coarse program and plot odor 19

    % create filters for plotting outcome by odor
    filtered_odor17_trials = (db_trials.programNum == program &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials &...
        db_trials.odorID == 17);
    filtered_odor18_trials = (db_trials.programNum == program &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials &...
        db_trials.odorID == 18);

    % calculate outcomes by odor (if you have data)
    if size(db_trials(filtered_odor17_trials,:),1)>=2*trials_per_bin & size(db_trials(filtered_odor18_trials,:),1)>=2*trials_per_bin

        total_trials_odor17 = size(db_trials(filtered_odor17_trials,:),1);
        total_hit_odor17 = sum(matches(db_trials(filtered_odor17_trials,:).outcome, "hit")) / total_trials_odor17;
        total_fc_odor17 = sum(matches(db_trials(filtered_odor17_trials,:).outcome, "false choice")) / total_trials_odor17;
        total_miss_odor17 = sum(matches(db_trials(filtered_odor17_trials,:).outcome, "miss")) / total_trials_odor17;
    
        total_trials_odor18 = size(db_trials(filtered_odor18_trials,:),1);
        total_hit_odor18 = sum(matches(db_trials(filtered_odor18_trials,:).outcome, "hit")) / total_trials_odor18;
        total_fc_odor18 = sum(matches(db_trials(filtered_odor18_trials,:).outcome, "false choice")) / total_trials_odor18;
        total_miss_odor18 = sum(matches(db_trials(filtered_odor18_trials,:).outcome, "miss")) / total_trials_odor18;
   
        % plot bar graph of outcomes by odor
        figure('name', strcat(db_trials(filtered_odor17_trials,:).mouse(1),...
            '_', dateToPlot,'_program_',num2str(program),'_',...
            db_trials(filtered_odor17_trials,:).programName(1),...
            '_outcome by odor'));
        hold on;
    
        x = [17 18];
        y = 100*[total_hit_odor17 total_fc_odor17 total_miss_odor17; total_hit_odor18 total_fc_odor18 total_miss_odor18];
    
        bh = bar(x,y,'stacked');

        % outcome color-coding 
        bh(1).FaceColor = 'flat';
        bh(1).CData = hit_color;        
        bh(2).FaceColor = 'flat';
        bh(2).CData = false_choice_color;       
        bh(3).FaceColor = 'flat';
        bh(3).CData = miss_color;

        % beautifying
        legend('Hit', 'False choice', 'Miss', 'Location','eastoutside')
        legend('boxoff')
        xlabel('Odor');
        ylabel('Outcome %')
        yticks([0,100]);
        xticks([17,18]);
        title([strcat(db_trials(filtered_odor17_trials,:).mouse(1),...
            '_', dateToPlot),...
            'outcomes by odor',...
            strcat('program_',num2str(program)),...
            db_trials(filtered_odor17_trials,:).programName(1),...
            ''],...
            'Interpreter','none');
        for bar_plotted = 1:3
            % bh(bar_plotted).LineStyle = 'none';
            bh(bar_plotted).BarWidth = 0.5;
        end
        set(gca,'FontName','Arial');
        set(gca,'TickLength', [0.025, 0.025]);
        set(gca,'LineWidth', 0.75);
        set(findall(gcf,'-property','FontSize'),'FontSize',12)
        set(gcf,'OuterPosition',[0 100 300 350]); 
        hold off;


%% FIG (lick raster (R vs L) and performance vs odor, color-coded by outcome)

        % create plot of licks & outcomes by odor in chronological order
        % iterate through relative trials
        figure('name', strcat(db_trials(filtered_odor17_trials,:).mouse(1),...
            '_', dateToPlot,'_program_',num2str(program),'_',...
            db_trials(filtered_odor17_trials,:).programName(1),...
            '_lick rasters'));

        % 3 columns: L lick raster, R lick raster, performance
        % 2 rows: odor 17, odor 18 
        tl = tiledlayout(2,3);
        title(tl, strcat(db_trials(filtered_odor17_trials,:).mouse(1),...
            '_', dateToPlot, '_program_',num2str(program),'_',...
            db_trials(filtered_odor17_trials,:).programName(1)),...
            'Interpreter','none');

        % plot L licks for odor 17
        nexttile
        hold on;
        for trial = 1:total_trials_odor17
            if db_trials(filtered_odor17_trials,:).outcome(trial) == "hit"
                color_to_use = hit_color;
            elseif db_trials(filtered_odor17_trials,:).outcome(trial) == "false choice"
                color_to_use = false_choice_color;
            elseif db_trials(filtered_odor17_trials,:).outcome(trial) == "miss"
                color_to_use = miss_color;
            end
            plot(db_trials(filtered_odor17_trials,:).all_lick_L_latency_per_trial{trial}, ...
                trial, '.', 'LineWidth', 1, 'Color', color_to_use)
        end

            % beautifying
            axis([0,totalDur_s,0,total_trials_odor17+1])  
            title('L licks');
            xlabel('Time from odor onset (s)');
            ylabel('Trials (Odor 17)')
            yticks([1,total_trials_odor17]);
            xticks([0,odorDur_s,totalDur_s]);            
            set(gca,'FontName','Arial');
            set(gca,'TickLength', [0.025, 0.025]);
            set(gca,'LineWidth', 0.75);
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
            set(gcf,'OuterPosition',[0 100 300 350]);
            set(gca,'Ydir','reverse');
            hold off;

        % plot R licks for odor 17
        nexttile
        hold on;
        for trial = 1:total_trials_odor17
            if db_trials(filtered_odor17_trials,:).outcome(trial) == "hit"
                color_to_use = hit_color;
            elseif db_trials(filtered_odor17_trials,:).outcome(trial) == "false choice"
                color_to_use = false_choice_color;
            elseif db_trials(filtered_odor17_trials,:).outcome(trial) == "miss"
                color_to_use = miss_color;
            end
            plot(db_trials(filtered_odor17_trials,:).all_lick_R_latency_per_trial{trial}, ...
                trial, '.', 'LineWidth', 1, 'Color', color_to_use)
        end

            % beautifying
            axis([0,totalDur_s,0,total_trials_odor17+1])  
            title('R licks');
            xlabel('Time from odor onset (s)');
            ylabel('Trials (Odor 17)')
            yticks([1,total_trials_odor17]);
            xticks([0,odorDur_s,totalDur_s]);            
            set(gca,'FontName','Arial');
            set(gca,'TickLength', [0.025, 0.025]);
            set(gca,'LineWidth', 0.75);
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
            set(gcf,'OuterPosition',[0 100 300 350]);
            set(gca,'Ydir','reverse');
            hold off;       

        % plot performance for odor 17
        % to match the raster, I will make a rotated performance plot, with
        % outcome on the x-axis and trials on the y-axis
        nexttile
        hold on;
        % can use round to discard trials in the last bin if we have less
        % trials than 1/2 of trials_per_bin, but I will instead just
        % throw away the last trials if they don't fit neatily in a bin (ie
        % I will use "floor" instead of "round").
        % trialBins = round(total_trials_odor17/trials_per_bin);
        trialBins = floor(total_trials_odor17/trials_per_bin);
        filtered_data = db_trials(filtered_odor17_trials,:); 
        first_trial = 0;
        if trialBins >= 2
            for trialBin = 1:trialBins            
                last_trial = trials_per_bin * trialBin;
                first_trial = last_trial - trials_per_bin + 1;
                binned_hit_odor17(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "hit")) / trials_per_bin;
                binned_fc_odor17(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "false choice")) / trials_per_bin;
                binned_miss_odor17(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "miss")) / trials_per_bin;
            end       
            % create a vector of trialBins that will be used to make a rotated
            % outcome x session plot
            x_bins = 1:trialBins;        
            plot(binned_hit_odor17,x_bins,'Color',hit_color,'LineWidth',1,'DisplayName','Hit');
            plot(binned_fc_odor17,x_bins,'Color',false_choice_color,'LineWidth',1,'DisplayName','False choice');
            plot(binned_miss_odor17,x_bins,'Color',miss_color,'LineWidth',1,'DisplayName','Miss');
            % FYI if you want to add markers to this plot, do this
            % plot(binned_hit_odor17,x_bins,'-o','Color',hit_color,'LineWidth',1,'DisplayName','Hit');
                
                % beautifying
                % note the rotation of x and y axis
                xline(50,'--')
                axis([0 100 1 trialBins])
                xticks([0,100]);
                yticks([1,trialBins]);
                ylabel(strcat("Session (", num2str(trials_per_bin), " trials)"));
                xlabel('Outcome (%)');
                title('Performance');
                legend('Hit', 'False choice', 'Miss','Location','eastoutside')
                legend('boxoff')
                set(gca, 'YDir', 'reverse'); 
                set(gca,'FontName','Arial');
                set(gca,'TickLength', [0.025, 0.025]);
                set(gca,'LineWidth', 0.75);
                set(findall(gcf,'-property','FontSize'),'FontSize',12)
                set(gcf,'OuterPosition',[0 100 750 600]);
                hold off; 
        end

        % plot L licks for odor 18
        nexttile
        hold on;
        for trial = 1:total_trials_odor18
            if db_trials(filtered_odor18_trials,:).outcome(trial) == "hit"
                color_to_use = hit_color;
            elseif db_trials(filtered_odor18_trials,:).outcome(trial) == "false choice"
                color_to_use = false_choice_color;
            elseif db_trials(filtered_odor18_trials,:).outcome(trial) == "miss"
                color_to_use = miss_color;
            end
            plot(db_trials(filtered_odor18_trials,:).all_lick_L_latency_per_trial{trial}, ...
                trial, '.', 'LineWidth', 1, 'Color', color_to_use)
        end

            % beautifying
            axis([0,totalDur_s,0,total_trials_odor18+1])  
            title('L licks');
            xlabel('Time from odor onset (s)');
            ylabel('Trials (Odor 18)')
            yticks([1,total_trials_odor18]);
            xticks([0,odorDur_s,totalDur_s]);            
            set(gca,'FontName','Arial');
            set(gca,'TickLength', [0.025, 0.025]);
            set(gca,'LineWidth', 0.75);
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
            set(gcf,'OuterPosition',[0 100 300 350]);
            set(gca,'Ydir','reverse');
            hold off;

        % plot R licks for odor 18
        nexttile
        hold on;
        for trial = 1:total_trials_odor18
            if db_trials(filtered_odor18_trials,:).outcome(trial) == "hit"
                color_to_use = hit_color;
            elseif db_trials(filtered_odor18_trials,:).outcome(trial) == "false choice"
                color_to_use = false_choice_color;
            elseif db_trials(filtered_odor18_trials,:).outcome(trial) == "miss"
                color_to_use = miss_color;
            end
            plot(db_trials(filtered_odor18_trials,:).all_lick_R_latency_per_trial{trial}, ...
                trial, '.', 'LineWidth', 1, 'Color', color_to_use)
        end

            % beautifying
            axis([0,totalDur_s,0,total_trials_odor18+1])  
            title('R licks');
            xlabel('Time from odor onset (s)');
            ylabel('Trials (Odor 18)')
            yticks([1,total_trials_odor18]);
            xticks([0,odorDur_s,totalDur_s]);            
            set(gca,'FontName','Arial');
            set(gca,'TickLength', [0.025, 0.025]);
            set(gca,'LineWidth', 0.75);
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
            set(gcf,'OuterPosition',[0 100 500 600]);
            set(gca,'Ydir','reverse');
            hold off;    

        % plot performance for odor 18
        nexttile
        hold on;
        trialBins = floor(total_trials_odor18/trials_per_bin);
        filtered_data = db_trials(filtered_odor18_trials,:); 
        first_trial = 0;
        if trialBins >= 2
            for trialBin = 1:trialBins            
                last_trial = trials_per_bin * trialBin;
                first_trial = last_trial - trials_per_bin + 1;
                binned_hit_odor18(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "hit")) / trials_per_bin;
                binned_fc_odor18(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "false choice")) / trials_per_bin;
                binned_miss_odor18(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "miss")) / trials_per_bin;
            end       
            x_bins = 1:trialBins;        
            plot(binned_hit_odor18,x_bins,'Color',hit_color,'LineWidth',1,'DisplayName','Hit');
            plot(binned_fc_odor18,x_bins,'Color',false_choice_color,'LineWidth',1,'DisplayName','False choice');
            plot(binned_miss_odor18,x_bins,'Color',miss_color,'LineWidth',1,'DisplayName','Miss');
                
                % beautifying
                xline(50,'--')
                axis([0 100 1 trialBins])
                xticks([0,100]);
                yticks([1,trialBins]);
                ylabel(strcat("Session (", num2str(trials_per_bin), " trials)"));
                xlabel('Outcome (%)');
                title('Performance');
                legend('Hit', 'False choice', 'Miss','Location','eastoutside')
                legend('boxoff')
                set(gca, 'YDir', 'reverse'); 
                set(gca,'FontName','Arial');
                set(gca,'TickLength', [0.025, 0.025]);
                set(gca,'LineWidth', 0.75);
                set(findall(gcf,'-property','FontSize'),'FontSize',12)
                set(gcf,'OuterPosition',[0 100 750 600]);
                hold off; 
        end
    end

%%  FIG

    % create filters for plotting lick raster and performance
    % chronologically, regardless of odor
    filtered_trials = (db_trials.programNum == program &...
        matches(db_trials.date, dateToPlot) &...
        db_trials.trialNum_rel > exclude_First_N_Trials);

    % create plot of licks & outcomes in chronological order
    % iterate through relative trials
    figure('name', strcat(db_trials(filtered_trials,:).mouse(1),...
        '_', dateToPlot,'_program_',num2str(program),'_',...
        db_trials(filtered_trials,:).programName(1),...
        '_lick rasters_joint odors'));

    % 3 columns: L lick raster, R lick raster, performance
    % 1 rows: all odors
    tl = tiledlayout(1,3);
    title(tl, strcat(db_trials(filtered_trials,:).mouse(1),...
        '_', dateToPlot, '_program_',num2str(program),'_',...
        db_trials(filtered_trials,:).programName(1)),...
        'Interpreter','none');

    % plot L licks
    nexttile
    hold on;
    total_trials = size(db_trials(filtered_trials,:),1);
    for trial = 1:total_trials
        if db_trials(filtered_trials,:).outcome(trial) == "hit"
            color_to_use = hit_color;
        elseif db_trials(filtered_trials,:).outcome(trial) == "false choice"
            color_to_use = false_choice_color;
        elseif db_trials(filtered_trials,:).outcome(trial) == "miss"
            color_to_use = miss_color;
        end
        plot(db_trials(filtered_trials,:).all_lick_L_latency_per_trial{trial}, ...
            trial, '.', 'LineWidth', 1, 'Color', color_to_use)
    end

        % beautifying
        axis([0,totalDur_s,0,total_trials+1])  
        title('L licks');
        xlabel('Time from odor onset (s)');
        ylabel('Trials')
        yticks([1,total_trials]);
        xticks([0,odorDur_s,totalDur_s]);            
        set(gca,'FontName','Arial');
        set(gca,'TickLength', [0.025, 0.025]);
        set(gca,'LineWidth', 0.75);
        set(findall(gcf,'-property','FontSize'),'FontSize',12)
        set(gcf,'OuterPosition',[0 100 300 350]);
        set(gca,'Ydir','reverse');
        hold off;

    % plot R licks
    nexttile
    hold on;
    total_trials = size(db_trials(filtered_trials,:),1);
    for trial = 1:total_trials
        if db_trials(filtered_trials,:).outcome(trial) == "hit"
            color_to_use = hit_color;
        elseif db_trials(filtered_trials,:).outcome(trial) == "false choice"
            color_to_use = false_choice_color;
        elseif db_trials(filtered_trials,:).outcome(trial) == "miss"
            color_to_use = miss_color;
        end
        plot(db_trials(filtered_trials,:).all_lick_R_latency_per_trial{trial}, ...
            trial, '.', 'LineWidth', 1, 'Color', color_to_use)
    end

        % beautifying
        axis([0,totalDur_s,0,total_trials+1])  
        title('R licks');
        xlabel('Time from odor onset (s)');
        ylabel('Trials')
        yticks([1,total_trials]);
        xticks([0,odorDur_s,totalDur_s]);            
        set(gca,'FontName','Arial');
        set(gca,'TickLength', [0.025, 0.025]);
        set(gca,'LineWidth', 0.75);
        set(findall(gcf,'-property','FontSize'),'FontSize',12)
        set(gcf,'OuterPosition',[0 100 300 350]);
        set(gca,'Ydir','reverse');
        hold off;   

    % plot performance 
    % to match the raster, I will make a rotated performance plot, with
    % outcome on the x-axis and trials on the y-axis
    nexttile
    hold on;
    % can use round to discard trials in the last bin if we have less
    % trials than 1/2 of trials_per_bin, but I will instead just
    % throw away the last trials if they don't fit neatily in a bin (ie
    % I will use "floor" instead of "round").
    % trialBins = round(total_trials_odor17/trials_per_bin);
    trialBins = floor(total_trials/trials_per_bin);
    filtered_data = db_trials(filtered_trials,:); 
    first_trial = 0;
    if trialBins >= 2
        for trialBin = 1:trialBins            
            last_trial = trials_per_bin * trialBin;
            first_trial = last_trial - trials_per_bin + 1;
            binned_hit(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "hit")) / trials_per_bin;
            binned_fc(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "false choice")) / trials_per_bin;
            binned_miss(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "miss")) / trials_per_bin;
        end       
        % create a vector of trialBins that will be used to make a rotated
        % outcome x session plot
        x_bins = 1:trialBins;        
        plot(binned_hit,x_bins,'Color',hit_color,'LineWidth',1,'DisplayName','Hit');
        plot(binned_fc,x_bins,'Color',false_choice_color,'LineWidth',1,'DisplayName','False choice');
        plot(binned_miss,x_bins,'Color',miss_color,'LineWidth',1,'DisplayName','Miss');
            
            % beautifying
            % note the rotation of x and y axis
            xline(50,'--')
            axis([0 100 1 trialBins])
            xticks([0,100]);
            yticks([1,trialBins]);
            ylabel(strcat("Session (", num2str(trials_per_bin), " trials)"));
            xlabel('Outcome (%)');
            title('Performance');
            legend('Hit', 'False choice', 'Miss','Location','eastoutside')
            legend('boxoff')
            set(gca, 'YDir', 'reverse'); 
            set(gca,'FontName','Arial');
            set(gca,'TickLength', [0.025, 0.025]);
            set(gca,'LineWidth', 0.75);
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
            set(gcf,'OuterPosition',[0 100 750 600]);
            hold off; 
    end
end


%% SAVE PLOTS

% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% 
% % save all open figs
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName = FigList(iFig).Name;
%   set(0, 'CurrentFigure', FigHandle);
%   % forces matlab to save fig as a vector
%   FigHandle.Renderer = 'painters';  
%   % actually saves a vector file
%   saveas(FigHandle,fullfile(saveDir, [FigName '.svg']));
% end 
% disp('saved all figs')
% close all

% end


%% ARCHIVE

% previous version of performance code, which puts outcome in the y-axis
% and sessions in the x-axis
        % % plot performance for odor 17
        % nexttile
        % hold on;
        % % can use round to discard trials in the last bin if we have less
        % % trials than 1/2 of trials_per_bin, but I will instead just
        % % throw away the last trials if they don't fit neatily in a bin (ie
        % % I will use "floor" instead of "round").
        % % trialBins = round(total_trials_odor17/trials_per_bin);
        % trialBins = floor(total_trials_odor17/trials_per_bin);
        % filtered_data = db_trials(filtered_odor17_trials,:); 
        % for trialBin = 1:trialBins            
        %     last_trial = trials_per_bin * trialBin;
        %     first_trial = last_trial - trials_per_bin + 1;
        %     binned_hit_odor17(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "hit")) / trials_per_bin;
        %     binned_fc_odor17(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "false choice")) / trials_per_bin;
        %     binned_miss_odor17(trialBin) = 100 * sum(matches(filtered_data(first_trial:last_trial,:).outcome, "miss")) / trials_per_bin;
        % end        
        % plot(binned_hit_odor17,'Color',hit_color,'LineWidth',1,'DisplayName','Hit');
        % plot(binned_fc_odor17,'Color',false_choice_color,'LineWidth',1,'DisplayName','False choice');
        % plot(binned_miss_odor17,'Color',miss_color,'LineWidth',1,'DisplayName','Miss');
        % 
        %     % beautifying
        %     yline(50,'--')
        %     axis([1 trialBins 0 100])
        %     yticks([0,100]);
        %     xticks([1,trialBins]);
        %     xlabel('Session (10 trials)');
        %     ylabel('Outcome (%)');
        %     title('Performance');
        %     legend('Hit', 'False choice', 'Miss','Location','eastoutside')
        %     legend('boxoff')
        %     set(gca,'FontName','Arial');
        %     set(gca,'TickLength', [0.025, 0.025]);
        %     set(gca,'LineWidth', 0.75);
        %     set(findall(gcf,'-property','FontSize'),'FontSize',12)
        %     set(gcf,'OuterPosition',[0 100 300 350]);
        %     hold off;  