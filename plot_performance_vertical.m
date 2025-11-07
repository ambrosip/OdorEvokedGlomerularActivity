function plot_performance_vertical(filtered_trials, trials_per_bin, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels)

title_prefix = strcat(filtered_trials(1,:).mouse, '_',...
    filtered_trials(1,:).date, '_program_',...
    num2str(filtered_trials(1,:).programNum));   

% if standard trials_per_bin = 5, set trials_per_bin to 10.
trials_per_bin = 2 * trials_per_bin;

% set number of trials to skip from beginning (to discount the first 10
% pavlovian trials and the following 10 trials)
trials_to_skip = 20;

% get total number of trials and cap if beyond min number
total_trials = size(filtered_trials,1);
if total_trials >= 160
    total_trials = 140;
end

% adjust total_trials to account for skipped trials from beginning
total_trials = total_trials - trials_to_skip;

avg_trial_dur_s = filtered_trials(1,:).avg_trial_dur_s;
total_dur_plotted_min = round((total_trials * avg_trial_dur_s)/60);

hold on;
trialBins = floor(total_trials/trials_per_bin);
if trialBins >= 2
    for trialBin = 1:trialBins            
        last_trial = trials_to_skip + trials_per_bin * trialBin;
        first_trial = last_trial - trials_per_bin + 1;
        binned_hit(trialBin) = 100 * sum(matches(filtered_trials(first_trial:last_trial,:).outcome, "hit")) / trials_per_bin;
        binned_fc(trialBin) = 100 * sum(matches(filtered_trials(first_trial:last_trial,:).outcome, "false choice")) / trials_per_bin;
        binned_miss(trialBin) = 100 * sum(matches(filtered_trials(first_trial:last_trial,:).outcome, "miss")) / trials_per_bin;
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
        filtered_trials(1,:).programName,...
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