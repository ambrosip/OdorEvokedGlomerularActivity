function plot_licks(lick_side, filtered_trials, odorDur_s, totalDur_s, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels)

total_trials = size(filtered_trials,1);
for trial = 1:total_trials
    if filtered_trials.outcome(trial) == "hit"
        color_to_use = hit_color;
    elseif filtered_trials.outcome(trial) == "false choice"
        color_to_use = false_choice_color;
    elseif filtered_trials.outcome(trial) == "miss"
        color_to_use = miss_color;
    end

    if lick_side == "left"
        plot(filtered_trials.all_lick_L_latency_per_trial{trial}, ...
            trial, '.', 'LineWidth', 1, 'Color', color_to_use)    
        hold on;
    elseif lick_side == "right"
        plot(filtered_trials.all_lick_R_latency_per_trial{trial}, ...
            trial, '.', 'LineWidth', 1, 'Color', color_to_use)
        hold on;
    end
end

% beautifying
if ~hide_title
    if lick_side == "left"
        title('L licks');
    elseif lick_side == "right"
        title('R licks');
    end
end
axis([0,totalDur_s,0,total_trials+1])  

if ~hide_x_labels
    xlabel('Time from odor onset (s)');
end

if ~hide_y_labels
    if contains(filtered_trials(1,:).programName, "fine")
        ylabel('Trials (Fine)')
    elseif contains(filtered_trials(1,:).programName, "coarse")
        ylabel('Trials (Coarse)')
    else
        ylabel('Trials')
    end
end
yticks([1,total_trials]);
xticks([0,odorDur_s,totalDur_s]); 
% set(gca,'TickDir','out');
set(gca,'FontName','Arial');
set(gca,'TickLength', [0.025, 0.025]);
set(gca,'LineWidth', 0.75);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf,'OuterPosition',[0 100 300 350]);
set(gca,'Ydir','reverse');
hold off;
box off;

end