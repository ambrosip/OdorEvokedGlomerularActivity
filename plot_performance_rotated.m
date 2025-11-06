function plot_performance_rotated(filtered_trials, trials_per_bin, hit_color, false_choice_color, miss_color)

% to match the raster, I will make a rotated performance plot, with
% outcome on the x-axis and trials on the y-axis
hold on;
total_trials = size(filtered_trials,1);
trialBins = floor(total_trials/trials_per_bin);
first_trial = 0;
if trialBins >= 2
    for trialBin = 1:trialBins            
        last_trial = trials_per_bin * trialBin;
        first_trial = last_trial - trials_per_bin + 1;
        binned_hit(trialBin) = 100 * sum(matches(filtered_trials(first_trial:last_trial,:).outcome, "hit")) / trials_per_bin;
        binned_fc(trialBin) = 100 * sum(matches(filtered_trials(first_trial:last_trial,:).outcome, "false choice")) / trials_per_bin;
        binned_miss(trialBin) = 100 * sum(matches(filtered_trials(first_trial:last_trial,:).outcome, "miss")) / trials_per_bin;
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