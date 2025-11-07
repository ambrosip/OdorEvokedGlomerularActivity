function plot_performance_bars(filtered_trials, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels)

% create filters for plotting outcome by odor
filtered_odor17_rows = filtered_trials.odorID == 17;
filtered_odor18_rows = filtered_trials.odorID == 18;
filtered_odor19_rows = filtered_trials.odorID == 19;

% filter by odor
filtered_odor17_trials = filtered_trials(filtered_odor17_rows,:);
filtered_odor18_trials = filtered_trials(filtered_odor18_rows,:);
filtered_odor19_trials = filtered_trials(filtered_odor19_rows,:);

total_trials_odor17 = size(filtered_odor17_trials,1);
total_trials_odor18 = size(filtered_odor18_trials,1);
total_trials_odor19 = size(filtered_odor19_trials,1);

total_hit_odor17 = sum(matches(filtered_odor17_trials.outcome, "hit")) / total_trials_odor17;
total_fc_odor17 = sum(matches(filtered_odor17_trials.outcome, "false choice")) / total_trials_odor17;
total_miss_odor17 = sum(matches(filtered_odor17_trials.outcome, "miss")) / total_trials_odor17;

total_hit_odor18 = sum(matches(filtered_odor18_trials.outcome, "hit")) / total_trials_odor18;
total_fc_odor18 = sum(matches(filtered_odor18_trials.outcome, "false choice")) / total_trials_odor18;
total_miss_odor18 = sum(matches(filtered_odor18_trials.outcome, "miss")) / total_trials_odor18;

total_hit_odor19 = sum(matches(filtered_odor19_trials.outcome, "hit")) / total_trials_odor19;
total_fc_odor19 = sum(matches(filtered_odor19_trials.outcome, "false choice")) / total_trials_odor19;
total_miss_odor19 = sum(matches(filtered_odor19_trials.outcome, "miss")) / total_trials_odor19;

hold on;
% x = [17 18 19];
x = ["α" "α'" "β"];
y = 100 *...
    [total_hit_odor17 total_fc_odor17 total_miss_odor17;...
    total_hit_odor18 total_fc_odor18 total_miss_odor18;...
    total_hit_odor19 total_fc_odor19 total_miss_odor19];

bh = bar(x,y,'stacked');

% outcome color-coding 
bh(1).FaceColor = 'flat';
bh(1).CData = hit_color;        
bh(2).FaceColor = 'flat';
bh(2).CData = false_choice_color;       
bh(3).FaceColor = 'flat';
bh(3).CData = miss_color;

% beautifying
bh_legend = legend('Hit', 'False choice', 'Miss', 'Location','eastoutside');
legend('boxoff')
legend('Direction','normal')
% set(bh_legend, 'Visible','off')

if ~hide_x_labels
    xlabel('Odor');
end

if ~hide_y_labels
    ylabel('Outcome %')
end

ylim([0 100]);
yticks([0,100]);
% xticks([17,18,19]);
xticks = ["α" "α'" "β"];

if ~hide_title
    title('Outcomes by odor');
end

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

end