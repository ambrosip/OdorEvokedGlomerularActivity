function plot_odorID_color(filtered_trials, hit_color, false_choice_color, miss_color, hide_title, hide_x_labels, hide_y_labels)

odor17_rows = filtered_trials.odorID == 17;
odor18_rows = filtered_trials.odorID == 18;
odor19_rows = filtered_trials.odorID == 19;

hit_rows = filtered_trials.outcome == "hit";
fc_rows = filtered_trials.outcome == "false choice";
miss_rows = filtered_trials.outcome == "miss";

odor17_hit_rows = odor17_rows .* hit_rows;
odor17_fc_rows = odor17_rows .* fc_rows;
odor17_miss_rows = odor17_rows .* miss_rows;

odor18_hit_rows = odor18_rows .* hit_rows;
odor18_fc_rows = odor18_rows .* fc_rows;
odor18_miss_rows = odor18_rows .* miss_rows;

odor19_hit_rows = odor19_rows .* hit_rows;
odor19_fc_rows = odor19_rows .* fc_rows;
odor19_miss_rows = odor19_rows .* miss_rows;

% assign numbers to each outtome
% 1: hit
% 2: false choice
% 3: miss
odor17_combined_rows = odor17_hit_rows + 2*odor17_fc_rows + 3*odor17_miss_rows;
odor18_combined_rows = odor18_hit_rows + 2*odor18_fc_rows + 3*odor18_miss_rows;
odor19_combined_rows = odor19_hit_rows + 2*odor19_fc_rows + 3*odor19_miss_rows;

all_cols = [odor17_combined_rows, odor18_combined_rows, odor19_combined_rows];
all_cols = double(all_cols);

h = heatmap(all_cols, "Colormap", [[1 1 1]; hit_color; false_choice_color; miss_color],'CellLabelColor','none');
clim([0 3]);

% remove x and y labels
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

% remove colorbar
h.ColorbarVisible = 'off';

% add custom x label
customXLabels = {"α", "α'", "β"};
Ax.XDisplayLabels = customXLabels;

if ~hide_title
    % add title
    title('Odor');
end

% if ~hide_x_labels
%     xlabel('Odor');
% end

% adjust size
% set(gcf,'OuterPosition',[100 100 50 350]);
set(gcf,'Position',[100 100 50 350]);

% adjust font
set(gca,'FontName','Arial');
set(findall(gcf,'-property','FontSize'),'FontSize',12)

end


%% ARCHIVE ALT

% imshow(filtered_trials.odorID, "Colormap", [odor17_color; odor18_color; odor19_color],'DisplayRange', [17,19],'Border', 'tight')
% set(gcf,'OuterPosition',[0 100 300 350]);

% h = heatmap(filtered_trials.odorID, "Colormap", [odor17_color; odor18_color; odor19_color],'CellLabelColor','none');
% set color limits to avoid mis-coloring when not all odors are present
% clim([17 19]);