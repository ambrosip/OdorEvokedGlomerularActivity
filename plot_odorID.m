function plot_odorID(filtered_trials, odor17_color, odor18_color, odor19_color)

col1 = filtered_trials.odorID == 17;
col2 = filtered_trials.odorID == 18;
col3 = filtered_trials.odorID == 19;

all_cols = [col1, col2, col3];
all_cols = double(all_cols);

h = heatmap(all_cols, "Colormap", [[1 1 1];[0 0 0]],'CellLabelColor','none');

% remove x and y labels
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

% add custom x label
customXLabels = {'A', 'B', 'C'};
Ax.XDisplayLabels = customXLabels;

% remove colorbar
h.ColorbarVisible = 'off';

% adjust size
% set(gcf,'OuterPosition',[100 100 50 350]);
set(gcf,'Position',[100 100 50 350]);

title('Odor');
set(gca,'FontName','Arial');
set(findall(gcf,'-property','FontSize'),'FontSize',12)

end


%% ARCHIVE ALT

% imshow(filtered_trials.odorID, "Colormap", [odor17_color; odor18_color; odor19_color],'DisplayRange', [17,19],'Border', 'tight')
% set(gcf,'OuterPosition',[0 100 300 350]);

% h = heatmap(filtered_trials.odorID, "Colormap", [odor17_color; odor18_color; odor19_color],'CellLabelColor','none');
% set color limits to avoid mis-coloring when not all odors are present
% clim([17 19]);