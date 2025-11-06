function plot_odorID(filtered_trials, odor17_color, odor18_color, odor19_color, hide_title, hide_x_labels, hide_y_labels)

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

% remove colorbar
h.ColorbarVisible = 'off';

% add custom x label
customXLabels = {'A', 'B', 'C'};
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