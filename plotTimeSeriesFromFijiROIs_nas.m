if olfactory_task ~= "passive_odor_presentations"
    error("check the olfactory_task")
end


%%

% Gets the first integer below min and first above max
ymax = ceil(max(structfun(@(x) max(x,[],'all'), ...
    fileZdFF,'UniformOutput',true)));
ymin = floor(min(structfun(@(x) min(x,[],'all'), ...
    fileZdFF,'UniformOutput',true)));

% Get the min and max value for x
xmin = round(xs(1));
xmax = round(xs(end));

% Get the max number of odors used in this experiment
maxOdorNumber = 0;
for iProgram = 1:length(programFieldNames)
    programFieldName = programFieldNames(iProgram);

    if s_olfactometer.(programFieldName).type ~= "ignore"
        odorNumber = length(s_olfactometer.(programFieldName).odorList);
        if maxOdorNumber < odorNumber
            maxOdorNumber = odorNumber;
        end
    end
end

% Number of rows and columns on the tiledlayout
rows = maxOdorNumber;
columns = programsToAnalyze + 1;


%% PLOT AVGS - nas only, color-coded by odor

allProgramFieldName = fieldnames(s_mean_zscore);

% (number of odors) x (number of programs + 1)
% The extra column is to put the plot showing the ROI
rows = maxOdorNumber;
columns = length(allProgramFieldName) + 1;

% % Struct with the odor colors, for easy reference
% odorColors = struct( ...
%     'hit', [56 83 163]/255, ...
%     'false_choice', [236 30 36]/255, ...
%     'miss',  [127 127 127]/255 ...
%     );

for iROI = 1:nROIs
    % Create one figure per ROI
    figName = strcat( ...
        firstAcquisitionName, '_to_', lastAcquisitionName, ...
        '_roi_', num2str(iROI), '_zdFF_AVG_nas');

    fig = figure('Name', figName);

    % Plot settings
    set(gca, 'FontName', 'Arial');
    set(gcf, 'OuterPosition', [100 100 1000 900]);
    set(gca, 'LineWidth', 0.75);

    % Create tiledlayout of shape odors by programs
    t = tiledlayout(rows, columns);
    title(t, figName, 'Interpreter','none');

    % Iterate through programs and odors, getting the fieldnames as we go
    for iProgram = 1:length(allProgramFieldName)
        % Get the struct to shorten the names in the next for loops
        programFieldName = allProgramFieldName{iProgram};        
        programStruct = s_mean_zscore.(programFieldName);
        allOdorFieldName = fieldnames(programStruct);

        % Get program type to make tile title
        programType = s_olfactometer.(programFieldName).type;

        for iOdor = 1:length(allOdorFieldName)
            odorFieldName = allOdorFieldName{iOdor};
            odorStruct = programStruct.(odorFieldName);
            allOutcomeFieldName = fieldnames(odorStruct);
            
            % Create the plot
            nexttile(iProgram + (iOdor - 1) * columns)
            hold on;
            
            % Getting the odorID like Priscilla did
            % ALERT - there is a better way to accomplish the next lines
            % that I should implement later... 
            odorList = s_olfactometer.(programFieldName).odorList;
            odorID = extractBetween(odorList(iOdor), "I ", " -");

            % More plot settings + highlighting odor presentation window
            rectangle('Position',[0 ymin odor_dur_s ymax-ymin], ...
                'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
            yline(0,'k--')

            axis([xmin xmax ymin ymax])
            title(strcat(odorFieldName, '_', programType), ...
                'Interpreter', 'none');
            xlabel('Time from odor onset (s)')
            ylabel('dF/F Z-score')
            
            % Iterate through outcomes, plotting both on the same graph
            for iOutcome = 1:length(allOutcomeFieldName)
                outcomeFieldName = allOutcomeFieldName{iOutcome};

                % This has all the ROIs so we get a slice below
                allFrames_allROIs = odorStruct.(outcomeFieldName);

                % Get the color from user input in preProcessing_v2
                color = odor_color.colorID(odor_color.odorID==str2double(odorID),:);

                % Lines are a little less transparent and thicker than
                % in other graphs (0.5 to 0.7 for both of them).
                plot(xs', allFrames_allROIs(:, iROI), ...
                    'Color', [color 0.7], 'LineWidth', 0.7);
            end

            hold off;
            disp(strcat("plot odor ", odorID, " done"))
        end
    end

    % Last tile showing the ROI
    nexttile(columns, [maxOdorNumber, 1])
    
    hold on
    imshow(imadjust(cast(zProjFileToAnalyze,'uint16')))

    thetas = linspace(0,2*pi,200);
    ellipseR1 = (ROIs{iROI}.vnRectBounds(4) - ROIs{iROI}.vnRectBounds(2))/2;
    ellipseR2 = (ROIs{iROI}.vnRectBounds(3) - ROIs{iROI}.vnRectBounds(1))/2;
    ellipseA = (ROIs{iROI}.vnRectBounds(4) + ROIs{iROI}.vnRectBounds(2))/2;
    ellipseB = (ROIs{iROI}.vnRectBounds(3) + ROIs{iROI}.vnRectBounds(1))/2;
    ellipseX = ellipseR1 * cos(thetas) + ellipseA;
    ellipseY = ellipseR2 * sin(thetas) + ellipseB; 
    plot(ellipseX, ellipseY, 'Color', 'y', 'LineWidth', 1);

    hold off  

    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    disp(strcat("plot roi ", num2str(iROI), " done"))
end