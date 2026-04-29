%% USER INPUT

saveFolder = 'C:\Temp\old\sid260\e1\processed\matlab\2026-04-28';

should_I_plot_dFF = 0;
visibility = 'off';

plotRoiSubset = 1;
roiSubset = [1:2];

% if using mcor files from odyn (caiman python)
convertFrom32bit = true;

saveFigs = 1;

manual_y_limits = 1;
ymax = 20;
ymin = -10;

manual_x_limits = 1;
xmin = -2;
xmax = 4;

ymin_cmPerSec = -50;
ymax_cmPerSec = 150;

filterData = 1;
smooth_span = 5;


% %% Organize wheel data into structure s_wheel
% % s_wheel.(programFieldName).(odorFieldName) = [xAxis_sec,cmPerSecondPerOdor];
% 
% % clear things
% clear s_wheel
% clear cmPerSecondPerOdor
% clear odorPulsePerOdor
% 
% % iterate through programs and odors
% totalOdorPresentations = 0;
% for programNum = 1:size(programFieldNames,1) % ALERT added ,1 
%     programFieldName = programFieldNames(programNum);
%     if s_olfactometer.(programFieldName).type ~= "ignore"
%         totalOdorPresentationsPerProgram = 0;
%         for odorNum = 1:length(s_olfactometer.(programFieldName).odorList)
%             odorID = extractBetween(s_olfactometer.(programFieldName).odorList(odorNum),"I "," -");
%             odorFieldName = s_olfactometer.(programFieldName).odorFieldNames(odorNum);
%             % find all acq idx that match this program and odor
%             % the idx of odor presentation will match the acq idx
%             odorPresentationNum = 0;
%             clear cmPerSecondPerOdor
%             clear odorPulsePerOdor
%             for acqIdx = s_olfactometer.(programFieldName).summary_by_trial.acqIdx(s_olfactometer.(programFieldName).summary_by_trial.odor==str2double(odorID))'
%                 % had to put this check to deal with problem file M:\ImagingData\20260316\m357\e1
%                 if acqIdx > 0
%                     odorPresentationNum = odorPresentationNum + 1; 
%                     totalOdorPresentationsPerProgram = totalOdorPresentationsPerProgram + 1;
%                     totalOdorPresentations = totalOdorPresentations + 1;
%                     odorOnset_dp = allOdorOnsets_dp(acqIdx);
%                     windowStart_dp = odorOnset_dp - round((baseline_dur_s - photobleaching_window_s) * si);
%                     windowEnd_dp = odorOnset_dp + round((odor_dur_s + baseline_dur_s - photobleaching_window_s) * si);
%                     cmPerSecondPerOdor(:,odorPresentationNum) = cmPerSecond_smooth(windowStart_dp:windowEnd_dp);
%                     odorPulsePerOdor(:,odorPresentationNum) = odorPulse(windowStart_dp:windowEnd_dp);
%                     xAxis_sec = linspace(0 - round(baseline_dur_s - photobleaching_window_s), round(odor_dur_s + baseline_dur_s - photobleaching_window_s), windowEnd_dp-windowStart_dp+1)';
%                 end
%             end
%             s_wheel.(programFieldName).(odorFieldName).allAcqs = cmPerSecondPerOdor;
%             s_wheel.(programFieldName).(odorFieldName).mean = mean(cmPerSecondPerOdor,2);
%             s_wheel.xAxis_sec = xAxis_sec;
%         end
%     end
% end


% %% Show ROIS
% 
% if plotRoiSubset == 1
%     roiRange = roiSubset;
% else
%     roiRange = 1:nROI;
% end
% 
% figName = strcat(firstAcqName(1:end), '_to_', lastAcqName(1:end), '_rois');
% fig = figure('Name',figName);
% 
% if convertFrom32bit
%     avgProjection = cast(avgProjection,'uint16');
%     imshow(imadjust(avgProjection))
% else
%     imshow(imadjust(avgProjection,[0.5 0.65])) 
% end
% hold on;
% finalMaskRGB = label2rgb(finalMask, 'jet', 'k', 'shuffle'); 
% hOverlay = imshow(finalMaskRGB);
% alphaMap = double(finalMask == roiRange) * 0.5; 
% set(hOverlay, 'AlphaData', alphaMap);
% hold off;
% 
% FigName = fig.Name;
% set(0, 'CurrentFigure', fig);
% % forces matlab to save fig as a vector
% fig.Renderer = 'painters';  
% % actually saves a vector file
% saveas(fig,fullfile(saveFolder, [FigName '.pdf']));
% close(fig);
% 
% hOverlay = imshow(labelRGB);
%     alphaMap = double(finalMask == iROI) * 0.5; 
%     set(hOverlay, 'AlphaData', alphaMap);
% 
%     hold off;  
% 
%     t.TileSpacing = 'compact';
%     t.Padding = 'compact';
%     disp(strcat("plot roi ", num2str(iROI), " done"))
% 



%% Prep for z plots

if manual_y_limits == 0
    % Gets the first integer below min and first above max
    ymax = ceil(max(structfun(@(x) max(x,[],'all'), ...
        fileZdFF,'UniformOutput',true)));
    ymin = floor(min(structfun(@(x) min(x,[],'all'), ...
        fileZdFF,'UniformOutput',true)));
end

if manual_x_limits == 0
    % Get the min and max value for x
    xmin = round(xs(1));
    xmax = round(xs(end));
end

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


%% PLOT AVGS - all outcomes (individual trials & mean), color-coded by outcome (or odor, if outcome is nan)

allProgramFieldName = fieldnames(s_zscore);

% (number of odors) x (number of programs + extra columns)
% The extra columns are to put the plot showing the ROI
extraColumns = 3;

% 2 rows per program - 1 for fluorescence, 1 for wheel data
rows = 2*length(allProgramFieldName);

% 1 column per odor + extra columns to show ROI
columns = maxOdorNumber + extraColumns;

% xsx is the palindromic version of xs (needed for the fill plot)
xsx = [xs flip(xs)];

if plotRoiSubset == 1
    roiRange = roiSubset;
else
    roiRange = 1:nROIs;
end

for iROI = roiRange
% for iROI = 1:nROIs
    % Create one figure per ROI
    figName = strcat( ...
        firstAcquisitionName, '_to_', lastAcquisitionName, ...
        '_roi_', num2str(iROI), '_z and wheel');

    fig = figure('Name', figName, 'Visible',visibility);
    % fig = figure('Name', figName, 'Visible',"on");

    % Create tiledlayout of shape odors by programs
    t = tiledlayout(rows, columns);
    title(t, figName, 'Interpreter','none');

    % Iterate through programs and odors, getting the fieldnames as we go
    for iProgram = 1:length(allProgramFieldName)
        % Get the struct to shorten the names in the next for loops
        programFieldName = allProgramFieldName{iProgram};        
        programStruct = s_zscore.(programFieldName);
        allOdorFieldName = fieldnames(programStruct);

        % Get program type to make tile title
        programType = s_olfactometer.(programFieldName).type;
        odorList = s_olfactometer.(programFieldName).odorList;

        for iOdor = 1:length(allOdorFieldName)
            odorFieldName = allOdorFieldName{iOdor};
            odorStruct = programStruct.(odorFieldName);
            allOutcomeFieldName = fieldnames(odorStruct);
            
            % plot fluorescence time series
            nexttile(iOdor + 2*(iProgram - 1) * columns)
            hold on;
            
            % Getting the odorID like Priscilla did
            odorID = extractBetween(odorList(iOdor), "I ", " -");

            % More plot settings + highlighting odor presentation window
            rectangle('Position',[0 ymin odor_dur_s ymax-ymin], ...
                'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
            % yline(0,'k--')
            
            % Iterate through outcomes, plotting both on the same graph
            for iOutcome = 1:length(allOutcomeFieldName)
                outcomeFieldName = allOutcomeFieldName{iOutcome};

                if olfactory_task == "passive_odor_presentations"
                    % Get the odor color from user input in preProcessing_v2
                    color = odor_color.colorID(odor_color.odorID==str2double(odorID),:);
                else
                    % Get the outcome color from struct in user input
                    color = outcomeColors.(outcomeFieldName);
                end

                % Gets data for this ROI and all acquisitions
                allFrames_allROIs_allAcq_z = odorStruct.(outcomeFieldName);     %%%%%%%%%%%%%%%
                allFrames_thisROI_allAcq = ...
                    allFrames_allROIs_allAcq_z(:,iROI,:);

                % Get the mean and along acquisitions
                allFrames_thisROI_meanAcq = ...
                    mean(allFrames_thisROI_allAcq, 3);

                % iterate through acquisitions
                for acqNum = 1:size(allFrames_thisROI_allAcq,3)
                    % Lines are a little less transparent and thicker than
                    % in other graphs (0.5 to 0.7 for both of them).
                    if filterData == 1
                        plot(xs', smooth(allFrames_thisROI_allAcq(:,1,acqNum),smooth_span), ...
                        'Color', [color 0.2], 'LineWidth', 0.5);
                    else
                        plot(xs', allFrames_thisROI_allAcq(:,1,acqNum), ...
                        'Color', [color 0.2], 'LineWidth', 0.5);
                    end
                    hold on;
                end
                plot(xs', allFrames_thisROI_meanAcq, ...
                        'Color', [color 1], 'LineWidth', 1.5);
            end

            hold off;
            disp(strcat("plot odor ", odorID, " done"))

            % Plot aesthetics     
            box off
            axis([xmin xmax ymin ymax])
            axis square
            if olfactory_task == "passive_odor_presentations"
                set(gcf, 'OuterPosition', [100 100 1600 900]);
            else 
                set(gcf, 'OuterPosition', [100 100 1600 900]);
            end
            set(gca, 'LineWidth', 0.75);
            set(gca, 'FontName', 'Arial');
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
            title(strcat(odorFieldName, '_', programType), ...
                'Interpreter', 'none');
            xlabel('Time from odor onset (s)')
            ylabel('z-score')
            xticks([xmin,0,1,xmax]);
            yticks([ymin,0,ymax]);   

            % plot wheel time series
            odorFieldName_alt = strcat("odor_", odorID);
            nexttile(iOdor + (2*iProgram - 1) * columns)
            plot(s_wheel.xAxis_sec,...
                s_wheel.(programFieldName).(odorFieldName_alt).allAcqs,...
                'Color', [color 0.2], 'LineWidth', 0.5);
            hold on;
            plot(s_wheel.xAxis_sec,...
                s_wheel.(programFieldName).(odorFieldName_alt).mean,...
                'Color', [color 1], 'LineWidth', 1.5);

            % beautification of wheel time series
            rectangle('Position',[0 ymin_cmPerSec odor_dur_s ymax_cmPerSec-ymin_cmPerSec], ...
                'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
            % yline(0,'k--')
            axis([xmin xmax ymin_cmPerSec ymax_cmPerSec])
            axis square
            box off

            if olfactory_task == "passive_odor_presentations"
                set(gcf, 'OuterPosition', [100 100 1600 700]);
            else 
                set(gcf, 'OuterPosition', [100 100 1600 900]);
            end
            set(gca, 'LineWidth', 0.75);
            set(gca, 'FontName', 'Arial');
            set(findall(gcf,'-property','FontSize'),'FontSize',12)
            xlabel('Time from odor onset (s)')
            ylabel('cm/s')
            xticks([xmin,0,1,xmax]);
            yticks([ymin_cmPerSec,0,ymax_cmPerSec]);   
        end
    end

    % Last tile showing the ROI
    nexttile(columns - extraColumns + 1, [rows, extraColumns])
    
    hold on;
    if convertFrom32bit
        avgProjection = cast(avgProjection,'uint16');
        imshow(imadjust(avgProjection));
    else
        imshow(imadjust(avgProjection,[0.5 0.65]));
    end
    
    hOverlay = imshow(labelRGB);
    alphaMap = double(finalMask == iROI) * 0.5; 
    set(hOverlay, 'AlphaData', alphaMap);

    hold off;  

    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    disp(strcat("plot roi ", num2str(iROI), " done"))

    if saveFigs
        FigName = fig.Name;
        set(0, 'CurrentFigure', fig);
        % forces matlab to save fig as a vector
        fig.Renderer = 'painters';  
        % actually saves a vector file
        saveas(fig,fullfile(saveFolder, [FigName '.svg']));
        close(fig);
    end

end


% %% Plot dFF (if user wants it)
% 
% if should_I_plot_dFF == 1
%     plot_dFF
% end


%% Save figs

% if saveFigs
%     FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% 
%     % save all open figs
%     for iFig = 1:length(FigList)
%       FigHandle = FigList(iFig);
%       FigName = FigList(iFig).Name;
%       set(0, 'CurrentFigure', FigHandle);
%       % forces matlab to save fig as a vector
%       FigHandle.Renderer = 'painters';  
%       % actually saves a vector file
%       saveas(FigHandle,fullfile(saveFolder, [FigName '.svg']));
%     end 
%     disp('saved all figs')
%     close all
% end
