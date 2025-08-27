%{ 
DOCUMENTATION
written by PA on Aug/2025

GOAL:
    Show dF/F in each ROI organized by odor (rows), program type (columns),
    and trial outcome (color). Hits are green, False choices are orange,
    Misses are gray

ASSUMPTIONS:
    !! search for ALERT and ASSUMPTION to read important info

DEPENDS on:
    mat file created by timeSeriesFromFijiROI.m
    From others > ReadImageJROI (https://www.mathworks.com/matlabcentral/fileexchange/32479-readimagejroi)
    From others > ScanImageTiffReader

TO DO: 
    re-calculate mean_dF
%}


%% PLOT data in ROIs organized by odor (rows) and program type (columns)

% get max and min value of dF/F to set y axis limits
ymax = round(max(structfun(@(x) max(x,[],'all'),s_dF,'UniformOutput',true)), TieBreaker='plusinf');
ymin = round(min(structfun(@(x) min(x,[],'all'),s_dF,'UniformOutput',true)), TieBreaker='minusinf');

% get max and min value of xAxis to set x acis limits
xmin = round(min(xAxisInSec),TieBreaker='minusinf');
xmax = round(max(xAxisInSec),TieBreaker='plusinf');

% adjust ymin to some negative number in case it's zero
if ymin == 0
   ymin = -1;
end

% get max number of odors used in this experiment
max_odor_num = 0;
for programNum = 1:size(programFieldNames)
    programFieldName = programFieldNames(programNum);
    if s_olfactometer.(programFieldName).type ~= "ignore"
        odor_num = size(s_olfactometer.(programFieldName).odorList,1);
        if max_odor_num < odor_num
            max_odor_num = odor_num;
        end
    end
end

% iterate through rois
ignoredPrograms = size(programFieldNames,1) - programsToAnalyze;
for roi=1:rois_numberOf
    figName = strcat(firstAcqName(2:end), '_to_', lastAcqName(2:end), '_roi_', num2str(roi), '_dF');
    fig = figure('Name',figName);
    set(gca,'FontName','Arial');
    set(gcf,'OuterPosition',[100 100 1200 900]); % [left bottom width height]
    set(gca,'LineWidth', 0.75);
    % make a layout with "odor" rows and "program" + 1 columns
    rows = max_odor_num;
    columns = programsToAnalyze+1;
    t = tiledlayout(rows,columns);
    title(t,figName,'Interpreter','none');
    for programNum = 1:size(programFieldNames)
        programFieldName = programFieldNames(programNum);
        if s_olfactometer.(programFieldName).type ~= "ignore"
            for odorNum = 1:length(s_olfactometer.(programFieldName).odorList)
                nexttile(programNum - ignoredPrograms + (odorNum-1)*columns)
                odorID = extractBetween(s_olfactometer.(programFieldName).odorList(odorNum),"I "," -");
                odorFieldName = s_olfactometer.(programFieldName).odorFieldNames(odorNum);
                % color = odor_color.colorID(odor_color.odorID==str2double(odorID),:);
                hold on;
                rectangle('Position',[0 ymin odor_dur_s ymax-ymin],'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
                yline(0,'k--')
                axis([xmin xmax ymin ymax])
                title(strcat(odorFieldName, '_', s_olfactometer.(programFieldName).type), 'Interpreter','none');
                xlabel('Time from odor onset (s)')
                ylabel('dF/F')

                for acqIdx = s_olfactometer.(programFieldName).summary_by_trial.acqIdx(s_olfactometer.(programFieldName).summary_by_trial.odor==str2double(odorID))'
                    % plot dF/F   
                    if ~isnan(acqIdx)
                        if s_olfactometer.(programFieldName).summary_by_trial.outcome(s_olfactometer.(programFieldName).summary_by_trial.acqIdx==acqIdx) == "hit"
                            color = bluish_green_color;
                        elseif s_olfactometer.(programFieldName).summary_by_trial.outcome(s_olfactometer.(programFieldName).summary_by_trial.acqIdx==acqIdx) == "false choice"
                            color = vermillion_color;
                        else
                            color = [0.5 0.5 0.5];
                        end
                        plot(xAxisInSec',s_dF.(fns{acqIdx})(:,roi),'Color',[color 0.5],'LineWidth',0.5);  
                    end
                end

                % plot mean dF/F
                % plot(xAxisInSec',s_mean_dF.(odorFieldNames(odor))(:,roi),'Color',color,'LineWidth',1);
                hold off;

                disp(strcat("plot odor ", odorID, " done"))
            end
        end
    end
    % show ROI location
    nexttile(columns,[max_odor_num,1])
    imshow(imadjust(zProjFileToAnalyze,[0.5 0.65])) 
    hold on
    thetas = linspace(0,2*pi,200);
    ellipseR1 = (rois{roi}.vnRectBounds(4) - rois{roi}.vnRectBounds(2))/2;
    ellipseR2 = (rois{roi}.vnRectBounds(3) - rois{roi}.vnRectBounds(1))/2;
    ellipseA = (rois{roi}.vnRectBounds(4) + rois{roi}.vnRectBounds(2))/2;
    ellipseB = (rois{roi}.vnRectBounds(3) + rois{roi}.vnRectBounds(1))/2;
    ellipseX = ellipseR1*cos(thetas)+ellipseA;
    ellipseY = ellipseR2*sin(thetas)+ellipseB; 
    plot(ellipseX,ellipseY,'Color','y','LineWidth',1);
    hold off  
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    disp(strcat("plot roi ", num2str(roi), " done"))
end

disp("plot done")