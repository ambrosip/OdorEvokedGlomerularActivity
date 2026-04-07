%% INPUTS

% user input
timeInSec_label = "Var1";
cmPerSecond_label = "Var4";
odorPulse_label = "Var5";
cmDisplacement_label = "Var3";
threshold_cmPerSec = 15;
ymin_cmPerSec = -50;
ymax_cmPerSec = 200;
smooth_span = 5;

% input from previous code
expDir;
baseline_dur_s;
analysisDate;


%% Pre-processing

% clear things
clear cmPerSecondPerOdor
clear odorPulsePerOdor

% get dir
wheelDataDir = dir(fullfile(expDir, '*.csv'));
wheelDataDir = remove_stupid_mac_hidden_files(wheelDataDir);
wheelDataDir = fullfile(wheelDataDir.folder,wheelDataDir.name);
[wheelDataFolder, wheelDataName, ~] = fileparts(wheelDataDir);
wheelData = readtable(wheelDataDir);

% get data
timeInSec = wheelData.(timeInSec_label);
cmPerSecond = wheelData.(cmPerSecond_label);
odorPulse = wheelData.(odorPulse_label);
cmDisplacement = wheelData.(cmDisplacement_label);

% exclude first data point (because of artifact)
timeInSec = timeInSec(2:end);
cmPerSecond = cmPerSecond(2:end);
odorPulse = odorPulse(2:end);
cmDisplacement = cmDisplacement(2:end);
timeInMin = timeInSec/60;
[cmPerSecond_envelope, ~] = envelope(cmPerSecond);

% filter cmPerSec
cmPerSecond_smooth = smooth(cmPerSecond, smooth_span);

% find odor onsets and offsets in data points (dp)
allOdorOnsets_dp = find(diff(odorPulse)>1);
allOdorOffsets_dp = find(diff(odorPulse)<-1);
firstOdorOnset_dp = allOdorOnsets_dp(1);
firstOdorOnset_sec = timeInSec(firstOdorOnset_dp);
firstOdorOnset_min = timeInMin(firstOdorOnset_dp);
lastOdorOnset_dp = allOdorOnsets_dp(end);
lastOdorOnset_sec = timeInSec(lastOdorOnset_dp);
lastOdorOnset_min = timeInMin(lastOdorOnset_dp);

% estimate inter-odor-interval
interOdorInterval_dp = round(mean(diff(allOdorOnsets_dp)));
interOdorInterval_sec = timeInSec(interOdorInterval_dp);
interOdorInterval_min = timeInMin(interOdorInterval_dp);

% estimate sampling interval
si = 1/median(diff(timeInSec));

% crop data based on odor presentation window
xmin_timeInMin = firstOdorOnset_min - baseline_dur_s/60;
xmax_timeInMin = lastOdorOnset_min + interOdorInterval_min;

% % plot envelope
% figure('name', strcat(wheelDataName, '_', analysisDate, '_treadmill env'));
% yyaxis left
% plot(timeInMin,cmPerSecond_envelope)
% hold on;
% yline(20)
% hold off;
% ylabel('Abs Velocity (cm/s)')
% yyaxis right
% plot(timeInMin,odorPulse)
% ylabel('Odor Pulse (V)')
% xlabel('Time (min)')
% xlim([xmin_min xmax_min]);

% % plot cm/s
% figure('name', strcat(wheelDataName, '_', analysisDate, '_treadmill'));
% yyaxis left
% plot(timeInMin,cmPerSecond)
% hold on;
% yline(threshold_cmPerSec, '--')
% hold off;
% ylim([ymin_cmPerSec ymax_cmPerSec]);
% ylabel('Velocity (cm/s)')
% yyaxis right
% plot(timeInMin,odorPulse)
% ylabel('Odor Pulse (V)')
% xlabel('Time (min)')
% xlim([xmin_timeInMin xmax_timeInMin]);
% title(wheelDataName,Interpreter="none");

% plot cm/s
figure('name', strcat(wheelDataName, '_', analysisDate, '_treadmill_smooth'));
yyaxis left
plot(timeInMin,cmPerSecond_smooth)
hold on;
yline(threshold_cmPerSec, '--')
hold off;
ylim([ymin_cmPerSec ymax_cmPerSec]);
ylabel('Velocity (cm/s)')
yyaxis right
plot(timeInMin,odorPulse)
ylabel('Odor Pulse (V)')
xlabel('Time (min)')
xlim([xmin_timeInMin xmax_timeInMin]);
title(wheelDataName,Interpreter="none");

% sanity check that treadmill data has as many odor presentations as the
% rest of the experiment
if size(allOdorOnsets_dp,1) == size(db_trials,1)
    disp("sanity check pass: matching number of odor presentations")
else
    disp("oops, we have mismatching number of odor presentations")
end


%% Plot locomotion during odor presentation 

% ALERT - approximating x axis
% xAxis_sec = linspace(0 - round(baseline_dur_s - photobleaching_window_s), round(baseline_dur_s), round(baseline_dur_s/median(diff(timeInSec))))';

% % plot raw data
% figure('name', strcat(wheelDataName, '_', analysisDate, '_odor psth'));
% for odorPresentation = 1:size(allOdorOnsets_dp,1)
%     odorOnset_dp = allOdorOnsets_dp(odorPresentation);
%     windowStart_dp = odorOnset_dp - round((baseline_dur_s - photobleaching_window_s) * si);
%     windowEnd_dp = odorOnset_dp + round(baseline_dur_s * si);
%     cmPerSecondPerOdor(:,odorPresentation) = cmPerSecond(windowStart_dp:windowEnd_dp);
%     odorPulsePerOdor(:,odorPresentation) = odorPulse(windowStart_dp:windowEnd_dp);
%     xAxis_sec = linspace(0 - round(baseline_dur_s - photobleaching_window_s), round(baseline_dur_s), windowEnd_dp-windowStart_dp+1)';
%     odorID = db_trials(odorPresentation,:).odorID;
%     color = odor_color.colorID(odor_color.odorID==odorID,:);
%     plot(xAxis_sec,cmPerSecondPerOdor(:,odorPresentation),'LineWidth',0.5,'Color',color)
%     hold on
%     % plot(xAxis_sec,odorPulsePerOdor(:,odorPresentation),'LineWidth',0.5,'Color',[1 0 0 0.5])
% end
% rectangle('Position',[0 ymin_cmPerSec odor_dur_s ymax_cmPerSec - ymin_cmPerSec],'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
% ylabel('Velocity (cm/s)')
% xlabel('Time from odor onset (s)')
% hold off

% plot smoothed data
figure('name', strcat(wheelDataName, '_', analysisDate, '_odor psth_smooth'));
for odorPresentation = 1:size(allOdorOnsets_dp,1)
    odorOnset_dp = allOdorOnsets_dp(odorPresentation);
    windowStart_dp = odorOnset_dp - round((baseline_dur_s - photobleaching_window_s) * si);
    windowEnd_dp = odorOnset_dp + round(baseline_dur_s * si);
    cmPerSecondPerOdor(:,odorPresentation) = cmPerSecond_smooth(windowStart_dp:windowEnd_dp);
    odorPulsePerOdor(:,odorPresentation) = odorPulse(windowStart_dp:windowEnd_dp);
    xAxis_sec = linspace(0 - round(baseline_dur_s - photobleaching_window_s), round(baseline_dur_s), windowEnd_dp-windowStart_dp+1)';
    odorID = db_trials(odorPresentation,:).odorID;
    color = odor_color.colorID(odor_color.odorID==odorID,:);
    plot(xAxis_sec,cmPerSecondPerOdor(:,odorPresentation),'LineWidth',0.5,'Color',color)
    hold on
    % plot(xAxis_sec,odorPulsePerOdor(:,odorPresentation),'LineWidth',0.5,'Color',[1 0 0 0.5])
end
rectangle('Position',[0 ymin_cmPerSec odor_dur_s ymax_cmPerSec - ymin_cmPerSec],'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
ylabel('Velocity (cm/s)')
xlabel('Time from odor onset (s)')
hold off


%% Plot smoothed data segregated by program and odor
% Also annotate trials with movement during and post odor presentation

% clear things
clear cmPerSecondPerOdor
clear odorPulsePerOdor

% add var to db_trials
% db_trials = addvars(db_trials,NaN(trialsToAnalyze,1),'NewVariableName','mov_dur_odor');
% db_trials = addvars(db_trials,NaN(trialsToAnalyze,1),'NewVariableName','mov_post_odor');

% get max number of programs used in this experiment
max_prog_num = size(programFieldNames,1);
fig_width = 200 * max_prog_num;

% get max number of odors used in this experiment
max_odor_num = 0;
for programNum = 1:size(programFieldNames,1) % ALERT added ",1" on Mar 19 2026
    programFieldName = programFieldNames(programNum);
    if s_olfactometer.(programFieldName).type ~= "ignore"
        odor_num = size(s_olfactometer.(programFieldName).odorList,1);
        if max_odor_num < odor_num
            max_odor_num = odor_num;
        end
    end
end
fig_height = 120 * max_odor_num;

% create fig and set general aesthetic parameters
fig = figure('Name',strcat(wheelDataName, '_', analysisDate, '_odor and program psth'));
set(gca,'FontName','Arial');
set(gcf,'OuterPosition',[100 100 fig_width fig_height]); % [left bottom width height]
set(gca,'LineWidth', 0.75);
xmin = 0 - round(baseline_dur_s - photobleaching_window_s);
xmax = round(odor_dur_s + baseline_dur_s - photobleaching_window_s);

% make a layout with "odor" rows and "program" columns
rows = max_odor_num;
columns = programsToAnalyze;
t = tiledlayout(rows,columns);
title(t,wheelDataName,'Interpreter','none');
relativeProgramNum = 0;

% iterate through programs and odors
for programNum = 1:size(programFieldNames,1) % ALERT added ,1 
    programFieldName = programFieldNames(programNum);
    if s_olfactometer.(programFieldName).type ~= "ignore"
        relativeProgramNum = relativeProgramNum + 1;
        for odorNum = 1:length(s_olfactometer.(programFieldName).odorList)
            nexttile(relativeProgramNum + (odorNum-1)*columns)
            odorID = extractBetween(s_olfactometer.(programFieldName).odorList(odorNum),"I "," -");
            odorFieldName = s_olfactometer.(programFieldName).odorFieldNames(odorNum);
            color = odor_color.colorID(odor_color.odorID==str2double(odorID),:);
            hold on;
            rectangle('Position',[0 ymin_cmPerSec odor_dur_s ymax_cmPerSec-ymin_cmPerSec],'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
            % yline(0,'k--')
            % yline(threshold_cmPerSec,'k--')
            axis([xmin xmax ymin_cmPerSec ymax_cmPerSec])
            title(strcat(odorFieldName, '_', s_olfactometer.(programFieldName).type), 'Interpreter','none');
            xlabel('Time from odor onset (s)')
            ylabel('Velocity (cm/s)')

            % find all acq idx that match this program and odor
            % the idx of odor presentation will match the acq idx
            for acqIdx = s_olfactometer.(programFieldName).summary_by_trial.acqIdx(s_olfactometer.(programFieldName).summary_by_trial.odor==str2double(odorID))'
                % had to put this check to deal with problem file M:\ImagingData\20260316\m357\e1
                if acqIdx > 0
                    odorPresentation = acqIdx; 
                    odorOnset_dp = allOdorOnsets_dp(odorPresentation);
                    windowStart_dp = odorOnset_dp - round((baseline_dur_s - photobleaching_window_s) * si);
                    windowEnd_dp = odorOnset_dp + round((odor_dur_s + baseline_dur_s - photobleaching_window_s) * si);
                    cmPerSecondPerOdor(:,odorPresentation) = cmPerSecond_smooth(windowStart_dp:windowEnd_dp);
                    odorPulsePerOdor(:,odorPresentation) = odorPulse(windowStart_dp:windowEnd_dp);
                    xAxis_sec = linspace(0 - round(baseline_dur_s - photobleaching_window_s), round(odor_dur_s + baseline_dur_s - photobleaching_window_s), windowEnd_dp-windowStart_dp+1)';
                    plot(xAxis_sec,cmPerSecondPerOdor(:,odorPresentation),'LineWidth',0.5,'Color',color)
                    hold on
                    % plot(xAxis_sec,50*odorPulsePerOdor(:,odorPresentation),'LineWidth',0.5,'Color',[1 0 0 0.5])
    
                    % detect movement during odor presentation and annotate db_trials
                    [pks,~] = findpeaks(cmPerSecond_smooth(allOdorOnsets_dp(odorPresentation):allOdorOffsets_dp(odorPresentation)),"MinPeakHeight",threshold_cmPerSec);
                    if ~isempty(pks)
                        db_trials.mov_dur_odor(odorPresentation) = 1;
                    else
                        db_trials.mov_dur_odor(odorPresentation) = 0;
                    end
    
                    % detect movement post odor presentation and annotate db_trials
                    [pks,~] = findpeaks(cmPerSecond_smooth(allOdorOffsets_dp(odorPresentation):round(allOdorOffsets_dp(odorPresentation)+odor_dur_s*si)),"MinPeakHeight",threshold_cmPerSec);
                    if ~isempty(pks)
                        db_trials.mov_post_odor(odorPresentation) = 1;
                    else
                        db_trials.mov_post_odor(odorPresentation) = 0;
                    end
                end
            end
            hold off;
            disp(strcat("plot odor ", odorID, " done"))
        end
    end
end


%% Detect movement onset irrespective of odor presentations


%% Save Figs

% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% 
% % save all open figs
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName = FigList(iFig).Name;
%   set(0, 'CurrentFigure', FigHandle);
%   % forces matlab to save fig as a vector
%   FigHandle.Renderer = 'painters';  
%   % actually saves a vector file
%   saveas(FigHandle,fullfile(saveDir, [FigName '.svg']));
% end 
% disp('I saved the figs')
% close all


%% Save workspace variables

% matFileName = strcat(analysisDate, '_treadmill');
% save(fullfile(saveDir,matFileName));     
% disp('I saved the mat file')


%% Troubleshooting code

% % check that db_trials labeled as "movement occurred" actually have movement
% figure
% for acqIdx = db_trials.acqIdx(db_trials.mov_post_odor==1)'
%     odorPresentation = acqIdx; 
%     odorOnset_dp = allOdorOnsets_dp(odorPresentation);
%     windowStart_dp = odorOnset_dp - round((baseline_dur_s - photobleaching_window_s) * si);
%     windowEnd_dp = odorOnset_dp + round((odor_dur_s + baseline_dur_s - photobleaching_window_s) * si);
%     xAxis_sec = linspace(0 - round(baseline_dur_s - photobleaching_window_s), round(odor_dur_s + baseline_dur_s - photobleaching_window_s), windowEnd_dp-windowStart_dp+1)';
%     plot(xAxis_sec, cmPerSecondPerOdor(:,odorPresentation),'LineWidth',0.5,'Color',color)
%     hold on
% end
% 
% figure
% for acqIdx = db_trials.acqIdx(db_trials.mov_dur_odor==1)'
%     odorPresentation = acqIdx; 
%     odorOnset_dp = allOdorOnsets_dp(odorPresentation);
%     windowStart_dp = odorOnset_dp - round((baseline_dur_s - photobleaching_window_s) * si);
%     windowEnd_dp = odorOnset_dp + round((odor_dur_s + baseline_dur_s - photobleaching_window_s) * si);
%     xAxis_sec = linspace(0 - round(baseline_dur_s - photobleaching_window_s), round(odor_dur_s + baseline_dur_s - photobleaching_window_s), windowEnd_dp-windowStart_dp+1)';
%     plot(xAxis_sec, cmPerSecondPerOdor(:,odorPresentation),'LineWidth',0.5,'Color',color)
%     hold on
% end


%% Abandoned code

% figure('name', strcat(wheelDataName, '_', analysisDate, '_odor psth segregated'));
% for odorID = unique(db_trials.odorID)
%     odorOnset_dp = allOdorOnsets_dp(odorPresentation);
%     windowStart_dp = odorOnset_dp - round((baseline_dur_s - photobleaching_window_s) * si);
%     windowEnd_dp = odorOnset_dp + round(baseline_dur_s * si);
%     cmPerSecondPerOdor(:,odorPresentation) = cmPerSecond(windowStart_dp:windowEnd_dp);
%     odorPulsePerOdor(:,odorPresentation) = odorPulse(windowStart_dp:windowEnd_dp);
%     xAxis_sec = linspace(0 - round(baseline_dur_s - photobleaching_window_s), round(baseline_dur_s), windowEnd_dp-windowStart_dp+1)';
%     odorID = db_trials(odorPresentation,:).odorID;
%     color = odor_color.colorID(odor_color.odorID==odorID,:);
%     plot(xAxis_sec,cmPerSecondPerOdor(:,odorPresentation),'LineWidth',0.5,'Color',color)
%     hold on
%     % plot(xAxis_sec,odorPulsePerOdor(:,odorPresentation),'LineWidth',0.5,'Color',[1 0 0 0.5])
% end
% ylabel('Velocity (cm/s)')
% xlabel('Time (s)')
% hold off
