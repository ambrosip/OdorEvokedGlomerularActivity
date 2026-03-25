% user input
timeInSec_label = "Var1";
cmPerSecond_label = "Var4";
odorPulse_label = "Var5";
cmDisplacement_label = "Var3";

% input from previous code
baseline_dur_s;
expDir;

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
xmin_min = firstOdorOnset_min - baseline_dur_s/60;
xmax_min = lastOdorOnset_min + interOdorInterval_min;

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

% plot cm/s
figure('name', strcat(wheelDataName, '_', analysisDate, '_treadmill'));
yyaxis left
plot(timeInMin,cmPerSecond)
hold on;
yline(20)
hold off;
ylabel('Velocity (cm/s)')
yyaxis right
plot(timeInMin,odorPulse)
ylabel('Odor Pulse (V)')
xlabel('Time (min)')
xlim([xmin_min xmax_min]);

% sanity check that treadmill data has as many odor presentations as the
% rest of the experiment
if size(allOdorOnsets_dp,1) == size(db_trials,1)
    disp("sanity check pass: matching number of odor presentations")
else
    disp("oops, we have mismatching number of odor presentations")
end

%% ALERT - approximating x axis
xAxis_sec = linspace(0 - round(baseline_dur_s - photobleaching_window_s), round(baseline_dur_s), round(baseline_dur_s/median(diff(timeInSec))))';

figure('name', strcat(wheelDataName, '_', analysisDate, '_odor psth'));
for odorPresentation = 1:size(allOdorOnsets_dp,1)
    odorOnset_dp = allOdorOnsets_dp(odorPresentation);
    windowStart_dp = odorOnset_dp - round((baseline_dur_s - photobleaching_window_s) * si);
    windowEnd_dp = odorOnset_dp + round(baseline_dur_s * si);
    cmPerSecondPerOdor(:,odorPresentation) = cmPerSecond(windowStart_dp:windowEnd_dp);
    odorPulsePerOdor(:,odorPresentation) = odorPulse(windowStart_dp:windowEnd_dp);
    xAxis_sec = linspace(0 - round(baseline_dur_s - photobleaching_window_s), round(baseline_dur_s), windowEnd_dp-windowStart_dp+1)';
    odorID = db_trials(odorPresentation,:).odorID;
    color = odor_color.colorID(odor_color.odorID==odorID,:);
    plot(xAxis_sec,cmPerSecondPerOdor(:,odorPresentation),'LineWidth',0.5,'Color',color)
    hold on
    % plot(xAxis_sec,odorPulsePerOdor(:,odorPresentation),'LineWidth',0.5,'Color',[1 0 0 0.5])
end
rectangle('Position',[0 -50 odor_dur_s 200],'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor', 'none');
ylabel('Velocity (cm/s)')
xlabel('Time from odor onset (s)')
hold off

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


%% ff

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    
% save all open figs
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName = FigList(iFig).Name;
  set(0, 'CurrentFigure', FigHandle);
  % forces matlab to save fig as a vector
  FigHandle.Renderer = 'painters';  
  % actually saves a vector file
  saveas(FigHandle,fullfile(saveDir, [FigName '.svg']));
end 
disp('I saved the figs')
close all
