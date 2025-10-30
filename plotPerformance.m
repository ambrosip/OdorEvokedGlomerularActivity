clear all
close all


%% USER INPUT

mainDir = '/Users/priscilla/OHSU Dropbox/Priscilla Ambrosi/Dropbox - Moss Lab/Lab - Behavior/m291/2025_10_29-11_50_02';
first_session = 1;
last_session = 8;


%% MAIN

% colorblind-safe colors (source:
% https://www.nature.com/articles/nmeth.1618)
black_color = [0 0 0]/255;
orange_color = [230 159 0]/255;
sky_blue_color = [86 180 233]/255;
bluish_green_color = [0 158 115]/255;
yellow_color = [240 228 66]/255;
blue_color = [0 114 178]/255;
vermillion_color = [213 94 0]/255;
reddish_purple_color = [204 121 167]/255;

% get today's date
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% Get dir and name of excel file inside mainDir
dataDir = dir(fullfile(mainDir, '*.xlsx'));
dataFolder = {dataDir.folder}';
dataName = {dataDir.name}';
dataToAnalyze = fullfile(dataFolder{1}, dataName{1});

% Get data from sheet/tab called "Results"
sheets = sheetnames(dataToAnalyze);
for sheetNum = 1:length(sheets)
    if sheets(sheetNum)=='Results'
        data = readtable(dataToAnalyze,'Sheet','Results');
    end
end

% Get performance data
hits_per_session = data.hits_per_trial(first_session:last_session)/100;
false_choices_per_session = data.false_per_trial(first_session:last_session)/100;
misses_per_session = data.miss_per_trial(first_session:last_session)/100;


%% PLOT

fig = figure('name', strcat(dataName{1}, '_', analysisDate, ' - performance'));
    hold on;
    plot(hits_per_session,'Color',bluish_green_color,'LineWidth',1,'DisplayName','Hits');
    plot(false_choices_per_session,'Color',vermillion_color,'LineWidth',1,'DisplayName','False choices');
    plot(misses_per_session,'Color',black_color,'LineWidth',1,'DisplayName','Misses');
    yline(0.5,'--')
    axis([1 length(first_session:last_session) 0 1])
    yticks([0,1]);
    xticks([0,length(first_session:last_session)]);
    xlabel('Session');
    ylabel('Events/Trials');
    title([dataName{1}, strcat("analyzed_on_", analysisDate)], 'Interpreter','none');
    hold off;
    legend('Hits', 'False choices', 'Misses','Location','northwest')
    legend('boxoff')
    set(fig, 'Position', [100 100 200 300])    % x y width height


%% SAVE

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

% save all open figs
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName = FigList(iFig).Name;
  set(0, 'CurrentFigure', FigHandle);
  % forces matlab to save fig as a vector
  FigHandle.Renderer = 'painters';  
  % actually saves a vector file
  saveas(FigHandle,fullfile(mainDir, [FigName '.svg']));
end 
disp('saved all figs')
close all