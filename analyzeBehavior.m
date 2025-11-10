%{ 
DOCUMENTATION
%}

% function analyzeBehavior(mouseDir,saveDir)

clear all

%% USER INPUT

mouseDir = '/Users/priscilla/OHSU Dropbox/Priscilla Ambrosi/Dropbox - Moss Lab/Lab - Behavior/m291';

% label relevant events from olfactometer program that are not
% automatically found
% trial_start_label = "Output 1"; % works for scope files, but not for lab files
trial_start_label = "Trial_";
minLicksToTriggerReward = 3;


%% Get dir of relevant files

% get today's date
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% create save dir if needed
if ~exist('saveDir', 'var')
    saveDir = fullfile(mouseDir,'matlab-behavior',analysisDate);
    disp('created saveDir')
end

% check if saveDir exists
if not(isfolder(saveDir))
    % create saveDir
    mkdir(saveDir);
end

% create one struct per "*Events.csv" file found in any mouseDir subfolder
olfactometer_event_files = dir(fullfile(mouseDir, '**','*Events.csv'));

% order olfactometer files chronologically
olfactometer_event_files = olfactometer_event_files(~[olfactometer_event_files.isdir]);
[~,idx] = sort([olfactometer_event_files.datenum]);
olfactometer_event_files = olfactometer_event_files(idx);

% get mouse number
[~,mouseNum,~] = fileparts(mouseDir);

disp('got dirs')


%% Get olfactometer data based on olfactory task type

% set variables to 0
trialsToAnalyze = 0;
programsToAnalyze = 0;

% pre-allocate "programFieldNames": 1 row/program, 1 column
programFieldNames = strings([size(olfactometer_event_files,1),1]);

% iterate through olfactometer programs
for programNum = 1:size(olfactometer_event_files,1)
    % common s_olfactometer structure building
    programFieldName = strcat("program_", num2str(programNum));
    programFieldNames(programNum) = programFieldName;
    s_olfactometer.(programFieldName).name = olfactometer_event_files(programNum).name;
    s_olfactometer.(programFieldName).folder = olfactometer_event_files(programNum).folder;
    s_olfactometer.(programFieldName).date = olfactometer_event_files(programNum).date;
    s_olfactometer.(programFieldName).dir = fullfile(s_olfactometer.(programFieldName).folder, s_olfactometer.(programFieldName).name);
    % load csv data and specify text data as string instead of char
    s_olfactometer.(programFieldName).file = readtable(s_olfactometer.(programFieldName).dir, TextType="string");
    s_olfactometer.(programFieldName).shortName = s_olfactometer.(programFieldName).name(end-29:end-4);    

    % let's focus on the 2AFC programs (3, 4 and 5 series)
    if contains(s_olfactometer.(programFieldName).name, '2afc', 'IgnoreCase', true)
        s_olfactometer.(programFieldName).type = "2afc";
    % let's add the mineral oil tests    
    elseif contains(s_olfactometer.(programFieldName).name, 'mineral', 'IgnoreCase', true)
        s_olfactometer.(programFieldName).type = "mineral";
    else
        s_olfactometer.(programFieldName).type = "ignore";
    end

    % get further info from programs NOT labeled "ignore"
    % also ignore stupid hidden files due to mac bs
    if s_olfactometer.(programFieldName).type ~= "ignore" &&...
            ~startsWith(s_olfactometer.(programFieldName).name, "._")

        % find rows inside Events csv with trial starts
        trial_start_rows = contains(s_olfactometer.(programFieldName).file.Events,trial_start_label);
        % find rows inside Events csv with licks
        lick_R_rows = contains(s_olfactometer.(programFieldName).file.Events,"Lick R");
        lick_L_rows = contains(s_olfactometer.(programFieldName).file.Events,"Lick L");
        % x axis in minutes
        % ALERT: comparing doubles can lead to errors 
        x_minutes = table2array(s_olfactometer.(programFieldName).file(:,1))/60/1000;
        x_seconds = table2array(s_olfactometer.(programFieldName).file(:,1))/1000;
        x_milliseconds = table2array(s_olfactometer.(programFieldName).file(:,1));    
        % timestamps for trial-starts
        s_olfactometer.(programFieldName).startMin_by_trial = x_minutes(trial_start_rows);
        % timestamps for licks
        s_olfactometer.(programFieldName).lickR_s = x_seconds(lick_R_rows);
        s_olfactometer.(programFieldName).lickL_s = x_seconds(lick_L_rows);

        % label odors used in each trial
        % get list of unique odors used in the program
        eventTypes = unique(s_olfactometer.(programFieldName).file.Events);
        pat = "Odor";
        odorIdx = contains(eventTypes,pat);
        s_olfactometer.(programFieldName).odorList = eventTypes(odorIdx);       
        % pre-allocate "odorFieldNames": 1 row/odor, 1 column
        s_olfactometer.(programFieldName).odorFieldNames = strings([length(s_olfactometer.(programFieldName).odorList),1]);
        % initiate empty array
        odor_start_ts_labeled_all = [];

        % iterate through odors
        for odorNum = 1:length(s_olfactometer.(programFieldName).odorList) 
            % ASSUMPTION: in the olfactometer odor setup, odors were named
            % "number - name".
            odorID = extractBetween(s_olfactometer.(programFieldName).odorList(odorNum),"I "," -");
            % deal with cases when odors were named in a different way
            if isempty(odorID)
                odorID = num2str(1000+odorNum);
                odorFieldName = strcat('odor_',odorID);
            else
                odorFieldName = strcat('odor_',odorID);
            end
            s_olfactometer.(programFieldName).odorFieldNames(odorNum) = odorFieldName;
            % find rows inside Events csv with odors
            odor_start_rows = matches(s_olfactometer.(programFieldName).file.Events,s_olfactometer.(programFieldName).odorList(odorNum));
            % timestamps for odor presentation starts
            s_olfactometer.(programFieldName).(odorFieldName).startMin_by_odor = x_minutes(odor_start_rows);
            % save array with timestamp (in min) and odor presented
            odor_start_ts_labeled = x_minutes(odor_start_rows);
            odor_start_ts_labeled(:,2) = str2num(odorID);
            odor_start_ts_labeled_all = [odor_start_ts_labeled_all; odor_start_ts_labeled];
        end

        % sort array with timestamp (in min) and odor presented in chronological order (ie sort by timestamp)
        s_olfactometer.(programFieldName).odor_start_ts_labeled = sortrows(odor_start_ts_labeled_all);

        % prep to label trial outcomes and collect lick data
        % get total number of trials in this program
        trialNum_total = length(s_olfactometer.(programFieldName).startMin_by_trial);
        % get info to find response window within each trial
        response_window_start_idx = find(strcmp(s_olfactometer.(programFieldName).file.Events,'Response'));
        trial_interval_start_idx = find(strcmp(s_olfactometer.(programFieldName).file.Events,'Trial Interval'));
        % get odor start within trial for calcularing lick latency from
        % odor onset
        odor_start_idx = find(contains(s_olfactometer.(programFieldName).file.Events,'Odor'));
        
        % iterate through trials 
        for trialNum = 1:trialNum_total 
            % had to add this in the case that the user ends the program
            % after the trial started (and odor was delivered) but before
            % response window started.
            if trialNum <= size(response_window_start_idx,1)         
                firstIdx_this_trial = response_window_start_idx(trialNum);
                odorIdx_this_trial = odor_start_idx(trialNum);
            elseif trialNum <= size(odor_start_idx,1)
                firstIdx_this_trial = size(x_minutes,1);
                odorIdx_this_trial = odor_start_idx(trialNum);
            else
                firstIdx_this_trial = size(x_minutes,1);
                odorIdx_this_trial = size(x_minutes,1);
            end
            % check if this is the last trial
            if trialNum == trialNum_total
                % it this is the last trial, set lastIdx_this_trial to last idx of x_minutes
                % I had to use this because the olfactometer does not
                % initiate a "trial interval" for the last trial in the
                % program
                lastIdx_this_trial = size(x_minutes,1);
            else
                % it this is NOT the last trial, use the  "trial 
                % interval" timestamp as lastIdx_this_trial
                lastIdx_this_trial = trial_interval_start_idx(trialNum);
            end
            % set a search subset from "Response" to "Trial Interval"
            % to find "Reward" and/or "Lick" within search subset to label trial outcomes
            search_subset = s_olfactometer.(programFieldName).file.Events(firstIdx_this_trial:lastIdx_this_trial);

            % label outcomes
            if ~isempty(find(contains(search_subset,'Reward'),1))
                % if this trial was rewarded, mouse did the right
                % action, and outcome is "hit"
                % ALERT this will errouneously label pavlovian trials as
                % "hits" no matter what the mouse did
                s_olfactometer.(programFieldName).outcome_by_trial(trialNum,1) = "hit"; 
                % label which side was rewarded
                if ~isempty(find(contains(search_subset,'Reward R'),1))
                    s_olfactometer.(programFieldName).rewarded_side(trialNum,1) = "R";
                else
                    s_olfactometer.(programFieldName).rewarded_side(trialNum,1) = "L";
                end
            elseif size(find(contains(search_subset,'Lick')),1) < minLicksToTriggerReward                                       
                % if trial is NOT rewarded and mouse licked less than 
                % minLicksToTriggerReward times, the outcome is "miss"
                s_olfactometer.(programFieldName).outcome_by_trial(trialNum,1) = "miss"; 
                s_olfactometer.(programFieldName).rewarded_side(trialNum,1) = "na";
            else
                % if trial is NOT rewarded, but mouse licked more than minLicksToTriggerReward
                % times, the mouse licked the wrong spout at some point, so the outcome
                % is "false choice".
                s_olfactometer.(programFieldName).outcome_by_trial(trialNum,1) = "false choice";   
                s_olfactometer.(programFieldName).rewarded_side(trialNum,1) = "na";
            end

            % set a search subset from "Odor onset" to "Trial Interval"
            % to find lick latency
            search_subset_from_odor = s_olfactometer.(programFieldName).file.Events(odorIdx_this_trial:lastIdx_this_trial);
            odor_ts_this_trial = x_seconds(odorIdx_this_trial);

            % measure lick latency and lick number (R licks) from odor
            % onset
            if ~isempty(find(contains(search_subset_from_odor,'Lick R'),1))
                all_lick_R_idx_in_search_subset = find(contains(search_subset_from_odor,'Lick R'));
                all_lick_R_idx_in_this_trial = odorIdx_this_trial + all_lick_R_idx_in_search_subset - 1;
                all_lick_R_ts_in_this_trial = x_seconds(all_lick_R_idx_in_this_trial);
                all_lick_R_latency_in_this_trial = all_lick_R_ts_in_this_trial - odor_ts_this_trial;
                first_lick_R_latency_in_this_trial = all_lick_R_latency_in_this_trial(1);
                lick_R_num = size(all_lick_R_idx_in_search_subset,1);
                s_olfactometer.(programFieldName).first_lick_R_latency_s(trialNum,1) = first_lick_R_latency_in_this_trial;
                s_olfactometer.(programFieldName).lick_R_total(trialNum,1) = lick_R_num;
                s_olfactometer.(programFieldName).all_lick_R_latency_per_trial(trialNum,1) = {all_lick_R_latency_in_this_trial};
            else
                s_olfactometer.(programFieldName).first_lick_R_latency_s(trialNum,1) = NaN;
                s_olfactometer.(programFieldName).lick_R_total(trialNum,1) = NaN;
                s_olfactometer.(programFieldName).all_lick_R_latency_per_trial(trialNum,1) = {NaN};
            end

            % measure lick latency and lick number (L licks) from odor
            % onset
            if ~isempty(find(contains(search_subset_from_odor,'Lick L'),1))
                all_lick_L_idx_in_search_subset = find(contains(search_subset_from_odor,'Lick L'));
                all_lick_L_idx_in_this_trial = odorIdx_this_trial + all_lick_L_idx_in_search_subset - 1;
                all_lick_L_ts_in_this_trial = x_seconds(all_lick_L_idx_in_this_trial);
                all_lick_L_latency_in_this_trial = all_lick_L_ts_in_this_trial - odor_ts_this_trial;
                first_lick_L_latency_in_this_trial = all_lick_L_latency_in_this_trial(1);
                lick_L_num = size(all_lick_L_idx_in_search_subset,1);
                s_olfactometer.(programFieldName).first_lick_L_latency_s(trialNum,1) = first_lick_L_latency_in_this_trial;
                s_olfactometer.(programFieldName).lick_L_total(trialNum,1) = lick_L_num;
                s_olfactometer.(programFieldName).all_lick_L_latency_per_trial(trialNum,1) = {all_lick_L_latency_in_this_trial};
            else
                s_olfactometer.(programFieldName).first_lick_L_latency_s(trialNum,1) = NaN;
                s_olfactometer.(programFieldName).lick_L_total(trialNum,1) = NaN;
                s_olfactometer.(programFieldName).all_lick_L_latency_per_trial(trialNum,1) = {NaN};
            end
        end

        % check if you have the same number of odor presentations and
        % trials. If not, user likely aborted the program after a trial
        % started, but before odor was presented. 
        if size(s_olfactometer.(programFieldName).odor_start_ts_labeled,1) == size(s_olfactometer.(programFieldName).outcome_by_trial,1)    
            % if the sizes match, keep all the info
            s_olfactometer.(programFieldName).summary_by_trial = table(...
                s_olfactometer.(programFieldName).odor_start_ts_labeled(:,1),...
                s_olfactometer.(programFieldName).odor_start_ts_labeled(:,2),...
                s_olfactometer.(programFieldName).outcome_by_trial,...
                'VariableNames',{'min','odor','outcome'});    
        else
            % if you have more trials than odor presentations, throw away
            % info about the last trial
            s_olfactometer.(programFieldName).summary_by_trial = table(...
                s_olfactometer.(programFieldName).odor_start_ts_labeled(:,1),...
                s_olfactometer.(programFieldName).odor_start_ts_labeled(:,2),...
                s_olfactometer.(programFieldName).outcome_by_trial(1:end-1),...
                'VariableNames',{'min','odor','outcome'});
        end
        
        % add to the number of trials I expect to analyze for this
        % experiment (i.e. the amount of acquisitions I expect to find
        % for this experiment)
        trialsToAnalyze = trialsToAnalyze + size(s_olfactometer.(programFieldName).summary_by_trial,1);
    end
end

disp('created s_olfactometer structure')


%% Store olfactometer data into database

% create trials database
db_trials = table(...
    repmat(convertCharsToStrings(mouseDir),trialsToAnalyze,1), ... % mouseDir
    repmat("",trialsToAnalyze,1), ... % mouse
    repmat("",trialsToAnalyze,1), ... % date
    repmat("",trialsToAnalyze,1), ... % start_time
    NaN(trialsToAnalyze,1), ... % programNum
    repmat("",trialsToAnalyze,1), ... % programType
    repmat("",trialsToAnalyze,1), ... % programName   
    (1:trialsToAnalyze)', ... % absolute trialNum     
    NaN(trialsToAnalyze,1), ... % relative trialNum
    NaN(trialsToAnalyze,1), ... % trial_start_min
    NaN(trialsToAnalyze,1), ... % odor_start_min
    NaN(trialsToAnalyze,1), ... % odorID
    repmat("",trialsToAnalyze,1), ... % outcome
    repmat("",trialsToAnalyze,1), ... % rewarded_side
    NaN(trialsToAnalyze,1), ... % R lick latency from odor onset
    NaN(trialsToAnalyze,1), ... % R lick total during odor and response window
    NaN(trialsToAnalyze,1), ... % L lick latency from odor onset
    NaN(trialsToAnalyze,1), ... % L lick total during odor and response window
    'VariableNames',{'mouseDir','mouse','date','start_time',...
    'programNum','programType','programName',...
    'trialNum_abs','trialNum_rel'...
    'trial_start_min','odor_start_min','odorID',...
    'outcome','rewarded_side',...
    'first_R_lick_latency_s_from_odor_onset',...
    'R_lick_total_during_odor_and_response',...
    'first_L_lick_latency_s_from_odor_onset',...
    'L_lick_total_during_odor_and_response'});    

% populate trials database
% abs: absolute
abs_trialNum_total = 0;
abs_trialNum_next = 1;
for programNum = 1:size(programFieldNames,1)
    programFieldName = programFieldNames(programNum);
    if s_olfactometer.(programFieldName).type ~= "ignore" &&...
            ~startsWith(s_olfactometer.(programFieldName).name, "._")
        trialNum_here = size(s_olfactometer.(programFieldName).summary_by_trial,1);
        abs_trialNum_total = abs_trialNum_total + trialNum_here;
        db_trials.trialNum_rel(abs_trialNum_next:abs_trialNum_total) = (1:size(s_olfactometer.(programFieldName).summary_by_trial,1))';
        db_trials.programNum(abs_trialNum_next:abs_trialNum_total) = programNum;
        db_trials.programType(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).type;
        db_trials.trial_start_min(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).startMin_by_trial;
        db_trials.odor_start_min(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).summary_by_trial.min;
        db_trials.odorID(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).summary_by_trial.odor;
        db_trials.outcome(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).summary_by_trial.outcome;
        db_trials.mouse(abs_trialNum_next:abs_trialNum_total) = convertCharsToStrings(mouseNum);
        db_trials.date(abs_trialNum_next:abs_trialNum_total) = convertCharsToStrings(s_olfactometer.(programFieldName).shortName(1:10));
        db_trials.start_time(abs_trialNum_next:abs_trialNum_total) = convertCharsToStrings(s_olfactometer.(programFieldName).shortName(12:19));
        db_trials.programName(abs_trialNum_next:abs_trialNum_total) = convertCharsToStrings(s_olfactometer.(programFieldName).name(1:end-31));        
        db_trials.first_R_lick_latency_s_from_odor_onset(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).first_lick_R_latency_s;
        db_trials.R_lick_total_during_odor_and_response(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).lick_R_total;
        db_trials.first_L_lick_latency_s_from_odor_onset(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).first_lick_L_latency_s;
        db_trials.L_lick_total_during_odor_and_response(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).lick_L_total;
        db_trials.rewarded_side(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).rewarded_side;
        db_trials.all_lick_R_latency_per_trial(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).all_lick_R_latency_per_trial;
        db_trials.all_lick_L_latency_per_trial(abs_trialNum_next:abs_trialNum_total) = s_olfactometer.(programFieldName).all_lick_L_latency_per_trial;
        db_trials.avg_trial_dur_s(abs_trialNum_next:abs_trialNum_total) = mean(diff(s_olfactometer.(programFieldName).startMin_by_trial))*60;
        abs_trialNum_next = abs_trialNum_total + 1;
    end
end

disp('created db_trials table (database)')


%% Save workspace

% save workspace variables
matFileName = strcat(analysisDate, '_', mouseNum, '_behavior');
save(fullfile(saveDir,matFileName));     
disp('I saved the mat file')


%% Get all dates to plot

datesToPlot = unique(db_trials.date);
for date = datesToPlot'
    try
        plotBehavior_taskSwitching(date, db_trials, saveDir, mouseNum)
    catch
        plotBehavior_v2(date, db_trials, saveDir, mouseNum)
        plotBehavior_appended(date, db_trials, saveDir, mouseNum)
    end   
    disp(strcat("plotting ", date))
end


%% ARCHIVE

            % % measure lick latency and lick number (R licks) during
            % % response window
            % if ~isempty(find(contains(search_subset,'Lick R'),1))
            %     first_lick_R_idx = find(contains(search_subset,'Lick R'),1);
            %     lick_R_num = size(find(contains(search_subset,'Lick R')),1);
            %     s_olfactometer.(programFieldName).lick_R_latency_sec(trialNum,1) = x_seconds(first_lick_R_idx);
            %     s_olfactometer.(programFieldName).lick_R_total(trialNum,1) = lick_R_num;
            % else
            %     s_olfactometer.(programFieldName).lick_R_latency_sec(trialNum,1) = NaN;
            %     s_olfactometer.(programFieldName).lick_R_total(trialNum,1) = NaN;
            % end
            % 
            % % measure lick latency and lick number (L licks) during
            % % response window
            % if ~isempty(find(contains(search_subset,'Lick L'),1))
            %     first_lick_L_idx = find(contains(search_subset,'Lick L'),1);
            %     lick_L_num = size(find(contains(search_subset,'Lick L')),1);
            %     s_olfactometer.(programFieldName).lick_L_latency_sec(trialNum,1) = x_seconds(first_lick_L_idx);
            %     s_olfactometer.(programFieldName).lick_L_total(trialNum,1) = lick_L_num;
            % else
            %     s_olfactometer.(programFieldName).lick_L_latency_sec(trialNum,1) = NaN;
            %     s_olfactometer.(programFieldName).lick_L_total(trialNum,1) = NaN;
            % end

            % archived code (before lick raster plot)
            % measure lick latency and lick number (R licks) from odor
            % onset
            % if ~isempty(find(contains(search_subset_from_odor,'Lick R'),1))
            %     first_lick_R_idx_search_subset = find(contains(search_subset_from_odor,'Lick R'),1);
            %     first_lick_R_idx_this_trial = odorIdx_this_trial + first_lick_R_idx_search_subset - 1;
            %     first_lick_R_ts_this_trial = x_seconds(first_lick_R_idx_this_trial);
            %     first_lick_R_latency = first_lick_R_ts_this_trial - odor_ts_this_trial;
            %     lick_R_num = size(find(contains(search_subset_from_odor,'Lick R')),1);
            %     s_olfactometer.(programFieldName).lick_R_latency_s(trialNum,1) = first_lick_R_latency;
            %     s_olfactometer.(programFieldName).lick_R_total(trialNum,1) = lick_R_num;
            % else
            %     s_olfactometer.(programFieldName).lick_R_latency_s(trialNum,1) = NaN;
            %     s_olfactometer.(programFieldName).lick_R_total(trialNum,1) = NaN;
            % end

% end