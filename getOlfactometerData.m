% pre-allocate "programFieldNames": 1 row/program, 1 column
programFieldNames = strings([size(olfactometer_event_files,1),1]);

% set variables to 0
isFirstFine = 0;
isFirstCoarse = 0;
trialsToAnalyze = 0;
programsToAnalyze = 0;


%% Get olfactometer data based on olfactory task type

if olfactory_task == "2afc_fine_coarse_fine"    
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
        % find rows inside Events csv with trial starts
        trial_start_rows = matches(s_olfactometer.(programFieldName).file.Events,trial_start_label);
        % x axis in minutes
        % ALERT: comparing doubles can lead to errors 
        x_minutes = table2array(s_olfactometer.(programFieldName).file(:,1))/60/1000;
        % timestamps for trial-starts
        s_olfactometer.(programFieldName).startMin_by_trial = x_minutes(trial_start_rows);

        % label program types
        if contains(s_olfactometer.(programFieldName).name, 'fine', 'IgnoreCase', true) && ~contains(s_olfactometer.(programFieldName).name, 'pav', 'IgnoreCase', true)
            if isFirstFine == 0
                s_olfactometer.(programFieldName).type = "Fine 1";
                isFirstFine = 1;
                programsToAnalyze = programsToAnalyze + 1;
            else
                s_olfactometer.(programFieldName).type = "Fine 2";
                programsToAnalyze = programsToAnalyze + 1;
            end
        elseif contains(s_olfactometer.(programFieldName).name, 'coarse', 'IgnoreCase', true) && ~contains(s_olfactometer.(programFieldName).name, 'pav', 'IgnoreCase', true)
            if isFirstCoarse == 0
                s_olfactometer.(programFieldName).type = "Coarse 1";
                isFirstCoarse = 1;
                programsToAnalyze = programsToAnalyze + 1;
            else
                s_olfactometer.(programFieldName).type = "Coarse 2";
                programsToAnalyze = programsToAnalyze + 1;
            end
        else
            s_olfactometer.(programFieldName).type = "ignore";
        end

        % get further info from programs NOT labeled "ignore"
        if s_olfactometer.(programFieldName).type ~= "ignore"
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
                odorID = extractBetween(s_olfactometer.(programFieldName).odorList(odorNum),"I "," -");
                odorFieldName = strcat('odor_',odorID);
                s_olfactometer.(programFieldName).odorFieldNames(odorNum) = odorFieldName;
                % find rows inside Events csv with odors
                odor_start_rows = matches(s_olfactometer.(programFieldName).file.Events,s_olfactometer.(programFieldName).odorList(odorNum));
                % timestamps for odor presentation starts
                s_olfactometer.(programFieldName).(odorFieldName).startMin_by_odor = x_minutes(odor_start_rows);
                % save array with timestamp (in min) and odor presented
                odor_start_ts_labeled = x_minutes(odor_start_rows);
                odor_start_ts_labeled(:,2) = odorID;
                odor_start_ts_labeled_all = [odor_start_ts_labeled_all; odor_start_ts_labeled];
            end
            % sort array with timestamp (in min) and odor presented in chronological order (ie sort by timestamp)
            s_olfactometer.(programFieldName).odor_start_ts_labeled = sortrows(odor_start_ts_labeled_all);

            % label trial outcomes 
            % get total number of trials in this program
            trialNum_total = length(s_olfactometer.(programFieldName).startMin_by_trial);
            % get info to find response window within each trial
            response_window_start_idx = find(strcmp(s_olfactometer.(programFieldName).file.Events,'Response'));
            trial_interval_start_idx = find(strcmp(s_olfactometer.(programFieldName).file.Events,'Trial Interval'));
            % iterate through trials 
            for trialNum = 1:trialNum_total        
                firstIdx_this_trial = response_window_start_idx(trialNum);
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
                if ~isempty(find(contains(search_subset,'Reward'),1))
                    % if this trial was rewarded, mouse did the right
                    % action, and outcome is "hit"
                    s_olfactometer.(programFieldName).outcome_by_trial(trialNum,1) = "hit";                       
                elseif size(find(contains(search_subset,'Lick')),1) < minLicksToTriggerReward                                       
                    % if trial is NOT rewarded and mouse licked less than 
                    % minLicksToTriggerReward times, the outcome is "miss"
                    s_olfactometer.(programFieldName).outcome_by_trial(trialNum,1) = "miss";                    
                else
                    % if trial is NOT rewarded, but mouse licked more than minLicksToTriggerReward
                    % times, the mouse licked the wrong spout at some point, so the outcome
                    % is "false choice".
                    s_olfactometer.(programFieldName).outcome_by_trial(trialNum,1) = "false choice";                      
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
end


if olfactory_task == "passive_odor_presentations"
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
        % find rows inside Events csv with trial starts
        trial_start_rows = matches(s_olfactometer.(programFieldName).file.Events,trial_start_label);
        % x axis in minutes
        % ALERT: comparing doubles can lead to errors 
        x_minutes = table2array(s_olfactometer.(programFieldName).file(:,1))/60/1000;
        % timestamps for trial-starts
        s_olfactometer.(programFieldName).startMin_by_trial = x_minutes(trial_start_rows);

        % label program type as "na"
        s_olfactometer.(programFieldName).type = "na";
        programsToAnalyze = programsToAnalyze + 1;

        % label trial outcomes as "na"
        trialNum_total = length(s_olfactometer.(programFieldName).startMin_by_trial);
        s_olfactometer.(programFieldName).outcome_by_trial = strings(trialNum_total,1);
        s_olfactometer.(programFieldName).outcome_by_trial(:,:) = "na";

        % label odors used in each trial
        % get list of unique odors used in the program
        eventTypes = unique(s_olfactometer.(programFieldName).file.Events);
        pat = "Odor";
        odorIdx = contains(eventTypes,pat);
        s_olfactometer.(programFieldName).odorList = eventTypes(odorIdx);    
        % pre-allocate: 1 row/odor, 1 column
        s_olfactometer.(programFieldName).odorFieldNames = strings([length(s_olfactometer.(programFieldName).odorList),1]);
        % initiate empty array
        odor_start_ts_labeled_all = [];
        % iterate through odors
        for odorNum = 1:length(s_olfactometer.(programFieldName).odorList) 
            odorID = extractBetween(s_olfactometer.(programFieldName).odorList(odorNum),"I "," -");
            odorFieldName = strcat('odor_',odorID);
            s_olfactometer.(programFieldName).odorFieldNames(odorNum) = odorFieldName;
            % find rows inside Events csv with odors
            odor_start_rows = matches(s_olfactometer.(programFieldName).file.Events,s_olfactometer.(programFieldName).odorList(odorNum));
            % timestamps for odor presentation starts
            s_olfactometer.(programFieldName).(odorFieldName).startMin_by_odor = x_minutes(odor_start_rows);
            % save array with timestamp (in min) and odor presented
            odor_start_ts_labeled = x_minutes(odor_start_rows);
            odor_start_ts_labeled(:,2) = odorID;
            odor_start_ts_labeled_all = [odor_start_ts_labeled_all; odor_start_ts_labeled];
        end
        % sort array with timestamp (in min) and odor presented in chronological order (ie sort by timestamp)
        s_olfactometer.(programFieldName).odor_start_ts_labeled = sortrows(odor_start_ts_labeled_all);

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

            s_olfactometer.(programFieldName).startMin_by_trial = s_olfactometer.(programFieldName).startMin_by_trial(1:end-1);
        end
        
        % add to the number of trials I expect to analyze for this
        % experiment (i.e. the amount of acquisitions I expect to find
        % for this experiment)
        trialsToAnalyze = trialsToAnalyze + size(s_olfactometer.(programFieldName).summary_by_trial,1);
    end
end


%% Store olfactometer data into database

% create trials database
db_trials = table(...
    repmat(convertCharsToStrings(expDir),trialsToAnalyze,1), ...
    (1:trialsToAnalyze)', ...
    NaN(trialsToAnalyze,1), ...
    repmat("",trialsToAnalyze,1), ...
    NaN(trialsToAnalyze,1), ...
    NaN(trialsToAnalyze,1), ...
    NaN(trialsToAnalyze,1), ...
    repmat("",trialsToAnalyze,1), ...
    'VariableNames',{'expDir','trialNum','programNum','programType', ...
    'trial_start_min','odor_start_min','odorID','outcome'});    

% populate trials database
trialNum_total = 0;
trialNum_next = 1;
for programNum = 1:size(programFieldNames,1)
    programFieldName = programFieldNames(programNum);
    if s_olfactometer.(programFieldName).type ~= "ignore"
        trialNum_here = size(s_olfactometer.(programFieldName).summary_by_trial,1);
        trialNum_total = trialNum_total + trialNum_here;
        db_trials.programNum(trialNum_next:trialNum_total) = programNum;
        db_trials.programType(trialNum_next:trialNum_total) = s_olfactometer.(programFieldName).type;
        db_trials.trial_start_min(trialNum_next:trialNum_total) = s_olfactometer.(programFieldName).startMin_by_trial;
        db_trials.odor_start_min(trialNum_next:trialNum_total) = s_olfactometer.(programFieldName).summary_by_trial.min;
        db_trials.odorID(trialNum_next:trialNum_total) = s_olfactometer.(programFieldName).summary_by_trial.odor;
        db_trials.outcome(trialNum_next:trialNum_total) = s_olfactometer.(programFieldName).summary_by_trial.outcome;
        trialNum_next = trialNum_total + 1;
    end
end

disp('got olfactometer data')