%% Gather acq files

% analyze motion corrected (mcor) data or raw data based on user input
if analyzeMcorImgs == 1
    % get all tif file names in mcorImgDir
    imgsToAnalyzeDirs = dir(fullfile(mcorImgDir, '*.tif'));
    imgsToAnalyzeNames = {imgsToAnalyzeDirs.name}';
    imgsToAnalyzeFolder = mcorImgDir;
    imgsToAnalyze_numberOf = length(imgsToAnalyzeNames);
    fileNameIdxStart = 13;
    fileNameIdxEnd = 9;
else
    % get all tif file names in rawImgFileDirs
    imgsToAnalyzeDirs = rawImgFileDirs;
    imgsToAnalyzeNames = rawImgFileNames;
    imgsToAnalyzeFolder = rawImgDir;
    imgsToAnalyze_numberOf = length(imgsToAnalyzeNames); 
    fileNameIdxStart = 8;
    fileNameIdxEnd = 4;
end

% get acquisition # from img file name
% ASSUMPTION: file names end in "AAAAA_mcor.tif" where AAAAA is acquisition number
acq_list = [];
for file=1:imgsToAnalyze_numberOf
    img_file_name = cell2mat(imgsToAnalyzeNames(file));
    acq_list = [acq_list; img_file_name(end-fileNameIdxStart:end-fileNameIdxEnd)];
    acq_list = string(acq_list);
end

disp('got acquisition files')


%% Align Acqs to Trials

% adjust timing of acqs relative to start of scope loop (ie start of h5)
% the timing of the 1st acq is 0, so to align acqs to trials, I need to
% discount the time interval between the start of the loop and the start of
% the first trial
adjusted_rawImgStarts_min = rawImgStarts_min + trial_locs(1);

% plot for sanity checks
% FIG 1 - Scope h5 time series and peaks + Acqs starts
fig1 = figure('name', strcat(fileName_h5, '_', analysisDate, ' - scope events'));
hold on;
plot(x_minutes_h5,trial_start_TTL, 'Color','k')
plot(trial_locs,trial_pks,'o','Color','k')
plot(x_minutes_h5,odor_TTL, 'Color','m')
plot(odor_locs,odor_pks,'o','Color','m')
plot(odor_end_locs,odor_end_pks,'*','Color','m')
% show acq start times
xline(adjusted_rawImgStarts_min)
xlabel('Time (min)')
ylabel('TTLs (trial start, odor) and Events (acqs)','Interpreter','none')
hold off
disp('plot fig1 complete')


%% Map acqs to trials

% % get last idx of trials found with scope h5
% % if you have more trials than odor presentations, ignore the last trial
% if size(trial_locs,1) > size(odor_locs_idx,1)
%     trial_locs_idx = size(trial_locs,1) - 1;
% else
%     trial_locs_idx = size(trial_locs,1);
% end

% get last idx of acq list
acq_idx = size(acq_list,1);

% get last idx of trials found with scope h5
trial_locs_idx = size(trial_locs,1);

% if there is mismatch between the number of trials, annotate trials with missing acquisitions
% ALERT: potential for code breaking in case last trial was aborted mid-acquisition before odor was presented
if size(trial_locs,1) ~= size(rawImgStarts_min,1)  
    % iterate from last to 1st program
    for programNum = size(programFieldNames,1):-1:1         
        programFieldName = programFieldNames(programNum);   
        if s_olfactometer.(programFieldName).type ~= "ignore" 
            trialNum_total = size(s_olfactometer.(programFieldName).startMin_by_trial,1);
            % add columns pre-allocated with NaN where acq # per trial will go
            s_olfactometer.(programFieldName).summary_by_trial = addvars(s_olfactometer.(programFieldName).summary_by_trial,NaN(trialNum_total,1),'NewVariableName','acqNum');
            s_olfactometer.(programFieldName).summary_by_trial = addvars(s_olfactometer.(programFieldName).summary_by_trial,NaN(trialNum_total,1),'NewVariableName','acqIdx');
            s_olfactometer.(programFieldName).summary_by_trial = addvars(s_olfactometer.(programFieldName).summary_by_trial,NaN(trialNum_total,1),'NewVariableName','trial_locs_min');
            s_olfactometer.(programFieldName).summary_by_trial = addvars(s_olfactometer.(programFieldName).summary_by_trial,NaN(trialNum_total,1),'NewVariableName','odor_locs_min');
            % iterate from last to 1st trial
            for trialNum = trialNum_total:-1:1
                if abs(trial_locs(trial_locs_idx) - adjusted_rawImgStarts_min(acq_idx)) > tolerance
                    disp('oops we have trials without acq')
                    trial_locs_idx = trial_locs_idx - 1;
                else
                    s_olfactometer.(programFieldName).summary_by_trial.acqNum(trialNum) = str2double(acq_list(acq_idx));
                    s_olfactometer.(programFieldName).summary_by_trial.acqIdx(trialNum) = acq_idx;
                    s_olfactometer.(programFieldName).summary_by_trial.trial_locs_min(trialNum) = trial_locs(trial_locs_idx);
                    s_olfactometer.(programFieldName).summary_by_trial.odor_locs_min(trialNum) = odor_locs(trial_locs_idx);
                    trial_locs_idx = trial_locs_idx - 1;
                    acq_idx = acq_idx - 1;
                end
            end
        end
    end
else
    % iterate from last to 1st program
    for programNum = size(programFieldNames,1):-1:1         
        programFieldName = programFieldNames(programNum);   
        if s_olfactometer.(programFieldName).type ~= "ignore" 
            trialNum_total = size(s_olfactometer.(programFieldName).summary_by_trial,1);
            % add a column pre-allocated with NaN where acq # per trial will go
            s_olfactometer.(programFieldName).summary_by_trial = addvars(s_olfactometer.(programFieldName).summary_by_trial,NaN(trialNum_total,1),'NewVariableName','acqNum');
            s_olfactometer.(programFieldName).summary_by_trial = addvars(s_olfactometer.(programFieldName).summary_by_trial,NaN(trialNum_total,1),'NewVariableName','acqIdx');
            s_olfactometer.(programFieldName).summary_by_trial = addvars(s_olfactometer.(programFieldName).summary_by_trial,NaN(trialNum_total,1),'NewVariableName','trial_locs_min');
            s_olfactometer.(programFieldName).summary_by_trial = addvars(s_olfactometer.(programFieldName).summary_by_trial,NaN(trialNum_total,1),'NewVariableName','odor_locs_min');
            
            % if ignoreLastTrial == 1 % I'm not sure if this is still needed
            %     % had to add "if" statement here to handle case when the scope
            %     % loop was aborted mid-acquisition
            %     trialNum_total = trialNum_total - 1;
            % end
            
            % iterate from last to 1st trial
            for trialNum = trialNum_total:-1:1 

                % if acq_idx > 0 % I'm not sure if this is still needed
                %     % had to add "if" statement here to handle the case when
                %     % the scope loop started after the olfactometer, causing
                %     % the scope to miss the first trial
                %     s_olfactometer.(programFieldName).summary_by_trial.acqNum(trialNum) = str2double(acq_list(acq_idx));
                %     s_olfactometer.(programFieldName).summary_by_trial.acqIdx(trialNum) = acq_idx;
                %     acq_idx = acq_idx - 1;
                % end

                s_olfactometer.(programFieldName).summary_by_trial.acqNum(trialNum) = str2double(acq_list(acq_idx));
                s_olfactometer.(programFieldName).summary_by_trial.acqIdx(trialNum) = acq_idx;
                s_olfactometer.(programFieldName).summary_by_trial.trial_locs_min(trialNum) = trial_locs(trial_locs_idx);
                s_olfactometer.(programFieldName).summary_by_trial.odor_locs_min(trialNum) = odor_locs(trial_locs_idx);
                acq_idx = acq_idx - 1;
                trial_locs_idx = trial_locs_idx - 1;
            end
        end
    end
end


%% populate trials database

db_trials = addvars(db_trials,NaN(trialsToAnalyze,1),'NewVariableName','acqNum');
db_trials = addvars(db_trials,NaN(trialsToAnalyze,1),'NewVariableName','acqIdx');
db_trials = addvars(db_trials,NaN(trialsToAnalyze,1),'NewVariableName','trial_locs_min');
db_trials = addvars(db_trials,NaN(trialsToAnalyze,1),'NewVariableName','odor_locs_min');

trialNum_total = 0;
trialNum_next = 1;
for programNum = 1:size(programFieldNames,1)
    programFieldName = programFieldNames(programNum);
    if s_olfactometer.(programFieldName).type ~= "ignore"
        trialNum_here = size(s_olfactometer.(programFieldName).summary_by_trial,1);
        trialNum_total = trialNum_total + trialNum_here;
        db_trials.acqNum(trialNum_next:trialNum_total) = s_olfactometer.(programFieldName).summary_by_trial.acqNum;
        db_trials.acqIdx(trialNum_next:trialNum_total) = s_olfactometer.(programFieldName).summary_by_trial.acqIdx;
        db_trials.trial_locs_min(trialNum_next:trialNum_total) = s_olfactometer.(programFieldName).summary_by_trial.trial_locs_min;
        db_trials.odor_locs_min(trialNum_next:trialNum_total) = s_olfactometer.(programFieldName).summary_by_trial.odor_locs_min;
        trialNum_next = trialNum_total + 1;
    end
end


disp('finished matching olfactometer trials to acqs')