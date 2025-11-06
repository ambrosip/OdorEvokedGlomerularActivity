%% Get dir of relevant files

% get today's date
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% get img dirs
mcorImgDir = fullfile(expDir,'processed','mcor');
rawImgDir = fullfile(expDir,'raw');

% create save dir if needed
saveDir = fullfile(expDir,'processed','matlab',analysisDate);
% check if saveDir exists
if not(isfolder(saveDir))
    % create saveDir
    mkdir(saveDir);
end

% improvemens with AI
% create one struct per "*Events.csv" file found in any expDir subfolder
% olfactometer_event_files = dir(fullfile(expDir, '**','*Events.csv'));

% Original line to get all files ending in 'Events.csv'
olfactometer_event_files = dir(fullfile(expDir, '**', '*Events.csv'));
olfactometer_event_files = remove_stupid_mac_hidden_files(olfactometer_event_files);
    
% order olfactometer files chronologically
olfactometer_event_files = olfactometer_event_files(~[olfactometer_event_files.isdir]);
[~,idx] = sort([olfactometer_event_files.datenum]);
olfactometer_event_files = olfactometer_event_files(idx);

% get scope h5 file dir
h5_file_dirs = dir(fullfile(expDir, '*.h5'));
h5_file_dirs = remove_stupid_mac_hidden_files(h5_file_dirs);   
h5_file_name = h5_file_dirs.name;
h5_file_dir = fullfile(expDir,h5_file_name);

disp('got dirs')