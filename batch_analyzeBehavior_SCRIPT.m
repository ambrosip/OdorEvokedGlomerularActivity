% batch_analyzeBehavior SCRIPT

%% USER INPUTS 

mainDir = 'M:\ImagingData';
% matchingFolders = findDateTimeFolders(mainDir);


%% MAIN

% get today's date
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% create date subfolder for saveDir
saveDir = fullfile(mainDir,'matlab-behavior',analysisDate);

% check if saveDir exists
if not(isfolder(saveDir))
    % create saveDir
    mkdir(saveDir);
end

% run the downstream code
processMatchingFolders(matchingFolders,saveDir)