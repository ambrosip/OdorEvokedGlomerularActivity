% batch_analyzeBehavior SCRIPT

clear all

%% USER INPUTS 

mainDir = '/Users/priscilla/OHSU Dropbox/Priscilla Ambrosi/Dropbox - Moss Lab/Lab - Data/Olfactometer';
matchingFolders = findDateTimeFolders(mainDir);


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