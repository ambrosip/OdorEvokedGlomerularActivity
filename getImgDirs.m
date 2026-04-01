% getImgDirs

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

% analyze motion corrected (mcor) data or raw data based on user input
if analyzeMcorImgs == 1
    % get all tif file names in mcorImgDir
    imgsToAnalyzeDirs = dir(fullfile(mcorImgDir, '*.tif'));
    imgsToAnalyzeDirs = remove_stupid_mac_hidden_files(imgsToAnalyzeDirs);
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