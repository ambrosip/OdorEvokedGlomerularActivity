%{
photostim_SCRIPT


needs:
    scanimage > util > readTiffRoiData
%}


%% USER INPUT - experiment directory and others - EDIT ME

% experiment dir to be analyzed
tifFileDir = '/Users/ambrosi/Temp/photostim tests/raw/file_00004.tif';

% set img-specific inputs
photobleaching_window_s = 0; % duration of data in senconds that will be removed from baseline to account for photobleaching

% choose stim ROI to analyze
stimROItoAnalyze = 2;


%% GET STIM & SAVE DIRS

% get today's date
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% get expDir
[filepath,name,ext] = fileparts(tifFileDir);
expDir = filepath;
tifFileName = name;
% add sanity check for ext?

% get stimFileDir
stimFileDir = fullfile(expDir, strcat(tifFileName, '.stim'));
[filepath,name,ext] = fileparts(tifFileDir);
stimFileName = name;
% add sanity check for ext?
% add sanity check that file exists?

% create save dir if needed
saveDir = fullfile(expDir,'processed','matlab',analysisDate);
% check if saveDir exists
if not(isfolder(saveDir))
    % create saveDir
    mkdir(saveDir);
end


%% GET METADATA

% get ScanImage (SI) metadata WITHOUT ScanImageTiffReader
% SI data is stored in each frame
% ASSUMPTION: all frames in file have the same metadata
% ASSUMPTION CHECK further down
rawImgToAnalyze_metadata = imfinfo(tifFileDir);
rawImgToAnalyze_metadata_SI = rawImgToAnalyze_metadata(1).Software;
framesPerSlice = extractBetween(rawImgToAnalyze_metadata_SI,"SI.hStackManager.framesPerSlice = ",newline);
framesPerSlice = str2double(framesPerSlice{1});
scanFrameRate = extractBetween(rawImgToAnalyze_metadata_SI,"SI.hRoiManager.scanFrameRate = ",newline);
scanFrameRate = str2double(scanFrameRate{1});
laserPower = extractBetween(rawImgToAnalyze_metadata_SI,"SI.hBeams.powers = ",newline);
laserPower = str2num(laserPower{1});
laserPower920 = laserPower(1);
laserPower1040 = laserPower(2);
loopAcqInterval = extractBetween(rawImgToAnalyze_metadata_SI,"SI.loopAcqInterval = ",newline);
loopAcqInterval = str2num(loopAcqInterval{1});

% % get roi metadata WITH scanimage function "scanimage.util.readTiffRoiData"
% roiCenterInOpticalAngleXY = scanimage.util.readTiffRoiData(tifFileDir).imagingRois.scanfields.centerXY;
% roiSizeInOpticalAngleXY = scanimage.util.readTiffRoiData(tifFileDir).imagingRois.scanfields.sizeXY;
% roiSizeInPixelsXY = scanimage.util.readTiffRoiData(tifFileDir).imagingRois.scanfields.pixelResolutionXY;
% roiPixelsPerOpticalAngleXY = scanimage.util.readTiffRoiData(tifFileDir).imagingRois.scanfields.pixelRatio;

% get roi metadata WITHOUT scanimage functions
% 1. Store your raw character vector/string (assuming 'ans' holds the data)
jsonString = rawImgToAnalyze_metadata(1).Artist; 
% 2. Decode the JSON string into a MATLAB structure
jsonData = jsondecode(jsonString);
% 3. Get imaging roi metadata
roiCenterInOpticalAngleXY = jsonData.RoiGroups.imagingRoiGroup.rois.scanfields.centerXY';
roiSizeInOpticalAngleXY = jsonData.RoiGroups.imagingRoiGroup.rois.scanfields.sizeXY';
roiSizeInPixelsXY = jsonData.RoiGroups.imagingRoiGroup.rois.scanfields.pixelResolutionXY';
roiPixelsPerOpticalAngleXY = roiSizeInPixelsXY./roiSizeInOpticalAngleXY;
% 4. Get photostim metadata
% ASSUMPTIONs: photostim sequence contained, in order, 1 pause, 1 stim, 1
% park; photostim sequence was triggered by frame; park duration was set to
% a number larger than the grab duration, so that only one photostim would
% be delivered to the sample 
stimCenterInOpticalAngleXY =  jsonData.RoiGroups.photostimRoiGroups.rois(stimROItoAnalyze).scanfields.centerXY';
stimSizeInOpticalAngleXY = jsonData.RoiGroups.photostimRoiGroups.rois(stimROItoAnalyze).scanfields.sizeXY';
stimDurInSec = jsonData.RoiGroups.photostimRoiGroups.rois(stimROItoAnalyze).scanfields.duration;
stimType = split(jsonData.RoiGroups.photostimRoiGroups.rois(stimROItoAnalyze).scanfields.stimulusFunction,'.');
stimType = stimType{end};
if stimType == 'logspiral'
    logspiralRevolutions = jsonData.RoiGroups.photostimRoiGroups.rois(stimROItoAnalyze).scanfields.stimparams{2};
end
stimRepetitions = jsonData.RoiGroups.photostimRoiGroups.rois(stimROItoAnalyze).scanfields.repetitions;
stimPower = jsonData.RoiGroups.photostimRoiGroups.rois(stimROItoAnalyze).scanfields.powers;
stimCenterInPixels = stimCenterInOpticalAngleXY.*roiPixelsPerOpticalAngleXY;
stimSizeInPixels = stimSizeInOpticalAngleXY.*roiPixelsPerOpticalAngleXY;
% 5. Get pause duration before photostim
% ASSUMPTION: the same as above
pauseBeforeStimInSec = jsonData.RoiGroups.photostimRoiGroups.rois(1).scanfields.duration;

% need to make mask around stim area and plot F over time for that region
% also need to output image with roi overlayed


%{
"centerXY": [4.410888778,-2.396539569],
               "sizeXY": [12.71752604,16.42009691],

               pixelResolutionXY": [790,1020],
               "pixelToRefTransform": [
                 [0.01609813423,0,-1.955923308],
                 [0,0.01609813423,-10.61463709],
                 [0,0,1]
               ],
               "affine": [
                 [12.71752604,0,-1.947874241],
                 [0,16.42009691,-10.60658802],
                 [0,0,1]
%}
