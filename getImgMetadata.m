%% Get metadata from images

% get complete file dir and name for first image in rawImgDir
% ASSUMPTION: all figures in folder rawImgDir have the same number of frames
% and frame rate
rawImgFileDirs = dir(fullfile(rawImgDir, '*.tif'));
rawImgFileNames = {rawImgFileDirs.name}';
rawImgToAnalyzeFileDir = fullfile(rawImgFileDirs(1).folder, rawImgFileDirs(1).name);

% get ScanImage (SI) metadata WITHOUT ScanImageTiffReader
% SI data is stored in each frame
% ASSUMPTION: all frames in file have the same metadata
% ASSUMPTION CHECK further down
rawImgToAnalyze_metadata = imfinfo(rawImgToAnalyzeFileDir);
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

% get image dimensions from tiff metadata
% XResolution and YResolution are in pixels/cm. To convert resolution into
% pixels/um, I had to divide each by 10000.
width_px = rawImgToAnalyze_metadata(1).Width;
height_px = rawImgToAnalyze_metadata(1).Height;
width_um = rawImgToAnalyze_metadata(1).Width/(rawImgToAnalyze_metadata(1).XResolution/10000);
height_um = rawImgToAnalyze_metadata(1).Height/(rawImgToAnalyze_metadata(1).YResolution/10000);

% sanity check of image dimensions
if rawImgToAnalyze_metadata(1).ResolutionUnit == "Centimeter"
    disp('sanity check pass: width and height in um properly calculated')
else
    disp('ALERT: width and height in um is likely wrong')
end

% store variables in easier to understand names
frames_per_img = framesPerSlice;
frame_rate_hz = scanFrameRate;
time_between_acq_onset_in_s = loopAcqInterval;

% get acq start time and loop start time from raw img metadata
% "rawImgLoopStart" is orgaized as [year month day hour min sec]
% also get frames_per_img for all imgs for sanity checking
rawImgStarts_sec = zeros(1,size(rawImgFileNames,1))';
rawImgLoopStart = zeros(6,size(rawImgFileNames,1))';
rawImgFrames = zeros(1,size(rawImgFileNames,1))';
for rawImgIdx = 1:size(rawImgFileNames,1)
    currentImgDir = fullfile(rawImgFileDirs(rawImgIdx).folder, rawImgFileDirs(rawImgIdx).name);
    currentImgInfo = imfinfo(currentImgDir);
    % currentImgInfo(1) is the metadata from the first frame
    currentImgStart_sec = extractBetween(currentImgInfo(1).ImageDescription,"frameTimestamps_sec = ",newline);
    currentImgStart_sec = str2double(currentImgStart_sec{1});
    currentImgLoopStart = extractBetween(currentImgInfo(1).ImageDescription,"epoch = [","]");
    currentImgLoopStart = str2num(currentImgLoopStart{1});
    rawImgStarts_sec(rawImgIdx) = currentImgStart_sec;
    rawImgLoopStart(rawImgIdx,:) = currentImgLoopStart;
    % get number of frames for assumption check
    currentImgInfo_SI = currentImgInfo(1).Software;
    currentImg_framesPerSlice = extractBetween(rawImgToAnalyze_metadata_SI,"SI.hStackManager.framesPerSlice = ",newline);
    currentImg_framesPerSlice = str2double(currentImg_framesPerSlice{1});
    rawImgFrames(rawImgIdx) = currentImg_framesPerSlice;
end

rawImgStarts_min = rawImgStarts_sec/60;

% sanity check #1: raw imgs belong to same loop
for i=1:size(rawImgLoopStart,2)
    if numel(unique(rawImgLoopStart(:,i))==1)
        disp('sanity check passed: all raw imgs belong to same loop')
    else
        disp('ALERT: we have raw img files from different loops')
    end
end

% sanity check #2: raw imgs have the same number of frames
if numel(unique(rawImgFrames)==1)
    disp('sanity check passed: all raw imgs have the same number of frames')
else
    disp('ALERT: not all raw imgs have the same number of frames')
end

disp('got raw img metadata')