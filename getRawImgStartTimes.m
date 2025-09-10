% get directory of raw imgs
rawImgFileDirs = dir(fullfile(rawImgDir, '*.tif'));
rawImgFileNames = {rawImgFileDirs.name}';

% get acq start time and loop start time from raw img metadata
% "rawImgLoopStart" is orgaized as [year month day hour min sec]
rawImgStarts_sec = zeros(1,size(rawImgFileNames,1))';
rawImgLoopStart = zeros(6,size(rawImgFileNames,1))';
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
end

rawImgStarts_min = rawImgStarts_sec/60;

for i=1:size(rawImgLoopStart,2)
    if numel(unique(rawImgLoopStart(:,i))==1)
        disp('sanity check passed')
    else
        disp('we have raw img files from different loops')
    end
end
