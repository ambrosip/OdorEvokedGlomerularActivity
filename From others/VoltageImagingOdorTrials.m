%% Set common parameters

% Image dimensions

exp = 'Lineocells'
imgHeight = 240;
imgWidth = 120;
dirname = uigetdir('*', 'Select Experiment Directory');
filelist = dir(dirname);
filelist=filelist(~ismember({filelist.name}, {'.', '..','.DS_Store'}));

%% Extract ROI information
 img = zeros(1,imgWidth*imgHeight);
 for j = 1:size(spatial_comp_refined,2)
     [idxx{j}, idxy{j}, val{j}] = find(spatial_comp_refined(:,j)>0);
     idx = idxx{j};
     img(idx) = img(idx)+1;
     clear idx
 end

% Display all ROIs together
imgg = reshape(img,240,120);
figure
imagesc(imgg) 

%% Export to imagej format

%%% Convert to imagej ROIs

% Your cell array of ROIs, where each cell contains pixel indices for one ROI
roiPixelIndices = idxx;

% Create a directory to save individual ROI mask files
outputDir = sprintf('%s_ROIs',exp);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for i = 1:numel(roiPixelIndices)
    % Create a blank binary mask for each ROI
    binaryMaskVector = zeros(1,(imgHeight*imgWidth));
    
    % Get the pixel indices for the current ROI
    indices = roiPixelIndices{i};
    
    % Set the ROI pixels to 1 in the binary mask
    binaryMaskVector(indices) = 1;

    % Reshape to 2D image
    binaryMask = reshape(binaryMaskVector,imgHeight,imgWidth);
    
    % Save each binary mask as an individual TIFF file
    imwrite(binaryMask, fullfile(outputDir, sprintf('ROI_%d.tif', i)));
end

%% Apply ROIs to new images and extract traces 

% Initialize structure
TrialData = [];

for j = 1:numel(filelist)

    pathtofile = sprintf('%s/%s',dirname,filelist(j).name);

   info = imfinfo(pathtofile);   % Get metadata for each frame
numFrames = numel(info);    % Number of frames in the TIFF stack

% Preallocate a 3D array for image data
imageHeight = info(1).Height;
imageWidth = info(1).Width;
imageData = zeros(imageHeight, imageWidth, numFrames, 'uint16');  % Adjust type if needed

% Load each frame
for k = 1:numFrames
    imageData(:,:,k) = imread(pathtofile, k, 'Info', info);
end

for k = 1:numFrames
    imageData(:,:,k) ;
end

%Loop through ROIs and create F_raw_refined variable


%%% Invert and process voltage traces

% Invert
 F = F_raw_refined*-1;

% Calculate df/f and z score using overall mean as F0
 for i = 1:size(F,1)
     tmp = F(i,:);
     df = tmp-(nanmean(tmp));
     df_f(i,:) = df/nanmean(tmp);
     zscore(i,:) = df/std(tmp);
     zscore_stag(i,:) = zscore(i,:)+(2*i);
 end

TrialData(j).File = filelist(j).name;
TrialData(j).RawData = F_raw_refined;
TrialData(j).DF_F = df_f;
TrialData(j).Zscore = zscore;

 % Plot traces staggered
 figure
 plot(zscore_stag')

