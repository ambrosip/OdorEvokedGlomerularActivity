%{
Segmentation GUI
written by VA on Apr/2026

GOAL:
    Give a quick and dirty solution to the segmentation problem

TODO:
    1) Improve the erosion, dilation, and watershed steps
    2) Add a tool to draw ROIs
    3) Include the whole movie instead of just a frame
    4) Improve the strip bounds calculation
    5) Add load state after closing
    6) Add border exclusion zone

DEPENDS:
    - image_dF_* until "Iterate over Data" section.

%}


% %% Compute Z-score envelope
%
% imageSize = size(figures(1).na.image);
%
% minEnvelope = zeros(imageSize);
% maxEnvelope = zeros(imageSize);
%
% for iFigure = 1:length(figures)
%     minEnvelope = min(minEnvelope, figures(iFigure).na.image);
%     maxEnvelope = max(maxEnvelope, figures(iFigure).na.image);
% end
%
% zScoreEnvelope = ...
%         (abs(minEnvelope) >= abs(maxEnvelope)) .* minEnvelope + ...
%         (abs(minEnvelope) < abs(maxEnvelope)) .* maxEnvelope;


%% Compute Z-score envelope (PA)
%
% zScoreEnvelope = figures(4).na.image;


%% Compute Z-score envelope

moviesMatFile = "M:\ImagingData\20260309\m357\processed\matlab\movies\20260309_m357_e1_correlation_array.mat";
zScoreEnvelope = load(moviesMatFile,'correlation_image');
zScoreEnvelope = zScoreEnvelope.correlation_image;

%% Create GUI

% ALERT: All relevant data is stored in fig.UserData
% We need to store data on the figure to be able to update the plot

% Create a UI Figure
segmentationGUI = uifigure();
segmentationGUI.Position = [800 150 1000 800];

% Create UserData struct with default values
segmentationGUI.UserData.original = zScoreEnvelope;
segmentationGUI.UserData.image = zScoreEnvelope;
segmentationGUI.UserData.blurStd = 1.0;
segmentationGUI.UserData.threshold = 1.0;
segmentationGUI.UserData.mask = zeros(size(zScoreEnvelope));
segmentationGUI.UserData.maskAfterExclusion = zeros(size(zScoreEnvelope));
segmentationGUI.UserData.plotRange = [-1 1];
segmentationGUI.UserData.gradient = divergingGradient;
segmentationGUI.UserData.showThreshold = false;
segmentationGUI.UserData.showROIs = true;
segmentationGUI.UserData.currentTool = "Pick Strip Center";
segmentationGUI.UserData.excludedCenter = ...
    round(size(zScoreEnvelope, 1) / 2);
segmentationGUI.UserData.excludedHeight = 50;
segmentationGUI.UserData.minArea = 0;
segmentationGUI.UserData.leftExcluded = 0;
segmentationGUI.UserData.rightExcluded = 0;
segmentationGUI.UserData.topExcluded = 0;
segmentationGUI.UserData.bottomExcluded = 0;
segmentationGUI.UserData.excludedROIs = [];

% Compute first mask
updatePlot(segmentationGUI);
updateMask(segmentationGUI);

% Create grid to organize widgets and plots
gl = uigridlayout(segmentationGUI, [3, 1]);
gl.RowHeight = {60, '1x', 120};

glButtons = uigridlayout(gl, [1, 3]);
glButtons.Layout.Row = 1;
glButtons.Layout.Column = 1;
glButtons.ColumnWidth = {'1x', '1x', '3x'};

showThreshold = uibutton(glButtons, "state");
showThreshold.Text = "Hide Threshold";
showThreshold.ValueChangedFcn = ...
    @(src, event) thresholdButton(src, segmentationGUI);

showROIs = uibutton(glButtons, "state");
showROIs.Text = "Hide ROIs";
showROIs.ValueChangedFcn = ...
    @(src, event) ROIButton(src, segmentationGUI);

tools = uibuttongroup(glButtons);
tools.SelectionChangedFcn = ...
    @(src, event) toolsCallback(src, segmentationGUI);

excludeRegionButton = uitogglebutton(tools);
excludeRegionButton.Position = [5 5 120 30];
excludeRegionButton.Text = "Pick Strip Center";

deleteROIButton = uitogglebutton(tools);
deleteROIButton.Position = [130 5 120 30];
deleteROIButton.Text = "Delete ROI";

addROIButton = uitogglebutton(tools);
addROIButton.Position = [260 5 120 30];
addROIButton.Text = "Add ROI";

segmentationGUI.UserData.axis = uiaxes(gl);
segmentationGUI.UserData.axis.Layout.Row = 2;
segmentationGUI.UserData.axis.Layout.Column = 1;

segmentationGUI.UserData.plotHandle = ...
    imshow( segmentationGUI.UserData.image ...
          , segmentationGUI.UserData.plotRange ...
          , 'Parent', segmentationGUI.UserData.axis ...
          );

colormap( segmentationGUI.UserData.axis ...
        , segmentationGUI.UserData.gradient);

% Add excluded region to the plot

stripCenter = segmentationGUI.UserData.excludedCenter;
stripWidth = size(segmentationGUI.UserData.image, 2);
stripHeight = segmentationGUI.UserData.excludedHeight;

stripXmin = 1;
stripYmin = stripCenter - round(stripHeight / 2);

excludedPosition = [ stripXmin stripYmin stripWidth stripHeight ];

segmentationGUI.UserData.excluded = ...
    rectangle( "Parent", segmentationGUI.UserData.axis ...
             , "FaceColor", [.8 .8 0.0] ...
             , "FaceAlpha", 0.5 ...
             , "Position", excludedPosition ...
             );

% Hold to plot overlay later
hold(segmentationGUI.UserData.axis, 'on');

% Add sliders
glSliders = uigridlayout(gl, [2, 8]);
glSliders.Layout.Row = 3;
glSliders.Layout.Column = 1;
glSliders.ColumnWidth = { 'fit', '1x' ...
                        , 'fit', '1x' ...
                        , 'fit', '1x' ...
                        , 'fit', '1x'
                        };

thresholdLabel = uilabel(glSliders);
thresholdLabel.Text = {'Min', 'Z-score'};
thresholdLabel.Layout.Row = 1;
thresholdLabel.Layout.Column = 1;

thresholdSlider = uislider(glSliders);
thresholdSlider.Layout.Row = 1;
thresholdSlider.Layout.Column = 2;
thresholdSlider.Limits = [0 1];
thresholdSlider.Value = segmentationGUI.UserData.threshold;
thresholdSlider.MajorTicks = 0:1;
thresholdSlider.MinorTicks = 0:.1:1;

thresholdSlider.ValueChangedFcn = ...
    @(src, event) updateThreshold(event, segmentationGUI);

blurLabel = uilabel(glSliders);
blurLabel.Text = 'Blur';
blurLabel.Layout.Row = 1;
blurLabel.Layout.Column = 3;

blurSlider = uislider(glSliders);
blurSlider.Layout.Row = 1;
blurSlider.Layout.Column = 4;
blurSlider.Limits = [0 2];
blurSlider.Value = segmentationGUI.UserData.blurStd;
blurSlider.MajorTicks = 0:2;
blurSlider.MinorTicks = 0:.1:2;

blurSlider.ValueChangedFcn = ...
    @(src, event) updateBlur(event, segmentationGUI);

excludedLabel = uilabel(glSliders);
excludedLabel.Text = {'Excluded', 'Strip Size'};
excludedLabel.Layout.Row = 1;
excludedLabel.Layout.Column = 5;

excludedSlider = uislider(glSliders);
excludedSlider.Layout.Row = 1;
excludedSlider.Layout.Column = 6;
excludedSlider.Limits = [0 400];
excludedSlider.Value = segmentationGUI.UserData.excludedHeight;
excludedSlider.MajorTicks = 0:100:400;
excludedSlider.MinorTicks = 0:50:400;

excludedSlider.ValueChangedFcn = ...
    @(src, event) updateExcludedRegion( segmentationGUI ...
                                      , height=event.Value ...
                                      );

minAreaLabel = uilabel(glSliders);
minAreaLabel.Text = {'Min', 'Area'};
minAreaLabel.Layout.Row = 1;
minAreaLabel.Layout.Column = 7;

minAreaSlider = uislider(glSliders);
minAreaSlider.Layout.Row = 1;
minAreaSlider.Layout.Column = 8;
minAreaSlider.Limits = [0 50];
minAreaSlider.Value = segmentationGUI.UserData.minArea;
minAreaSlider.MajorTicks = 0:10:50;
minAreaSlider.MinorTicks = 0:2.5:50;

minAreaSlider.ValueChangedFcn = ...
    @(src, event) updateExcludedRegion( segmentationGUI ...
                                      , minArea = event.Value ...
                                      );

leftExclusionLabel = uilabel(glSliders);
leftExclusionLabel.Text = {'Excluded', 'Left'};
leftExclusionLabel.Layout.Row = 2;
leftExclusionLabel.Layout.Column = 1;

leftExclusionSlider = uislider(glSliders);
leftExclusionSlider.Layout.Row = 2;
leftExclusionSlider.Layout.Column = 2;
leftExclusionSlider.Limits = [0 100];
leftExclusionSlider.Value = segmentationGUI.UserData.leftExcluded;
leftExclusionSlider.MajorTicks = 0:20:100;
leftExclusionSlider.MinorTicks = 0:5:100;

leftExclusionSlider.ValueChangedFcn = ...
    @(src, event) updateExcludedRegion( segmentationGUI ...
                                      , leftExcluded = event.Value ...
                                      );

rightExclusionLabel = uilabel(glSliders);
rightExclusionLabel.Text = {'Excluded', 'Right'};
rightExclusionLabel.Layout.Row = 2;
rightExclusionLabel.Layout.Column = 3;

rightExclusionSlider = uislider(glSliders);
rightExclusionSlider.Layout.Row = 2;
rightExclusionSlider.Layout.Column = 4;
rightExclusionSlider.Limits = [0 100];
rightExclusionSlider.Value = segmentationGUI.UserData.rightExcluded;
rightExclusionSlider.MajorTicks = 0:20:100;
rightExclusionSlider.MinorTicks = 0:5:100;

rightExclusionSlider.ValueChangedFcn = ...
    @(src, event) updateExcludedRegion( segmentationGUI ...
                                      , rightExcluded = event.Value ...
                                      );

topExclusionLabel = uilabel(glSliders);
topExclusionLabel.Text = {'Excluded', 'Top'};
topExclusionLabel.Layout.Row = 2;
topExclusionLabel.Layout.Column = 5;

topExclusionSlider = uislider(glSliders);
topExclusionSlider.Layout.Row = 2;
topExclusionSlider.Layout.Column = 6;
topExclusionSlider.Limits = [0 100];
topExclusionSlider.Value = segmentationGUI.UserData.topExcluded;
topExclusionSlider.MajorTicks = 0:20:100;
topExclusionSlider.MinorTicks = 0:5:100;

topExclusionSlider.ValueChangedFcn = ...
    @(src, event) updateExcludedRegion( segmentationGUI ...
                                      , topExcluded = event.Value ...
                                      );

bottomExclusionLabel = uilabel(glSliders);
bottomExclusionLabel.Text = {'Excluded', 'Bottom'};
bottomExclusionLabel.Layout.Row = 2;
bottomExclusionLabel.Layout.Column = 7;

bottomExclusionSlider = uislider(glSliders);
bottomExclusionSlider.Layout.Row = 2;
bottomExclusionSlider.Layout.Column = 8;
bottomExclusionSlider.Limits = [0 100];
bottomExclusionSlider.Value = segmentationGUI.UserData.bottomExcluded;
bottomExclusionSlider.MajorTicks = 0:20:100;
bottomExclusionSlider.MinorTicks = 0:5:100;

bottomExclusionSlider.ValueChangedFcn = ...
    @(src, event) updateExcludedRegion( segmentationGUI ...
                                      , bottomExcluded = event.Value ...
                                      );

% Add click callback
event.listener( segmentationGUI, 'WindowMousePress' ...
              , @(src, event) clickCallback(src, event, segmentationGUI));

% Add mask to plot
updateExcludedRegion(segmentationGUI);

segmentationGUI.UserData.overlayHandle = ...
    imshow( segmentationGUI.UserData.maskAfterExclusion ...
          , segmentationGUI.UserData.plotRange ...
          , 'Parent', segmentationGUI.UserData.axis);

colormap( segmentationGUI.UserData.axis ...
        , segmentationGUI.UserData.gradient);

redraw(segmentationGUI);

hold(segmentationGUI.UserData.axis, 'off');

%% Callbacks

function clickCallback(~, event, fig)
    if ~isequal(event.HitObject, fig.UserData.overlayHandle)
        return
    end

    pos = round(fig.UserData.axis.CurrentPoint(1,1:2));
    s = size(fig.UserData.mask');

    if all(pos > 0) && all(pos <= s)
        fprintf('Figure clicked at: X=%.2f, Y=%.2f\n', pos(:));

        if strcmp(fig.UserData.currentTool, "Delete ROI") ...
                && all(pos > 0) && all(pos <= s)
            iROI = fig.UserData.mask(pos(2), pos(1));
            fprintf('Clicked in ROI: %d\n', iROI);

            fig.UserData.mask(fig.UserData.mask == iROI) = 0;
            fig.UserData.maskAfterExclusion = ...
                min(fig.UserData.mask, fig.UserData.maskAfterExclusion);

            redraw(fig);

        elseif strcmp(fig.UserData.currentTool, "Pick Strip Center")
            fig.UserData.excludedCenter = pos(2);
            updateExcludedRegion(fig);

        elseif strcmp(fig.UserData.currentTool, "Add ROI")
            fprintf("Not implemented yet!\n");
        end
    end
end

function toolsCallback(src, fig)
    fig.UserData.currentTool = src.SelectedObject.Text;
end

function updateThreshold(event, fig)
    fig.UserData.threshold = event.Value;

    updatePlot(fig);
    updateMask(fig);
    updateExcludedRegion(fig);
    redraw(fig);
end

function updateBlur(event, fig)
    fig.UserData.blurStd = event.Value;

    updatePlot(fig);
    updateMask(fig);
    updateExcludedRegion(fig);
    redraw(fig);
end

function updateExcludedRegion(fig, options)
    arguments
        fig
        options.height = fig.UserData.excludedHeight
        options.minArea = fig.UserData.minArea
        options.leftExcluded = fig.UserData.leftExcluded
        options.rightExcluded = fig.UserData.rightExcluded
        options.topExcluded = fig.UserData.topExcluded
        options.bottomExcluded = fig.UserData.bottomExcluded
    end

    height = options.height;
    minArea = options.minArea;
    
    [imageHeight, imageWidth] = size(fig.UserData.mask);

    leftExcluded = max(1, round(options.leftExcluded));
    rightExcluded = imageWidth - round(options.rightExcluded);
    topExcluded = max(1, round(options.topExcluded));
    bottomExcluded = imageHeight - round(options.bottomExcluded);

    % Compute rectangle bounds
    xmin = 1;
    ymin = fig.UserData.excludedCenter - round(height / 2);

    width = size(fig.UserData.image, 2);
    height = round(height);

    % Update stored values
    fig.UserData.excludedHeight = height;
    fig.UserData.minArea = minArea;
    fig.UserData.leftExcluded = options.leftExcluded;
    fig.UserData.rightExcluded = options.rightExcluded;
    fig.UserData.topExcluded = options.topExcluded;
    fig.UserData.bottomExcluded = options.bottomExcluded;

    % ymin and height need to be such that the excluded region
    % is within the plot bounds.
    ymin = min(max(1, ymin), size(fig.UserData.mask, 1) - height);
    height = min(height, abs(size(fig.UserData.mask, 1) - ymin));

    % Update rectangle on the plot
    fig.UserData.excluded.Position = [ xmin ymin width height ];

    % Remove ROIs touching the strip
    U = unique(fig.UserData.mask(ymin : ymin + height, :));

    U = union(U, fig.UserData.mask(1:topExcluded, :));
    U = union(U, fig.UserData.mask(bottomExcluded:end, :));
    U = union(U, fig.UserData.mask(:, 1:leftExcluded));
    U = union(U, fig.UserData.mask(:, rightExcluded:end));

    % Exclude ROI that have area < minArea
    for iROI = 1:fig.UserData.nROIs
        if sum(fig.UserData.mask == iROI) < minArea
            U(end + 1) = iROI;
        end
    end

    fig.UserData.excludedROIs = U;

    fig.UserData.maskAfterExclusion = fig.UserData.mask;
    fig.UserData.maskAfterExclusion( ...
        ismember(fig.UserData.mask, fig.UserData.excludedROIs)) = 0;
    redraw(fig);
end

function thresholdButton(src, fig)
    if src.Value
        src.Text = "Show Threshold";
        fig.UserData.showThreshold = false;
    else
        src.Text = "Hide Threshold";
        fig.UserData.showThreshold = true;
    end

    updatePlot(fig);
    redraw(fig);
end

function ROIButton(src, fig)
    if src.Value
        src.Text = "Show ROIs";
        fig.UserData.showROIs = false;
    else
        src.Text = "Hide ROIs";
        fig.UserData.showROIs = true;
    end

    redraw(fig);
end

%% Plot update function

function updatePlot(fig)
    fig.UserData.image = fig.UserData.original;

    % Blur original image
    fig.UserData.image = ...
        imgaussfilt(fig.UserData.image, fig.UserData.blurStd);

    % Flatten low z-score regions
    if fig.UserData.showThreshold
        fig.UserData.image = ...
            (abs(fig.UserData.image) > fig.UserData.threshold) ...
            .* fig.UserData.image;
    end
end

function redraw(fig)
    % Update the Z-scores
    fig.UserData.plotHandle.CData = fig.UserData.image;
    fig.UserData.overlayHandle.CData = fig.UserData.maskColors;

    if fig.UserData.showROIs
        alphaMap = double(fig.UserData.maskAfterExclusion > 0) * 0.6;
    else
        alphaMap = zeros(size(fig.UserData.mask));
    end

    fig.UserData.overlayHandle.AlphaData = alphaMap;
end

function updateMask(fig)
    % Get latest plotted figure
    image = abs(fig.UserData.image) > fig.UserData.threshold;

    % Distance from boundary
    fig.UserData.mask = - bwdist(~image);

    % Make regions more smooth to split less
    fig.UserData.mask = medfilt2(fig.UserData.mask);

    % Split regions
    fig.UserData.mask = watershed(fig.UserData.mask, 4);
    fig.UserData.mask(~image) = 0.0;

    % Make regions rounder
    erosionDisk = strel("disk", 5);

    fig.UserData.mask = imerode(fig.UserData.mask, erosionDisk);
    fig.UserData.mask = imdilate(fig.UserData.mask, erosionDisk);

    [fig.UserData.mask, fig.UserData.nROIs] = ...
        bwlabel(fig.UserData.mask > 0, 4);

    fig.UserData.maskAfterExclusion = fig.UserData.mask;
    fig.UserData.excludedROIs = [];

    fig.UserData.maskColors = ...
        label2rgb(fig.UserData.mask, 'jet', 'k', 'shuffle');
end
