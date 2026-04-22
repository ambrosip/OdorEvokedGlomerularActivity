%{
Odor Difference
written by VA on Apr/2026

GOAL:
    Plots the difference between the z-scores of dF/F for two odors

DEPENDS:
    - image_dF_* figures struct.

%}

close all
clear all

%% USER INPUT

% Save folder
saveFolder = "/Volumes/T7 Shield/PA/Poster";

% Odor numbers (order doesn't matter)
firstOdorId = 19;
secondOdorId = 20;

% hit, miss, false, or na?
outcome = "na";

% Type of image to compare
% Options: "image_dFF", "image_zScore_v1", "image_zScore_v2"
imageType = "image_zScore_v1";

% Range for the z-scores
zScoreMax = 5;

% Path to .mat files (could be more than one)
% USE DOUBLE QUOTES
pathMatFiles = ["/Volumes/T7 Shield/PA/Poster/20260316m357/2026-04-16/20260316_m357_e1_00001_to_20260316_m357_e1_00001_sequentialZScore.mat"];

% Program numbers
programNumbers = [1, 2];

% Define the colors
max_df_color = [103 0 31] / 255;
min_df_color = [5 48 97] / 255;

%% Creates Diverging Colormap

% Creates 512 color steps between the two colors (512 is arbitrary)
% The number just need to be high enough for the gradient to be smooth
colorpoints = linspace(0.0, 1.0, 512);

% This function is from the following website:
% https://www.kennethmoreland.com/color-maps/
divergingGradient = divergingMap(colorpoints, min_df_color, max_df_color);

%% Plot the figures

function isIt = isRightFig(s, odorID, programNumber)
    isIt = s.odor == odorID && s.programNumber == programNumber;
end

% Create figure
fig = figure;

nRows = length(pathMatFiles) * length(programNumbers);

tl = tiledlayout(fig, nRows, 3);
title(tl, "Z-score Change");

iRow = 1;

% Plot ranges
bidirectionalRange = [-zScoreMax zScoreMax];
onlyPositiveRange = [0 2*zScoreMax];

% One new row for each combination of experiment and program
for path = pathMatFiles
    load(path, 'figures', 'firstAcqName');

    if iRow == 1
        figName = split(firstAcqName, '_');
        figName = join(figName(1:3), '_');
        figName = figName{1};
        startFigName = true;

    else
        addToFigName = split(firstAcqName, '_');
        figName = sprintf("%s_%s", figName, addToFigName{3});

    end


    for programNumber = programNumbers
        % First odor plot
        ax = nexttile;
        findFun = arrayfun( ...
            @(s) isRightFig(s, firstOdorId, programNumber), figures);
        firstOdorData = figures(find(findFun, 1));

        imshow(firstOdorData.(outcome).(imageType), bidirectionalRange);

        clim(ax, bidirectionalRange);
        colormap(ax, divergingGradient);

        if iRow == 1
            title(sprintf("Odor %d", firstOdorId));

        elseif iRow == nRows
            colorbar(ax, 'southoutside');

        end

        % Second odor plot
        ax = nexttile;
        findFun = arrayfun( ...
            @(s) isRightFig(s, secondOdorId, programNumber), figures);
        secondOdorData = figures(find(findFun, 1));

        imshow(secondOdorData.(outcome).(imageType), bidirectionalRange);

        clim(ax, bidirectionalRange);
        colormap(ax, divergingGradient);

        if iRow == 1
            title(sprintf("Odor %d", secondOdorId));

        elseif iRow == nRows
            colorbar(ax, 'southoutside');

        end

        % Difference between images
        ax = nexttile;

        imshow( abs(secondOdorData.(outcome).(imageType) ...
                  - firstOdorData.(outcome).(imageType)), onlyPositiveRange);

        clim(ax, onlyPositiveRange);
        colormap(ax, flipud(gray));

        if iRow == 1
            title("Absolute Difference")

        elseif iRow == nRows
            colorbar(ax, 'southoutside');

        end

        iRow = iRow + 1;
    end
end

figName = strcat(figName, "_", num2str(firstOdorId), '_OdorDiff');
fig.Name = figName;

tl.TileSpacing = "tight";
tl.Padding = "compact";

set(gcf, 'Position', [2000 100 900 700])    % x y width height  


%% Save figure

figPath = fullfile(saveFolder, sprintf("%s.pdf", figName));
fig.Theme = 'light';
saveas(fig, figPath);
close all
