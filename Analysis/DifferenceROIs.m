%{
Difference ROIs
written by VA on Apr/2026

GOAL:
    Find which ROIs touch the high points in the difference image

DEPENDS:
    - image_dF_* figures struct.
    - finalMask from saveGUIwork

%}

%% USER INPUT

firstProgramNumber = 1;
firstOdorId = 17;

secondProgramNumber = 1;
secondOdorId = 18;

% hit, miss, false, or na?
outcome = "na";

% We will check if any ROI touches the region where
% the difference of z-scores >= threshold
threshold = 1.0;


%% Load relevant figures

function isIt = isRightFig(s, odorID, programNumber)
    isIt = s.odor == odorID && s.programNumber == programNumber;
end

findFun = arrayfun( ...
    @(s) isRightFig(s, firstOdorId, firstProgramNumber), figures);
firstOdorData = figures(find(findFun, 1));

findFun = arrayfun( ...
    @(s) isRightFig(s, secondOdorId, secondProgramNumber), figures);
secondOdorData = figures(find(findFun, 1));

diffZScore = abs(secondOdorData.(outcome).image ...
                 - firstOdorData.(outcome).image);

ROIs = utils.intersectingROIs(finalMask, diffZScore, threshold);

if isempty(ROIs)
    fprintf("No ROIs found that match the parameters!")
else
    msg = "ROIs where absolute difference z-score is > %f:\n";
    fprintf(msg, threshold);
    
    disp(ROIs);
end

