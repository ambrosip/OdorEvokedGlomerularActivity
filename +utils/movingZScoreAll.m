function [FZScore, dFF, dFFZScoreTime, dFFZScoreAll] = ...
                movingZScoreAll(inputMovie, windowSize)
%COMPUTEZSCORE Returns a movie with the z-score of dF/F for each pixel.
%
%   This function creates a movie by computing dF/F of the average of a
%   sliding window, compared with the average of a baseline window of the
%   same size just before odor onset.

arguments
    inputMovie
    windowSize {mustBeInteger, mustBePositive}
end

% ASSUMPTIONS: 1) inputMovie has 3 dimensions.
%              2) Time dimension is the last one.

movingMean = movmean(inputMovie, windowSize, 3);

% movmean uses a smaller window when the frame is closer to the start/end
% The first frame that has a window with full size is the frame below
baselineMean = movingMean(:, :, ceil((windowSize + 1) / 2));
baselineStd = std(inputMovie(:, :, 1:windowSize), 0, 3);

% ALERT: 1) movingMean is now dF/F
%        2) ./ computes the entrywise division (don't use /)

% Compute the z-score of F
FZScore = (movingMean - baselineMean) ./ baselineStd;

% Compute dF/F
dFF = (movingMean - baselineMean) ./ baselineMean;

% Compute the z-score of dF/F (with respect to time)
dFFZScoreTime = (dFF - mean(dFF, 3)) / std(dFF, 0, 3);

% Compute the z-score of dF/F (with respect to time and space)
dFFZScoreAll = (dFF - mean(dFF(:))) / std(dFF(:));
end

