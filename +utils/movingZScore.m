function zScoreMovie = movingZScore(inputMovie, windowSize, frameOnset)
%COMPUTEZSCORE Returns a movie with the z-score of dF/F for each pixel.
%
%   This function creates a movie by computing dF/F of the average of a
%   sliding window, compared with the average of a baseline window of the
%   same size just before odor onset.

arguments
    inputMovie
    windowSize {mustBeInteger, mustBePositive}
    frameOnset {mustBeInteger, ...
                mustBeGreaterThan(frameOnset, windowSize)}
end

% ASSUMPTIONS: 1) inputMovie has 3 dimensions.
%              2) Time dimension is the last one.

zScoreMovie = movmean(inputMovie, windowSize, 3);

% movmean uses a smaller window when the frame is closer to the start/end
% The first frame that has a window with full size is the frame below
baseline = zScoreMovie(:, :, ceil((windowSize + 1) / 2));

% ASSUMPTION: Input was originally int16.
% 
% Thus, we will shift entries up to get positive values.
% This will make dF/F well-defined and finite (no division by zero).

% ALERT: 1) zScoreMovie is now dF/F
%        2) We didn't shift the arrays up in the numerator 
%           because it cancels out after computing the difference.
%        3) ./ computes the entrywise division (don't use /)
zScoreMovie = (zScoreMovie - baseline) ./ (baseline + 32769);

% Compute the z-score (with respect to time and space)
zScoreMovie = (zScoreMovie - mean(zScoreMovie(:))) / std(zScoreMovie(:));
end

