function ROIs = intersectingROIs(mask, image, threshold)
%INTERSECTINGROIS Get indices of ROIs touching points above the threshold.

ROIs = unique(mask .* (image >= threshold));

% Remove zero from the list (zero means that intersects the background).
ROIs = ROIs(ROIs > 0);

end