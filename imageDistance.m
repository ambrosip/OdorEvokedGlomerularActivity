function imageDist = imageDistance(img1, img2)
% IMAGEDISTANCE Computes the normalized euclidean distance between two images.
% The normalization makes the distance independent of the image size.

if isstruct(img1) && isstruct(img2) && ...
        isfield(img1, 'image') && isfield(img2, 'image')
    if size(img1.image) ~= size(img2.image)
        error('Images must be the same size.');
    end
    
    nPixels = numel(img1.image);
    imageDist = norm(img1.image - img2.image, 'fro') / sqrt(nPixels);
elseif ismatrix(img1) && ismatrix(img2)
    if size(img1) ~= size(img2)
        error('Images must be the same size.');
    end

    nPixels = numel(img1);
    imageDist = norm(img1 - img2, 'fro') / sqrt(nPixels);
else
    error('Inputs must be both outcome structs or images.');
end

end