function [corrected, sinogram, center] = simple_preprocess(proj, flat, dark, angles)
%SIMPLE_PREPROCESS  Basic preprocessing helper for MUFO examples.
%
%   [corrected, sinogram, center] = SIMPLE_PREPROCESS(proj, flat, dark, angles)
%   performs a simple flat-field correction on the input projection stack,
%   converts the result to sinogram orientation and estimates the center of
%   rotation from the middle slice.  ``proj``, ``flat`` and ``dark`` must be
%   numeric arrays with matching detector dimensions and ``angles`` provides
%   the projection angles in radians.
%
%   The implementation mirrors the logic of
%   ``mufo.tasks.FlatFieldCorrect.previewCorrection`` and uses
%   ``mufo.tasks.CenterOfRotation.detectFromSinogram`` for the center
%   estimation.

    narginchk(4,4);

    % verify array sizes
    if ndims(proj) ~= 3
        error('simple_preprocess:InvalidInput', 'Projections must be a 3-D array.');
    end

    if ndims(flat) == 3
        flatMean = mean(flat, 3);
    else
        flatMean = flat;
    end

    if ndims(dark) == 3
        darkMean = mean(dark, 3);
    else
        darkMean = dark;
    end

    if ~isequal(size(flatMean), size(darkMean))
        error('simple_preprocess:SizeMismatch', 'Flat and dark fields must have the same size.');
    end

    if any(size(proj,1:2) ~= size(flatMean))
        error('simple_preprocess:SizeMismatch', 'Projection and reference dimensions must agree.');
    end

    numAngles = numel(angles);
    dims = size(proj);
    if ~any(dims == numAngles)
        warning('simple_preprocess:AngleCount', ...
                'Number of angles does not match data dimensions. Proceeding anyway.');
    end

    % --- Flat-field correction (absorption mode) ---
    darkScale = 1.0;
    flatScale = 1.0;

    dRef = darkScale * darkMean;
    fRef = flatScale * flatMean - dRef;

    % replicate reference frames over all projections
    dRefRep = repmat(dRef, 1, 1, size(proj,3));
    fRefRep = repmat(fRef, 1, 1, size(proj,3));

    ratio = (proj - dRefRep) ./ fRefRep;
    corrected = -log(ratio);

    % replace invalid values with zero
    corrected(isnan(corrected)) = 0;
    corrected(isinf(corrected)) = 0;

    % --- Generate sinograms ---
    if dims(3) == numAngles
        % data ordered as [height,width,angles]
        sinogram = permute(corrected, [2 3 1]); % -> [width,angles,height]
    elseif dims(2) == numAngles
        % data ordered as [width,angles,height]
        sinogram = corrected;                   % already oriented
    elseif dims(1) == numAngles
        % data ordered as [angles,width,height]
        sinogram = permute(corrected, [2 1 3]); % -> [width,angles,height]
    else
        % default assumption
        sinogram = permute(corrected, [2 3 1]);
    end

    % --- Detect center of rotation from mid-slice ---
    midSlice = round(size(sinogram,3)/2);
    corTask = mufo.tasks.CenterOfRotation();
    angleStep = mean(diff(angles));
    center = corTask.detectFromSinogram(sinogram(:,:,midSlice), angleStep);
end
