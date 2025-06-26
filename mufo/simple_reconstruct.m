function volume = simple_reconstruct(sinogram, angles, center)
%SIMPLE_RECONSTRUCT Reconstruct volume from sinograms using FBP.
%
%   VOLUME = SIMPLE_RECONSTRUCT(SINOGRAM, ANGLES, CENTER) performs filtered
%   backprojection of the input SINOGRAM stack. ANGLES are in radians and
%   CENTER specifies the rotation axis position in pixels. The output is a
%   3-D volume with size [N N numSlices], where N is the detector width.
%
%   This helper relies on MATLAB's IRADON implementation and is intended
%   for small examples.

    narginchk(2,3);
    if nargin < 3 || isempty(center)
        center = size(sinogram,1)/2;
    end

    angleDeg = angles * 180 / pi;
    numSlices = size(sinogram,3);
    N = size(sinogram,1);
    volume = zeros(N, N, numSlices, 'like', sinogram);

    for k = 1:numSlices
        slice = sinogram(:,:,k);
        shift = round(center - N/2);
        slice = circshift(slice, [shift, 0]);
        volume(:,:,k) = iradon(slice, angleDeg, 'linear', 'Ram-Lak', 1, N);
    end
end
