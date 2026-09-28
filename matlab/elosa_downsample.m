function [LesionsNew, BrainVoxNew, DimsNew] = elosa_downsample(Lesions, BrainVox, Dims, Factor, MinFraction)
%ELOSA_DOWNSAMPLE Coarser voxels, so the overlap space gets smaller.
%   [LesionsNew, BrainVoxNew, DimsNew] = elosa_downsample(Lesions, BrainVox, Dims, Factor)
%   Lesions:  patients x brain voxels (as LesionsBinary1D)
%   BrainVox: linear indices of those voxels in a volume of size Dims
%             (e.g. [91 109 91] for 2 mm MNI)
%   Factor:   integer block size, e.g. 2 turns 2 mm voxels into 4 mm voxels
%   MinFraction: a new voxel is lesioned if at least this fraction of its
%             brain voxels are (default 0.5)
%   Returns binary lesions over the new brain voxels (blocks containing any
%   brain voxel), their indices, and the new volume size.
%
%   Pure MATLAB; replaces ReduceResolution / ChangeLesionResolutios, which
%   needed SPM, SaveVolume and write access to a hard-coded .\ReduceResolution\
%   folder, flipped the x axis after resampling, and kept the interpolated
%   (non-binary) values, so any partial overlap counted as a lesion.
if nargin < 5, MinFraction = 0.5; end
Dims = Dims(:)';
DimsNew = ceil(Dims / Factor);
[x, y, z] = ind2sub(Dims, BrainVox(:));
block = sub2ind(DimsNew, ceil(x / Factor), ceil(y / Factor), ceil(z / Factor));
[BrainVoxNew, ~, b] = unique(block);
nb = numel(BrainVoxNew);
% brain voxels per block, and lesioned voxels per block for every patient
M = sparse(1:numel(b), b, 1, numel(b), nb);
brain_count = full(sum(M, 1));
lesion_count = double(Lesions ~= 0) * M;          % patients x blocks
LesionsNew = full(lesion_count ./ brain_count >= MinFraction);
end
