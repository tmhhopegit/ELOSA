function elosa_write_input(file, Lesions, Use)
%ELOSA_WRITE_INPUT Write lesions in the format the elosa program reads.
%   elosa_write_input(file, Lesions)       every patient can be grouped
%   elosa_write_input(file, Lesions, Use)  only patients with Use true are
%       grouped; the others are only counted when their lesion covers a
%       region (e.g. Use = impaired patients, the rest give the PPV).
%   Lesions: patients x voxels, logical (or 0/1; anything > 0 counts as lesioned).
n = size(Lesions, 1);
if nargin < 3 || isempty(Use), Use = true(n, 1); end
if numel(Use) ~= n, error('elosa:use', 'Use must have one entry per patient (row of Lesions).'); end
if n >= 2^32 || size(Lesions, 2) >= 2^32, error('elosa:size', 'Too many patients or voxels.'); end
fid = fopen(file, 'w');
if fid < 0, error('elosa:write', 'Cannot write %s.', file); end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, 'ELOSAIN1', 'char');
fwrite(fid, [n, size(Lesions, 2)], 'uint32');
fwrite(fid, uint8(Use(:) ~= 0), 'uint8');
% MATLAB is column-major, so an n x V matrix is already stored voxel by voxel;
% write in blocks of voxels to limit memory for large inputs
block = max(1, floor(5e7 / max(n, 1)));
for v = 1:block:size(Lesions, 2)
    cols = v:min(v + block - 1, size(Lesions, 2));
    fwrite(fid, uint8(full(Lesions(:, cols)) > 0), 'uint8');
end
end
