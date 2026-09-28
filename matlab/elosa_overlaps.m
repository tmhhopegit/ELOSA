function Overlaps = elosa_overlaps(Regions, Lesions, chunk)
%ELOSA_OVERLAPS Fraction of each region that each patient's lesion covers.
%   Overlaps = elosa_overlaps(Regions, Lesions)
%   Regions: groups x voxels (from elosa_regions)
%   Lesions: patients x voxels (any patients - e.g. all, or a test set)
%   Overlaps: single, groups x patients, in [0, 1]
%   Replaces OverlapLesions (and the local copies in ELOSA_v1a/v2). You only
%   need this for a threshold below 1: with threshold 1, run elosa_run with
%   'Covering', true, which gives the same answer without this matrix.
if nargin < 3, chunk = 5000; end
if size(Regions, 2) ~= size(Lesions, 2)
    error('elosa:size', 'Regions has %d voxels but Lesions has %d.', size(Regions, 2), size(Lesions, 2));
end
L = sparse(double(Lesions ~= 0))';
G = size(Regions, 1);
Overlaps = zeros(G, size(Lesions, 1), 'single');
for c = 1:ceil(G / chunk)
    rows = (c - 1) * chunk + 1:min(c * chunk, G);
    Rg = sparse(double(Regions(rows, :)));
    sizes = full(sum(Rg, 2));
    Overlaps(rows, :) = single(full(Rg * L) ./ sizes);
end
end
