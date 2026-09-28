function Regions = elosa_regions(Groups, Lesions, chunk)
%ELOSA_REGIONS Voxels lesioned in every patient of each group.
%   Regions = elosa_regions(Groups, Lesions)
%   Groups:  groups x patients logical (from elosa_run / elosa_read_groups)
%   Lesions: patients x voxels, the SAME patients (columns of Groups)
%   Regions: sparse logical, groups x voxels
%   Groups are processed in chunks (default 5000) with sparse products, so
%   this is fast even for many groups; the original built each region with
%   prod(Lesions(group,:)) one group at a time.
if nargin < 3, chunk = 5000; end
if size(Groups, 2) ~= size(Lesions, 1)
    error('elosa:size', 'Groups has %d columns but Lesions has %d patients.', size(Groups, 2), size(Lesions, 1));
end
L = sparse(double(Lesions ~= 0));
k = full(sum(Groups, 2));
G = size(Groups, 1);
parts = cell(ceil(G / chunk), 1);
for c = 1:numel(parts)
    rows = (c - 1) * chunk + 1:min(c * chunk, G);
    counts = sparse(double(Groups(rows, :))) * L;          % lesioned group members per voxel
    [i, j, v] = find(counts);
    full_overlap = v == k(rows(i));
    parts{c} = sparse(i(full_overlap), j(full_overlap), true, numel(rows), size(L, 2));
end
if isempty(parts)
    Regions = sparse(false(0, size(L, 2)));
else
    Regions = vertcat(parts{:});
end
end
