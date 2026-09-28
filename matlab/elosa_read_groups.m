function [Groups, RegionVoxels, NumOther, Covering] = elosa_read_groups(prefix)
%ELOSA_READ_GROUPS Read the groups written by "elosa ... --out PREFIX".
%   [Groups, RegionVoxels, NumOther, Covering] = elosa_read_groups(prefix)
%   Groups        sparse logical, groups x patients: the patients in each group
%                 (every grouped patient whose lesion covers the whole region)
%   RegionVoxels  voxels in each group's region
%   NumOther      number of non-grouped patients whose lesions cover the region
%   Covering      sparse logical, groups x patients: group members plus those
%                 other patients (only if elosa was run with --covering;
%                 otherwise [])
fid = fopen([prefix '.groups'], 'r');
if fid < 0, error('elosa:read', 'Cannot read %s.groups.', prefix); end
magic = fread(fid, [1 8], '*char');
if ~strcmp(magic, 'ELOSAGR2'), fclose(fid); error('elosa:read', '%s.groups is not an ELOSA groups file.', prefix); end
n = fread(fid, 1, 'uint32');
flags = fread(fid, 1, 'uint32');
H = fread(fid, [4 inf], 'uint32=>double');
fclose(fid);
k = H(1, :)';
NumOther = H(2, :)';
RegionVoxels = H(3, :)' + H(4, :)' * 2^32;
G = numel(k);

fid = fopen([prefix '.members'], 'r');
if fid < 0, error('elosa:read', 'Cannot read %s.members.', prefix); end
if bitand(flags, 2), prec = 'uint16=>double'; else, prec = 'uint32=>double'; end
M = fread(fid, inf, prec);
fclose(fid);

has_cov = bitand(flags, 1) ~= 0;
if has_cov, len = k + NumOther; else, len = k; end
if numel(M) ~= sum(len), error('elosa:read', '%s.members does not match %s.groups.', prefix, prefix); end
gid = repelem((1:G)', len);
if has_cov
    pos = (1:numel(M))' - repelem(cumsum([0; len(1:end - 1)]), len);   % position within its record
    in_group = pos <= repelem(k, len);
    Groups = sparse(gid(in_group), M(in_group), true, G, n);
    Covering = sparse(gid, M, true, G, n);
else
    Groups = sparse(gid, M, true, G, n);
    Covering = [];
end
end
