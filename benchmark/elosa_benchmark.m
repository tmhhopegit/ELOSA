function T = elosa_benchmark(Lesions, sizes, seed)
%ELOSA_BENCHMARK Time the new program against the original MATLAB search.
%   T = elosa_benchmark(Lesions, [20 40 60 80])
%   For each size, a random subset of patients is analysed by
%     new: elosa_run (all cores and 1 thread)
%     old: ExhaustiveRegionFinder.CalcGroupsByVoxel + FindRecursive (the copy
%          in this folder: the original with only its class name fixed)
%   and the number of groups each finds is compared. Stop increasing the sizes
%   once the old code takes minutes - its time grows very quickly.
if nargin < 3, seed = 1; end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'matlab'));
stream = RandStream('mt19937ar', 'Seed', seed);
out = zeros(0, 7);
for s = sizes(:)'
    L = logical(Lesions(randperm(stream, size(Lesions, 1), s), :));
    t = tic; R = elosa_run(L); t_new = toc(t);
    t = tic; R1 = elosa_run(L, 'Threads', 1); t_new1 = toc(t); %#ok<NASGU>
    t = tic;
    E = ExhaustiveRegionFinder(L);
    E = E.CalcGroupsByVoxel();
    E = E.FindRecursive(1);
    t_old = toc(t);
    old_groups = size(E.Groups, 1);
    % how many of the old groups are real groups (non-empty region, closed)?
    real = 0;
    for g = 1:old_groups
        members = E.Groups(g, :);
        region = all(L(members, :), 1);
        real = real + (any(region) && isequal(all(L(:, region), 2)', members));
    end
    out(end + 1, :) = [s, R.count, old_groups, real, t_new1, t_new, t_old]; %#ok<AGROW>
    fprintf('%d patients: new %d groups in %.2f s (1 thread) / %.2f s; old %d groups (%d valid) in %.1f s\n', ...
        s, R.count, t_new1, t_new, old_groups, real, t_old);
end
T = array2table(out, 'VariableNames', {'Patients', 'Groups', 'OldGroups', 'OldValid', ...
    'Seconds_1thread', 'Seconds', 'OldSeconds'});
end
