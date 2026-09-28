function T = elosa_growth_curve(Lesions, varargin)
%ELOSA_GROWTH_CURVE How the number of distinct regions grows with the sample.
%   T = elosa_growth_curve(Lesions, 'Sizes', [10 20 40 80], 'Reps', 10)
%   For each sample size and repeat, counts the groups among a random subset
%   of patients (exactly, with elosa). Plot mean(groups) against size to see
%   how fast the information ceiling rises and to extrapolate to sizes too big
%   to enumerate.
%   Options: 'Sizes' (default: 10 steps up to all patients), 'Reps' (10),
%   'Seed' (1), 'MaxSeconds' (stop a size once one run takes longer, default
%   Inf), plus any elosa_run option ('MinVoxels', 'MaxPatients', 'Threads', ...).
%
%   T: table with Size, Rep, Groups, Seconds.
%
%   Replaces MakeEOLSAGrowthCurve, which counted the unique patterns of the
%   ALL-patient groups restricted to each subset. That gives the same counts
%   plus one (the empty pattern), but needs the full enumeration first, so it
%   cannot be used when the full space is too big; this version only needs
%   the subsets.
n = size(Lesions, 1);
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('Sizes', unique(round(linspace(max(2, round(n / 10)), n, 10))));
p.addParameter('Reps', 10);
p.addParameter('Seed', 1);
p.addParameter('MaxSeconds', Inf);
p.parse(varargin{:});
o = p.Results;
pass = namedargs(p.Unmatched);
if ~isfield(p.Unmatched, 'Progress'), pass = [pass, {'Progress', 'none'}]; end   % one line per run instead
stream = RandStream('mt19937ar', 'Seed', o.Seed);
out = zeros(0, 4);
for s = o.Sizes(:)'
    slow = false;
    for r = 1:o.Reps
        idx = randperm(stream, n, s);
        t = tic;
        R = elosa_run(Lesions(idx, :), pass{:});
        secs = toc(t);
        out(end + 1, :) = [s, r, R.count, secs]; %#ok<AGROW>
        fprintf('size %d, repeat %d: %d groups (%.1f s)\n', s, r, R.count, secs);
        slow = slow || secs > o.MaxSeconds;
        if s == n, break; end              % every repeat would be identical
    end
    if slow, fprintf('stopping: runs are taking longer than %g s\n', o.MaxSeconds); break; end
end
T = array2table(out, 'VariableNames', {'Size', 'Rep', 'Groups', 'Seconds'});
end

function c = namedargs(s)
f = fieldnames(s);
c = cell(1, 2 * numel(f));
c(1:2:end) = f;
c(2:2:end) = struct2cell(s);
end
