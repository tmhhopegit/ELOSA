function [Solution, Summary] = elosa_greedy_cls(Covered, Impaired, varargin)
%ELOSA_GREEDY_CLS Greedy search for critical lesion sites among ELOSA regions.
%   [Solution, Summary] = elosa_greedy_cls(Covered, Impaired, ...)
%   Covered:  groups x patients logical: patient's lesion covers the group's
%             region (R.Covering from elosa_run(..., 'Covering', true), or
%             elosa_overlaps(...) >= threshold)
%   Impaired: patients x 1 logical
%
%   Repeatedly picks the region that best explains the impaired patients not
%   yet explained, then removes the impaired patients it covers.
%   Options:
%     'Rule'          'min_ppv' (default; was GetCLS): among regions with
%                     PPV >= MinPPV, take the one covering most unexplained
%                     impaired patients; stop when no region reaches MinPPV.
%                     'weighted' (was FindGreedyCLS): maximise PPV^W * SEN;
%                     stop when no region covers an unexplained impaired patient.
%     'MinPPV'        default 0.9
%     'W'             default 1
%     'Train'         logical, patients x 1: choose regions on these patients
%                     only, and report performance on the rest (default: all)
%     'RegionVoxels'  groups x 1 region sizes, for tie-breaks. Ties go to the
%                     LARGEST region for 'min_ppv' (the original weighted size
%                     by the unimpaired patients' overlap margin; pass
%                     'Overlaps' to get exactly that) and to the SMALLEST
%                     region for 'weighted', as the originals did.
%     'Overlaps'      groups x patients overlap fractions (elosa_overlaps), for
%                     the original 'min_ppv' tie-break: size x (Threshold - max
%                     overlap among unimpaired patients)
%     'Threshold'     the coverage threshold used to build Covered (default 1;
%                     only used with 'Overlaps')
%
%   Solution: table, one row per chosen region:
%     Group         row of Covered
%     Right         unexplained impaired training patients it covers (now explained)
%     Wrong         unimpaired training patients it covers
%     PPV, SEN      at the time it was chosen (on the unexplained patients)
%     Solution_PPV, Solution_SEN   of all regions chosen so far, on all training patients
%     Test_*        the same on held-out patients, if 'Train' leaves any out
%   Summary: struct with the number of impaired training patients left unexplained.
%
%   Fixes relative to GetCLS/FindGreedyCLS: ELOSA_v4 passed Imp as indices into
%   one subset while indexing another; FindGreedyCLS looped forever when no
%   region covered the remaining patients, and divided held-out sensitivity by
%   the number of held-out patients instead of held-out impaired patients.
p = inputParser;
p.addParameter('Rule', 'min_ppv');
p.addParameter('MinPPV', 0.9);
p.addParameter('W', 1);
p.addParameter('Train', []);
p.addParameter('RegionVoxels', []);
p.addParameter('Overlaps', []);
p.addParameter('Threshold', 1);
p.parse(varargin{:});
o = p.Results;
Covered = logical(Covered);
Impaired = logical(Impaired(:));
n = size(Covered, 2);
if numel(Impaired) ~= n, error('elosa:size', 'Impaired must have one entry per column of Covered.'); end
Train = o.Train;
if isempty(Train), Train = true(n, 1); end
Train = logical(Train(:));
Test = ~Train;
G = size(Covered, 1);
sizes = o.RegionVoxels(:);
if isempty(sizes), sizes = zeros(G, 1); end

active = Train;                       % training patients not yet explained
chosen = false(G, 1);
solution_cover = false(1, n);
first_score = [];
rows = {};
while true
    imp_left = active & Impaired;
    n_imp = nnz(imp_left);
    if n_imp == 0, break; end
    right = full(sum(Covered(:, imp_left), 2));
    total = full(sum(Covered(:, active), 2));
    ppv = right ./ total;              % NaN where a region covers nobody still active
    sen = right / n_imp;
    switch o.Rule
        case 'min_ppv'
            ok = ppv >= o.MinPPV & right > 0;
            if ~any(ok), break; end
            best_sen = max(sen(ok));
            cand = find(ok & sen == best_sen);
            if numel(cand) > 1
                if ~isempty(o.Overlaps)
                    unimp = ~Impaired & Train;
                    margin = o.Threshold - max(double(o.Overlaps(cand, unimp)), [], 2);
                    score = sizes(cand) .* margin;
                else
                    score = sizes(cand);
                end
                cand = cand(score == max(score));
            end
        case 'weighted'
            q = ppv .^ o.W .* sen;
            q(isnan(q)) = 0;
            if isempty(first_score), first_score = q; end
            if max(q) <= 0, break; end   % nothing covers the unexplained patients
            cand = find(q == max(q));
            if numel(cand) > 1
                f = first_score(cand);
                cand = cand(f == max(f));
            end
            if numel(cand) > 1
                cand = cand(sizes(cand) == min(sizes(cand)));
            end
        otherwise
            error('elosa:rule', 'Rule must be ''min_ppv'' or ''weighted''.');
    end
    g = cand(1);
    chosen(g) = true;
    cov = full(Covered(g, :))';
    explained = cov & imp_left;
    solution_cover = solution_cover | cov';
    row = struct('Group', g, 'Right', nnz(explained), 'Wrong', nnz(cov & active & ~Impaired), ...
        'PPV', ppv(g), 'SEN', sen(g), ...
        'Solution_PPV', nnz(solution_cover' & Train & Impaired) / max(nnz(solution_cover' & Train), 1), ...
        'Solution_SEN', nnz(solution_cover' & Train & Impaired) / max(nnz(Train & Impaired), 1));
    if any(Test)
        row.Test_Right = nnz(cov & Test & Impaired);
        row.Test_Wrong = nnz(cov & Test & ~Impaired);
        row.Test_PPV = row.Test_Right / max(row.Test_Right + row.Test_Wrong, 1);
        row.Test_SEN = row.Test_Right / max(nnz(Test & Impaired), 1);
        row.Solution_Test_PPV = nnz(solution_cover' & Test & Impaired) / max(nnz(solution_cover' & Test), 1);
        row.Solution_Test_SEN = nnz(solution_cover' & Test & Impaired) / max(nnz(Test & Impaired), 1);
    end
    rows{end + 1} = row; %#ok<AGROW>
    active = active & ~explained;
end
if isempty(rows)
    Solution = table();
else
    Solution = struct2table([rows{:}]');
end
Summary.unexplained = nnz(active & Impaired);
Summary.impaired_train = nnz(Train & Impaired);
end
