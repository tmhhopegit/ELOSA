classdef ElosaTest < matlab.unittest.TestCase
%ELOSATEST Checks the compiled program and the MATLAB wrappers.
%   From the refactored folder:  addpath matlab tests; runtests('ElosaTest')
%   Needs bin/elosa(.exe) (run elosa_compile first).

    methods (TestClassSetup)
        function paths(tc) %#ok<MANU>
            root = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(root, 'matlab'));
        end
    end

    methods (Test)
        function testMatchesBruteForce(tc)
            s = RandStream('mt19937ar', 'Seed', 1);
            for trial = 1:60
                n = randi(s, [2 9]);
                L = rand(s, n, randi(s, [5 60])) < 0.2 + 0.5 * rand(s);
                use = true(n, 1);
                if mod(trial, 3) == 0, use = rand(s, n, 1) < 0.7; use(1) = true; end
                minv = 1 + (mod(trial, 4) == 0) * 2;
                [G, sup, other] = ElosaTest.brute(L, use, minv);
                R = elosa_run(L, 'Use', use, 'Covering', true, 'MinVoxels', minv, 'Threads', 2);
                tc.verifyEqual(R.count, size(G, 1));
                [A, ia] = sortrows(full(R.Groups));
                [B, ib] = sortrows(G);
                tc.verifyEqual(A, B);
                tc.verifyEqual(R.RegionVoxels(ia), sup(ib));
                tc.verifyEqual(R.NumOther(ia), other(ib));
                % Covering = group members + covering non-grouped patients
                C = full(R.Covering(ia, :));
                tc.verifyEqual(C(:, use), A(:, use));
                tc.verifyEqual(sum(C(:, ~use), 2), other(ib));
            end
        end

        function testRegionsAndOverlaps(tc)
            s = RandStream('mt19937ar', 'Seed', 2);
            L = rand(s, 8, 50) < 0.4;
            R = elosa_run(L, 'Groups', true);
            Reg = elosa_regions(R.Groups, L);
            for g = 1:size(R.Groups, 1)
                tc.verifyEqual(full(Reg(g, :)), all(L(R.Groups(g, :), :), 1));
            end
            tc.verifyEqual(full(sum(Reg, 2)), R.RegionVoxels);
            O = elosa_overlaps(Reg, L);
            % a patient covers a region completely exactly when it is in the group
            tc.verifyEqual(O >= 1, full(R.Groups));
        end

        function testEstimateIsPlausible(tc)
            s = RandStream('mt19937ar', 'Seed', 3);
            L = rand(s, 12, 200) < 0.3;
            R = elosa_run(L);
            E = elosa_run(L, 'Estimate', 20000);
            tc.verifyLessThan(abs(E.estimate - R.count), max(5 * E.se, 0.2 * R.count));
        end

        function testGreedyMinPPV(tc)
            % patients 1-4 impaired, 5-6 not. Region 1 covers 1,2,5 (PPV .67);
            % region 2 covers 1,2,3 (PPV 1); region 3 covers 4 (PPV 1).
            Covered = logical([1 1 0 0 1 0; 1 1 1 0 0 0; 0 0 0 1 0 0]);
            Impaired = logical([1 1 1 1 0 0])';
            [S, sum_] = elosa_greedy_cls(Covered, Impaired, 'MinPPV', 0.9);
            tc.verifyEqual(S.Group, [2; 3]);
            tc.verifyEqual(S.Right, [3; 1]);
            tc.verifyEqual(sum_.unexplained, 0);
            tc.verifyEqual(S.Solution_SEN(end), 1);
        end

        function testGreedyWeightedStops(tc)
            % nothing covers patient 3: the original looped forever here
            Covered = logical([1 1 0 0; 1 0 0 1]);
            Impaired = logical([1 1 1 0])';
            [S, sum_] = elosa_greedy_cls(Covered, Impaired, 'Rule', 'weighted', 'W', 1);
            tc.verifyEqual(S.Group, 1);
            tc.verifyEqual(sum_.unexplained, 1);
        end

        function testGreedyTrainTest(tc)
            Covered = logical([1 1 0 1 1 0]);
            Impaired = logical([1 1 0 1 0 1])';
            Train = logical([1 1 1 0 0 0])';
            S = elosa_greedy_cls(Covered, Impaired, 'Train', Train, 'MinPPV', 0.5);
            tc.verifyEqual(S.Test_Right, 1);
            tc.verifyEqual(S.Test_Wrong, 1);
            tc.verifyEqual(S.Test_SEN, 1 / 2);    % 2 impaired test patients, not 3 test patients
        end

        function testDownsample(tc)
            dims = [4 4 2];
            vol = false(dims); vol(1:2, 1:2, 1) = true; vol(3, 3, 2) = true;
            brain = find(true(dims));
            L = [vol(brain)'; false(1, numel(brain))];
            [Ln, bv, dn] = elosa_downsample(L, brain, dims, 2);
            tc.verifyEqual(dn, [2 2 1]);
            tc.verifyEqual(numel(bv), 4);
            tc.verifyEqual(Ln(1, :), logical([1 0 0 0]));   % one of 8 voxels is below 0.5
            tc.verifyFalse(any(Ln(2, :)));
        end

        function testSelectPatients(tc)
            D.LesionsBinary1D = logical([1 0; 0 0; 1 1; 1 0]);
            D.Behaviours = [1; 2; NaN; 4];
            D.ImpairmentLabels = false(4, 7);
            D.TimePost = [3; 5; 10; 20];
            D.Order = [1; 1; 1; 2];
            D.Hands = {'Right'; 'Left'; 'Right'; 'Now Right'};
            D.L1 = {'English'; 'ENGLISH'; 'French'; 'English'};
            D.VolLeft = [1; 0; 1; 1];
            D.VolRight = [0; 0; 0; 0];
            S = elosa_select_patients(D, {'L1==English', 'MinTime==4'}, 1);
            tc.verifyEqual(S.rows, logical([0 0 0 1])');   % 1 too early, 2 empty lesion, 3 French/NaN
            S = elosa_select_patients(D, {'RightHanded_Only', 'FirstScanOnly'});
            tc.verifyEqual(S.rows, logical([1 0 1 0])');
            tc.verifyError(@() elosa_select_patients(D, {'Typo'}), 'elosa:criterion');
        end
    end

    methods (Static)
        function [G, sup, other] = brute(L, use, minv)
            % every group by definition: closures of intersections of voxel patterns
            L = logical(L); use = logical(use(:));
            cols = L(:, any(L(use, :), 1))';                 % voxels involving a grouped patient
            [atoms, ~, ic] = unique(cols, 'rows');
            w = accumarray(ic, 1);
            fam = atoms(:, use);
            frontier = fam;
            while ~isempty(frontier)
                new = false(0, nnz(use));
                for i = 1:size(frontier, 1)
                    new = [new; frontier(i, :) & fam]; %#ok<AGROW>
                end
                new = unique(new(any(new, 2), :), 'rows');
                new = new(~ismember(new, fam, 'rows'), :);
                fam = [fam; new]; %#ok<AGROW>
                frontier = new;
            end
            fam = unique(fam(any(fam, 2), :), 'rows');
            n = size(L, 1);
            G = false(0, n); sup = zeros(0, 1); other = zeros(0, 1);
            for i = 1:size(fam, 1)
                g = false(1, n); g(use) = fam(i, :);
                inside = all(atoms(:, g), 2);
                s = sum(w(inside));
                if s < minv, continue; end
                cover = all(atoms(inside, :), 1);
                G(end + 1, :) = g; %#ok<AGROW>
                sup(end + 1, 1) = s; %#ok<AGROW>
                other(end + 1, 1) = nnz(cover & ~use'); %#ok<AGROW>
            end
        end
    end
end
