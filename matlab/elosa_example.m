%ELOSA_EXAMPLE A complete ELOSA analysis, section by section.
% Needs: the compiled program (run elosa_compile once) and a data object D
% (CLS / PLORAS_parallel object, or a load_ploras struct) in the workspace.
addpath(fileparts(mfilename('fullpath')));

%% 1. Choose patients and the impairment to explain
S = elosa_select_patients(D, {'L1==English', 'RightHanded_Only', 'FirstScanOnly'});
label = 3;                                   % column of ImpairmentLabels
Impaired = logical(S.Labels(:, label));
Lesions = S.Lesions;

%% 2. How big is the space? (fast; nothing stored)
R = elosa_run(Lesions, 'Use', Impaired);
fprintf('%d distinct overlap regions among %d impaired patients\n', R.count, nnz(Impaired));
disp(R.by_size)

% If that takes too long: estimate it, or restrict it
% E = elosa_run(Lesions, 'Use', Impaired, 'Estimate', 20000);
% R = elosa_run(Lesions, 'Use', Impaired, 'MinVoxels', floor(800 / 8));   % >= 800 mm^3 at 2 mm

%% 3. Growth curve: how fast the ceiling rises with sample size
T = elosa_growth_curve(Lesions(Impaired, :), 'Reps', 5, 'MaxSeconds', 60);
G = groupsummary(T, 'Size', {'mean', 'std'}, 'Groups');
errorbar(G.Size, G.mean_Groups, G.std_Groups); set(gca, 'YScale', 'log');
xlabel('patients'); ylabel('distinct overlap regions');

%% 4. Critical lesion sites: regions whose full coverage predicts impairment
R = elosa_run(Lesions, 'Use', Impaired, 'Covering', true);
[Solution, Summary] = elosa_greedy_cls(R.Covering, Impaired, 'MinPPV', 0.9, ...
    'RegionVoxels', R.RegionVoxels);
disp(Solution)
fprintf('%d of %d impaired patients unexplained\n', Summary.unexplained, Summary.impaired_train);

% The chosen regions as voxel masks:
Regions = elosa_regions(R.Groups(Solution.Group, :), Lesions);

%% 5. With a train/test split (e.g. later scans train, earlier scans test)
Train = S.Time > 12;
R = elosa_run(Lesions, 'Use', Impaired & Train, 'Covering', true);
Solution = elosa_greedy_cls(R.Covering, Impaired, 'Rule', 'weighted', 'W', 2, ...
    'Train', Train, 'RegionVoxels', R.RegionVoxels);
disp(Solution)
