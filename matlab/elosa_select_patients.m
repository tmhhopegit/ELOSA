function S = elosa_select_patients(D, criteria, behaviours, varargin)
%ELOSA_SELECT_PATIENTS Choose the patients for an ELOSA analysis (was GetGoodData).
%   S = elosa_select_patients(D, criteria, behaviours)
%   D:          a CLS object, a PLORAS_parallel object, or a struct from the
%               Confidence project's load_ploras (fields LesionsBinary1D,
%               Behaviours, ImpairmentLabels, TimePost, Order, Hands, L1,
%               VolLeft, VolRight)
%   criteria:   cell array, any of
%                 'L1==English'       first language English
%                 'MinTime==x'        time post-stroke >= x (months)
%                 'MaxTime==x'        time post-stroke <= x
%                 'LH_Only'           left-hemisphere lesions only
%                 'RightHanded_Only'
%                 'GoodCog'           unimpaired on the "cognitive" label columns
%                                     (option 'GoodCogColumns', default [1 2 4 5 6 7]
%                                     as in the original - check these for your data)
%                 'FirstScanOnly'     Order == 1
%               Unknown criteria are an error (the original silently ignored them).
%   behaviours: score columns that must not be missing (default: none). The
%               original referred to an undefined NoNans here and crashed.
%
%   S.rows (logical, patients of D kept), S.Lesions, S.Labels, S.Behaviours,
%   S.Time. Patients with empty lesions are always dropped.
if nargin < 2, criteria = {}; end
if nargin < 3, behaviours = []; end
p = inputParser;
p.addParameter('GoodCogColumns', [1 2 4 5 6 7]);
p.parse(varargin{:});
if ischar(criteria), criteria = {criteria}; end

Lesions = D.LesionsBinary1D;
n = size(Lesions, 1);
keep = full(sum(Lesions ~= 0, 2)) > 0;
for b = behaviours(:)'
    keep = keep & ~isnan(D.Behaviours(:, b));
end
for i = 1:numel(criteria)
    c = criteria{i};
    if strcmp(c, 'L1==English')
        keep = keep & text_has(D.L1, n, {'nglish', 'NGLISH'});
    elseif startsWith(c, 'MinTime==')
        keep = keep & D.TimePost(:) >= str2double(c(10:end));
    elseif startsWith(c, 'MaxTime==')
        keep = keep & D.TimePost(:) <= str2double(c(10:end));
    elseif strcmp(c, 'LH_Only')
        keep = keep & D.VolLeft(:) > 0 & D.VolRight(:) == 0;
    elseif strcmp(c, 'RightHanded_Only')
        % as ConvertDemographicsToPredictors: Left, Ambi and "...ow Right" are not right-handed
        hands = cellstr(string(D.Hands(:)));
        nonright = strcmp(hands, 'Left') | strcmp(hands, 'Ambi') | contains(hands, 'ow Right');
        keep = keep & ~nonright;
    elseif strcmp(c, 'GoodCog')
        keep = keep & sum(D.ImpairmentLabels(:, p.Results.GoodCogColumns), 2) == 0;
    elseif strcmp(c, 'FirstScanOnly')
        keep = keep & D.Order(:) == 1;
    else
        error('elosa:criterion', 'Unknown inclusion criterion "%s".', c);
    end
end
S.rows = keep;
S.Lesions = Lesions(keep, :);
S.Labels = D.ImpairmentLabels(keep, :);
S.Behaviours = D.Behaviours(keep, :);
S.Time = D.TimePost(keep);
end

function tf = text_has(values, n, patterns)
values = cellstr(string(values(:)));
if numel(values) ~= n, error('elosa:field', 'Text field has %d entries for %d patients.', numel(values), n); end
tf = contains(values, patterns);
end
