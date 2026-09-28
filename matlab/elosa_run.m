function R = elosa_run(Lesions, varargin)
%ELOSA_RUN Find every distinct lesion-overlap region (group) in a dataset.
%   R = elosa_run(Lesions)            count the groups (fast; nothing stored)
%   R = elosa_run(Lesions, 'Groups', true)   also return every group
%
%   A group is a set of patients G whose lesions all overlap somewhere, such
%   that no other patient's lesion covers all of that overlap region. Each
%   group is one distinct region the dataset can tell apart, so the number of
%   groups measures the information ceiling of the lesion data.
%
%   Lesions: patients x voxels (logical or 0/1).
%   Options (name/value):
%     'Use'          logical, patients x 1: only these patients form groups
%                    (e.g. the impaired). Others count only as "covering"
%                    patients, for PPV. Default: all patients.
%     'Groups'       true: read every group back (can be very large). Default false.
%     'Covering'     true: also return which non-grouped patients cover each
%                    region (R.Covering). Implies 'Groups'. Default false.
%     'MinVoxels'    smallest region to count (default 1 = exhaustive). The old
%                    code meant to use floor(800 / prod(VoxSize)) but never did.
%     'MaxPatients'  largest group to count (default 0 = no limit).
%     'Threads'      default 0 = all cores.
%     'Estimate'     N > 0: don't enumerate; estimate the number of groups from
%                    N random walks (Knuth's estimator). Useful when the space
%                    is too big to enumerate; in tests it UNDERestimated large
%                    spaces (e.g. 1.6-2.1 million for a true 3.3 million), so
%                    treat it as an order of magnitude.
%     'Seed'         seed for 'Estimate' (default 1).
%     'Keep'         a file prefix: keep the input/output files there instead of
%                    deleting them (default: temporary files, deleted).
%     'Exe'          path to the elosa program (default: elosa_exe()).
%     'Progress'     'text' (default): a status line in the Command Window,
%                    updated every half second, for runs longer than a second;
%                    'bar': a progress bar with a Cancel button; 'none'.
%                    It shows groups found and time elapsed (exact) and a rough
%                    percentage and time left (see elosa_progress_display).
%                    Ctrl+C (or Cancel) stops the program too, and it stops
%                    by itself if MATLAB is closed or crashes. elosa_stop
%                    ends any copies still running from older versions.
%
%   R has fields:
%     count         number of groups (NaN in estimate mode)
%     estimate, se  estimated number of groups and its standard error (estimate mode)
%     summary       struct: patients, grouped_patients, lesioned_voxels,
%                   distinct_voxel_patterns, seconds, ...
%     by_size       table: groups by number of patients
%     by_region     table: groups by region size (voxel ranges)
%     Groups, RegionVoxels, NumOther, Covering   see elosa_read_groups
p = inputParser;
p.addParameter('Use', []);
p.addParameter('Groups', false);
p.addParameter('Covering', false);
p.addParameter('MinVoxels', 1);
p.addParameter('MaxPatients', 0);
p.addParameter('Threads', 0);
p.addParameter('Estimate', 0);
p.addParameter('Seed', 1);
p.addParameter('Keep', '');
p.addParameter('Exe', '');
p.addParameter('Progress', 'text');
p.parse(varargin{:});
o = p.Results;
if isempty(o.Exe), o.Exe = elosa_exe(); end
want_groups = o.Groups || o.Covering;

if isempty(o.Keep)
    prefix = tempname;
    cleanup = onCleanup(@() delete_files(prefix));
else
    prefix = o.Keep;
end
infile = [prefix '.in'];
sumfile = [prefix '.summary.txt'];
progfile = [prefix '.progress.txt'];
if ~strcmp(o.Progress, 'none') && numel(Lesions) > 5e7
    fprintf('ELOSA: saving the lesions for the program...\n');
end
elosa_write_input(infile, Lesions, o.Use);

args = {infile, '--summary', sumfile, '--progress', progfile, '--parent-pid', matlab_pid(), '--threads', num2str(o.Threads), ...
    '--min-voxels', num2str(o.MinVoxels), '--max-patients', num2str(o.MaxPatients)};
if o.Estimate > 0
    args = [args, {'--estimate', num2str(o.Estimate), '--seed', num2str(o.Seed)}];
elseif want_groups
    args = [args, {'--out', prefix}];
    if o.Covering, args{end + 1} = '--covering'; end
end
[poll, finish] = elosa_progress_display(progfile, o.Progress);
clear_display = onCleanup(finish);        % also on errors and Ctrl+C
started = tic;
[status, out, cmd] = elosa_system(o.Exe, args, poll);
finish();
if exist(progfile, 'file'), delete(progfile); end
if status == -2
    error('elosa:cancelled', 'Cancelled.');
end
if status ~= 0
    error('elosa:run', 'elosa failed (%d):\n%s', status, out);
end

[R.summary, R.by_size, R.by_region] = read_summary(sumfile);
if o.Estimate > 0
    R.count = NaN;
    R.estimate = R.summary.estimated_groups;
    R.se = R.summary.standard_error;
else
    R.count = R.summary.groups;
end
if want_groups && o.Estimate == 0
    [R.Groups, R.RegionVoxels, R.NumOther, R.Covering] = elosa_read_groups(prefix);
end
R.command = cmd;
if ~strcmp(o.Progress, 'none') && toc(started) >= 1
    if o.Estimate > 0
        fprintf('ELOSA: estimated %.3g groups (standard error %.2g) in %.1f s\n', R.estimate, R.se, toc(started));
    else
        fprintf('ELOSA: %d groups in %.1f s\n', R.count, toc(started));
    end
end
end

function [S, by_size, by_region] = read_summary(file)
txt = fileread(file);
lines = regexp(txt, '\r?\n', 'split');
S = struct();
by_size = table();
by_region = table();
section = 0;
rows = [];
for i = 1:numel(lines)
    L = strtrim(lines{i});
    if isempty(L), continue; end
    if startsWith(L, '#')
        [by_size, by_region] = store(section, rows, by_size, by_region);
        section = section + 1;
        rows = [];
        continue
    end
    parts = strsplit(L, sprintf('\t'));
    if section == 0
        S.(parts{1}) = str2double(parts{2});
    else
        v = str2double(parts);
        if all(~isnan(v)), rows(end + 1, :) = v; end %#ok<AGROW>
    end
end
[by_size, by_region] = store(section, rows, by_size, by_region);
end

function [by_size, by_region] = store(section, rows, by_size, by_region)
if section == 1 && ~isempty(rows)
    by_size = array2table(rows, 'VariableNames', {'patients', 'groups'});
elseif section == 2 && ~isempty(rows)
    by_region = array2table(rows, 'VariableNames', {'voxels_from', 'voxels_to', 'groups'});
end
end

function delete_files(prefix)
for ext = {'.in', '.summary.txt', '.progress.txt', '.groups', '.members'}
    f = [prefix ext{1}];
    try
        if exist(f, 'file'), delete(f); end
    catch
    end
end
end

function pid = matlab_pid()
% elosa stops by itself if this MATLAB session ends (e.g. is killed)
try
    pid = num2str(feature('getpid'));
catch
    pid = '0';          % unknown: don't watch
end
end
