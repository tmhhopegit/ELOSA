function elosa_compile(compiler)
%ELOSA_COMPILE Build bin/elosa(.exe) from src/elosa.cpp.
%   elosa_compile            use MATLAB's MinGW-w64 compiler on Windows
%                            (the "MATLAB Support for MinGW-w64" add-on),
%                            or g++ on Linux/macOS
%   elosa_compile('C:\path\to\g++.exe')   use a specific g++ or clang++
%   With Visual Studio, build from a "Developer Command Prompt" instead:
%       cl /O2 /std:c++17 /EHsc /Fe:bin\elosa.exe src\elosa.cpp
%
%   It tries the fastest settings first and falls back to safer ones, and
%   writes everything the compiler printed to bin/compile_log.txt.
root = fileparts(fileparts(mfilename('fullpath')));
src = fullfile(root, 'src', 'elosa.cpp');
bindir = fullfile(root, 'bin');
if ~exist(bindir, 'dir'), mkdir(bindir); end
if ispc, out = fullfile(bindir, 'elosa.exe'); else, out = fullfile(bindir, 'elosa'); end
logfile = fullfile(bindir, 'compile_log.txt');

if nargin < 1 || isempty(compiler)
    if ispc
        bin = elosa_mingw_bin();
        if isempty(bin)
            error('elosa:compile', ['No MinGW-w64 compiler found. Install the "MATLAB Support for MinGW-w64 ' ...
                'C/C++ Compiler" add-on, or pass the path to g++.exe.']);
        end
        compiler = fullfile(bin, 'g++.exe');
    else
        compiler = 'g++';
    end
end

if ispc && ~exist(compiler, 'file')
    [folder, ~] = fileparts(compiler);
    listing = dir(fullfile(folder, '*g++*'));
    error('elosa:compile', 'The compiler %s does not exist. g++ files in %s: %s', compiler, folder, ...
        strjoin({listing.name}, ', '));
end

% check the compiler starts at all before trying to build
[status, msg, shown] = elosa_system(compiler, {'--version'});
lines = {sprintf('> %s\nexit status %d\n%s\n', shown, status, msg)};
if status ~= 0
    write_log(logfile, lines);
    error('elosa:compile', 'The compiler does not start:\n%s\n%s', shown, msg);
end
fprintf('Using %s\n', strtrim(strtok(msg, newline)));

attempts = {
    {'-march=native', '-pthread', '-static'},              'optimised for this CPU'
    {'-pthread', '-static'},                               'generic CPU'
    {'-static', '-DELOSA_WIN32_THREADS'},                  'Windows threads'
    {'-DELOSA_WIN32_THREADS'},                             'Windows threads, not static'
    };
if ~ispc
    attempts = attempts(1:2, :);
    attempts{1, 1} = {'-march=native', '-pthread'};
    attempts{2, 1} = {'-pthread'};
end

for a = 1:size(attempts, 1)
    if exist(out, 'file'), delete(out); end
    args = [{'-O3', '-std=c++17'}, attempts{a, 1}, {'-o', out, src}];
    fprintf('Compiling (%s)...\n', attempts{a, 2});
    [status, msg, shown] = elosa_system(compiler, args);
    lines{end + 1} = sprintf('> %s\nexit status %d\n%s\n', shown, status, msg); %#ok<AGROW>
    if status == 0 && exist(out, 'file')
        write_log(logfile, lines);
        % check that it runs
        [~, m2] = elosa_system(out, {'--help'});
        if ~contains(m2, 'usage')
            error('elosa:compile', 'Built %s, but it does not run:\n%s', out, m2);
        end
        fprintf('Built %s\n', out);
        return
    end
end
write_log(logfile, lines);
report = strjoin(lines, newline);
if numel(report) > 4000, report = [report(1:4000) sprintf('\n...')]; end
error('elosa:compile', 'Compilation failed with every setting. What the compiler printed:\n\n%s\nFull log: %s', ...
    report, logfile);
end

function write_log(file, lines)
fid = fopen(file, 'w');
if fid < 0, return; end
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end
