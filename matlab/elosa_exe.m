function exe = elosa_exe()
%ELOSA_EXE Path of the compiled elosa program (see elosa_compile).
%   Errors if the program is missing or older than src/elosa.cpp, so an
%   out-of-date build is never run by mistake.
root = fileparts(fileparts(mfilename('fullpath')));
if ispc, name = 'elosa.exe'; else, name = 'elosa'; end
exe = fullfile(root, 'bin', name);
if ~exist(exe, 'file')
    error('elosa:exe', 'Cannot find %s. Run elosa_compile first (see README.md).', exe);
end
src = dir(fullfile(root, 'src', 'elosa.cpp'));
bin = dir(exe);
if ~isempty(src) && bin.datenum < src.datenum
    error('elosa:stale', '%s is older than src%selosa.cpp, which has been updated. Run elosa_compile.', ...
        exe, filesep);
end
end
