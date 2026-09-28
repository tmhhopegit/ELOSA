function [status, out, shown] = elosa_system(exe, args, poll, interval)
%ELOSA_SYSTEM Run a program with arguments, without going through a shell.
%   [status, out, shown] = elosa_system(exe, {arg1, arg2, ...})
%   [status, out, shown] = elosa_system(exe, args, poll, interval)
%   status: exit code (-1 if the program could not be started, -2 if cancelled)
%   out:    everything it printed (stdout and stderr)
%   shown:  the command line, for messages
%   poll:   optional function called every `interval` seconds (default 0.5)
%           while the program runs; if it returns true the program is stopped.
%
%   MATLAB's system() runs commands through cmd.exe on Windows, which mangles
%   quoted paths and also runs any "AutoRun" command set up for the Command
%   Prompt (a common source of "The system cannot find the path specified").
%   So the program is started directly with Java's ProcessBuilder, with each
%   argument passed separately. The MinGW-w64 bin folder (if MATLAB has one)
%   is put on the program's PATH so a compiler finds its helper programs.
%   If MATLAB is interrupted (Ctrl+C) the program is stopped too.
if nargin < 2, args = {}; end
if nargin < 3, poll = []; end
if nargin < 4 || isempty(interval), interval = 0.5; end
args = cellfun(@char, args, 'UniformOutput', false);
shown = strjoin(cellfun(@quote, [{exe}, args], 'UniformOutput', false), ' ');
bin = '';
if ispc, bin = elosa_mingw_bin(); end

if ~usejava('jvm')
    % no Java: fall back to the shell (no progress, no cancelling)
    old_path = getenv('PATH');
    restore = onCleanup(@() setenv('PATH', old_path));
    if ~isempty(bin), setenv('PATH', [bin pathsep old_path]); end
    cmd = [shown ' 2>&1'];
    if ispc, cmd = ['"' cmd '"']; end   % cmd.exe strips one outer pair of quotes
    [status, out] = system(cmd);
    return
end

list = java.util.ArrayList();
list.add(java.lang.String(exe));
for i = 1:numel(args), list.add(java.lang.String(args{i})); end
pb = java.lang.ProcessBuilder(list);
pb.directory(java.io.File(pwd));            % relative paths as in MATLAB
% output goes to a file: a pipe that nobody reads during a long run can fill
% up and stall the program
logfile = [tempname '.log'];
pb.redirectErrorStream(true);
pb.redirectOutput(java.io.File(logfile));
if ~isempty(bin)
    env = pb.environment();
    old = env.get('PATH');
    if isempty(old), old = ''; else, old = char(old); end
    env.put('PATH', [bin pathsep old]);
end
try
    proc = pb.start();
catch err
    status = -1;
    out = sprintf('Could not start %s:\n%s', exe, err.message);
    return
end
stopper = onCleanup(@() stop_and_clean(proc, logfile));   % also runs on Ctrl+C

cancelled = false;
while proc.isAlive()
    if ~isempty(poll) && poll()
        cancelled = true;
        proc.destroyForcibly();
        break
    end
    % pause (not a blocking Java call) so Ctrl+C and the Cancel button are
    % noticed straight away; Ctrl+C then runs stop_and_clean via onCleanup
    pause(interval);
end
proc.waitFor();
status = double(proc.exitValue());
if cancelled, status = -2; end
if exist(logfile, 'file'), out = fileread(logfile); else, out = ''; end
end

function stop_and_clean(proc, logfile)
if proc.isAlive()
    proc.destroyForcibly();
    proc.waitFor(5, java.util.concurrent.TimeUnit.SECONDS);   % let Windows release the log file
end
try
    if exist(logfile, 'file'), delete(logfile); end
catch
end
end

function s = quote(s)
if any(s == ' ') || isempty(s), s = ['"' s '"']; end
end
