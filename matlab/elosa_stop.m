function elosa_stop()
%ELOSA_STOP End every running copy of the elosa program.
%   Use this if a copy was left running, e.g. by an older version of these
%   functions or after MATLAB was interrupted in a way that did not stop it.
%   (Current versions stop it on Ctrl+C, on Cancel, and when MATLAB exits.)
if ispc
    [status, out] = system('taskkill /F /IM elosa.exe 2>&1');
else
    [status, out] = system('pkill -x elosa 2>&1');
end
if status == 0
    fprintf('Stopped running elosa programs.\n%s', out);
else
    fprintf('No elosa program was running.\n');
end
end
