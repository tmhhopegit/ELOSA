function [poll, finish] = elosa_progress_display(file, style)
%ELOSA_PROGRESS_DISPLAY Show the progress the elosa program reports.
%   [poll, finish] = elosa_progress_display(progress_file, style)
%   style: 'text' (a line in the Command Window that updates in place),
%          'bar'  (a progress bar with a Cancel button), or 'none'.
%   poll() reads the file and updates the display; it returns true if the
%   user pressed Cancel. finish() clears the display.
%
%   How to read it: the percentage and time left are ROUGH. The search is a
%   tree whose branches differ enormously in size, and the program can only
%   guess how much is left from the branches it has finished. On test data
%   its guess of the total time was between 0.9x and 2x the real time, so
%   the time left is shown as a range: from half the program's guess to its
%   full guess. "Groups" and "elapsed" are exact.
if nargin < 2, style = 'text'; end
started = tic;
shown = 0;          % characters currently printed (text style)
hbar = [];
cancelled = false;
poll = @do_poll;
finish = @do_finish;

    function stop = do_poll()
        stop = false;
        if strcmp(style, 'none'), return; end
        if toc(started) < 1, return; end            % stay quiet for quick runs
        msg = message(read_progress(file));
        if isempty(msg), return; end
        if strcmp(style, 'bar')
            [txt, frac] = deal(msg{1}, msg{2});
            if isempty(hbar) || ~isvalid(hbar)
                if cancelled, stop = true; return; end
                hbar = waitbar(0, txt, 'Name', 'ELOSA', ...
                    'CreateCancelBtn', @(src, ~) set_cancel());
            end
            waitbar(max(0, min(1, frac)), hbar, txt);
            drawnow;
            stop = cancelled;
        else
            txt = msg{1};
            fprintf('%s%s', repmat(sprintf('\b'), 1, shown), txt);
            shown = numel(txt);
        end
    end

    function set_cancel()
        cancelled = true;
    end

    function do_finish()
        if shown > 0
            fprintf('%s', repmat(sprintf('\b'), 1, shown));
            shown = 0;
        end
        if ~isempty(hbar) && isvalid(hbar), delete(hbar); end
    end
end

function p = read_progress(file)
p = [];
try
    txt = strtrim(fileread(file));
catch
    return
end
parts = strsplit(txt, sprintf('\t'));
if numel(parts) < 5, return; end
v = str2double(parts(2:5));
if any(isnan(v)), return; end
p = struct('phase', parts{1}, 'fraction', v(1), 'groups', v(2), 'elapsed', v(3), 'left', v(4));
end

function msg = message(p)
msg = {};
if isempty(p), return; end
switch p.phase
    case 'reading'
        txt = sprintf('ELOSA: reading the lesions (%s)', dur(p.elapsed));
        frac = 0;
    case 'preparing'
        txt = sprintf('ELOSA: finding distinct voxel patterns (%s)', dur(p.elapsed));
        frac = 0;
    case 'estimating'
        txt = sprintf('ELOSA: estimating, %d%% of the random walks done (%s)', round(100 * p.fraction), dur(p.elapsed));
        if p.left >= 0, txt = [txt sprintf(', about %s left', dur(p.left))]; end
        frac = p.fraction;
    case 'searching'
        f = p.fraction;
        txt = sprintf('ELOSA: %s groups so far, %s elapsed, roughly %d%% done', ...
            int_text(p.groups), dur(p.elapsed), floor(100 * f));
        if p.left >= 0 && f > 0.01 && f < 1
            spent = p.left * f / (1 - f);          % time spent searching
            guess = spent + p.left;                % program's guess of the total search time
            lo = max(0, guess / 2 - spent);
            hi = guess - spent;
            if hi < 60 || lo > hi / 2
                txt = [txt sprintf(', about %s left', dur(hi))];
            else
                txt = [txt sprintf(', %s to %s left', dur(lo), dur(hi))];
            end
        else
            txt = [txt ', time left not known yet'];
        end
        frac = f;
    case {'finishing', 'done'}
        txt = sprintf('ELOSA: writing results (%s)', dur(p.elapsed));
        frac = 1;
    otherwise
        return
end
msg = {txt, frac};
end

function s = dur(sec)
if sec < 60
    s = sprintf('%ds', round(sec));
elseif sec < 3600
    s = sprintf('%dm %02ds', floor(sec / 60), floor(mod(sec, 60)));
elseif sec < 86400
    s = sprintf('%dh %02dm', floor(sec / 3600), floor(mod(sec, 3600) / 60));
else
    s = sprintf('%dd %dh', floor(sec / 86400), floor(mod(sec, 86400) / 3600));
end
end

function s = int_text(x)
s = sprintf('%d', round(x));
n = numel(s);
for k = n - 3:-3:1
    s = [s(1:k) ',' s(k + 1:end)];
end
end
