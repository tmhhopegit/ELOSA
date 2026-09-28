function bin = elosa_mingw_bin()
%ELOSA_MINGW_BIN The bin folder of the MinGW-w64 compiler MATLAB knows about ('' if none).
bin = '';
root = getenv('MW_MINGW64_LOC');       % set by the MinGW-w64 add-on
if isempty(root)
    try
        cc = mex.getCompilerConfigurations('C++', 'Selected');
        if ~isempty(cc) && contains(lower(cc.Name), 'mingw'), root = cc.Location; end
    catch
    end
end
if isempty(root)
    try
        cc = mex.getCompilerConfigurations('C++', 'Installed');
        for i = 1:numel(cc)
            if contains(lower(cc(i).Name), 'mingw'), root = cc(i).Location; break; end
        end
    catch
    end
end
if ~isempty(root) && exist(fullfile(root, 'bin'), 'dir')
    bin = fullfile(root, 'bin');
end
end
