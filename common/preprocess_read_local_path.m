function localPath = preprocess_read_local_path(repoRoot, key)
%PREPROCESS_READ_LOCAL_PATH Locate bundled PTMdecoder or read the input data root.

arguments
    repoRoot (1, 1) string
    key (1, 1) string = "ptmdecoder_root"
end

if ~ismember(key, ["ptmdecoder_root", "data_root"])
    error([mfilename, ':InvalidKey'], ...
        'Unsupported local path key: "%s".', key);
end
repoRoot = canonicalPath(repoRoot);
if key == "ptmdecoder_root"
    localPath = canonicalPath(fullfile( ...
        repoRoot, 'external', 'PTMdecoder'));
    return;
end
configFile = fullfile(repoRoot, 'local_paths.yaml');
if ~isfile(configFile)
    error([mfilename, ':MissingLocalPaths'], ...
        'Missing local path config: "%s".', configFile);
end

lines = readlines(configFile);
pattern = "^\s*" + key + "\s*:\s*(.*?)\s*$";
matches = regexp(lines, pattern, 'tokens', 'once');
matches = matches(~cellfun(@isempty, matches));
if numel(matches) ~= 1
    error([mfilename, ':InvalidLocalPath'], ...
        'Expected exactly one %s entry in "%s".', key, configFile);
end

value = strtrim(matches{1}{1});
if numel(value) >= 2 && ((value(1) == '"' && value(end) == '"') || ...
        (value(1) == '''' && value(end) == ''''))
    value = value(2:end - 1);
end
if isempty(value)
    error([mfilename, ':InvalidLocalPath'], ...
        '%s must be nonempty in "%s".', key, configFile);
end
localPath = canonicalPath(value);
end

function value = canonicalPath(value)
value = char(java.io.File(char(value)).getCanonicalPath());
end
