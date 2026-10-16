function preprocess_publish_candidate(stageRoot, candidateRoot, expectedFiles)
%PREPROCESS_PUBLISH_CANDIDATE Validate an exact tree and rename it in place.

arguments
    stageRoot (1, 1) string
    candidateRoot (1, 1) string
    expectedFiles string
end

stageRoot = canonicalPath(stageRoot);
candidateRoot = canonicalPath(candidateRoot);
if ~isfolder(stageRoot)
    error([mfilename, ':MissingStage'], ...
        'Staging directory does not exist: "%s".', stageRoot);
end
if ~strcmpi(canonicalPath(fileparts(stageRoot)), ...
        canonicalPath(fileparts(candidateRoot)))
    error([mfilename, ':DifferentParents'], ...
        'Staging and candidate roots must have the same parent.');
end
if isfile(candidateRoot) || isfolder(candidateRoot)
    error([mfilename, ':CandidateExists'], ...
        'Candidate root appeared before publish: "%s".', candidateRoot);
end

expectedFiles = normalizeExpectedFiles(expectedFiles);
expectedDirectories = expectedDirectoryList(expectedFiles);
[actualDirectories, actualFiles] = inspectTree(stageRoot);
if ~isequal(expectedDirectories, actualDirectories) || ...
        ~isequal(expectedFiles, actualFiles)
    error([mfilename, ':UnexpectedTree'], ...
        'Staging tree does not contain exactly the expected directories and files.');
end

[status, message] = movefile(stageRoot, candidateRoot);
if ~status
    error([mfilename, ':PublishFailed'], ...
        'Failed to publish candidate root "%s": %s', candidateRoot, message);
end
end

function expectedFiles = normalizeExpectedFiles(expectedFiles)
expectedFiles = expectedFiles(:);
if isempty(expectedFiles)
    error([mfilename, ':EmptyExpectedFiles'], ...
        'expectedFiles must contain at least one relative path.');
end
for idx = 1:numel(expectedFiles)
    value = replace(strtrim(expectedFiles(idx)), "\", "/");
    parts = split(value, "/");
    isAbsolute = startsWith(value, "/") || ...
        ~isempty(regexp(char(value), '^[A-Za-z]:', 'once'));
    if strlength(value) == 0 || isAbsolute || ...
            any(parts == "" | parts == "." | parts == "..")
        error([mfilename, ':InvalidExpectedPath'], ...
            'Expected file paths must be normalized relative paths: "%s".', ...
            expectedFiles(idx));
    end
    expectedFiles(idx) = value;
end
if numel(unique(expectedFiles)) ~= numel(expectedFiles)
    error([mfilename, ':DuplicateExpectedPath'], ...
        'expectedFiles contains duplicate paths.');
end
expectedFiles = sort(expectedFiles);
end

function directories = expectedDirectoryList(expectedFiles)
directories = strings(0, 1);
for idx = 1:numel(expectedFiles)
    parent = string(fileparts(char(replace( ...
        expectedFiles(idx), "/", filesep))));
    while strlength(parent) > 0 && parent ~= "."
        directories(end + 1, 1) = replace(parent, "\", "/"); %#ok<AGROW>
        parent = string(fileparts(char(parent)));
    end
end
directories = sort(unique(directories));
end

function [directories, files] = inspectTree(root)
items = dir(fullfile(root, '**', '*'));
directoryItems = items([items.isdir]);
fileItems = items(~[items.isdir]);

directories = strings(0, 1);
for idx = 1:numel(directoryItems)
    if any(strcmp(directoryItems(idx).name, {'.', '..'}))
        continue;
    end
    path = fullfile(directoryItems(idx).folder, directoryItems(idx).name);
    directories(end + 1, 1) = normalizeRelativePath( ...
        relativePath(root, path)); %#ok<AGROW>
end
directories = sort(directories);

files = strings(numel(fileItems), 1);
for idx = 1:numel(fileItems)
    path = fullfile(fileItems(idx).folder, fileItems(idx).name);
    if fileItems(idx).bytes == 0
        error([mfilename, ':EmptyFile'], ...
            'Generated file is empty: "%s".', path);
    end
    files(idx) = normalizeRelativePath(relativePath(root, path));
end
files = sort(files);
end

function relative = relativePath(root, path)
root = canonicalPath(root);
path = canonicalPath(path);
prefix = root;
if ~endsWith(prefix, filesep)
    prefix = [prefix, filesep];
end
if ~startsWith(lower(path), lower(prefix))
    error([mfilename, ':PathOutsideRoot'], ...
        'Path "%s" is outside root "%s".', path, root);
end
relative = path(numel(prefix) + 1:end);
end

function value = normalizeRelativePath(value)
value = replace(string(value), "\", "/");
end

function value = canonicalPath(value)
value = char(java.io.File(char(value)).getCanonicalPath());
end
