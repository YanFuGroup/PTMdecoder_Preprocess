function [candidateRoot, stageRoot] = preprocess_begin_candidate( ...
        candidateRoot, protectedRoots)
%PREPROCESS_BEGIN_CANDIDATE Validate paths and create deterministic staging.

arguments
    candidateRoot (1, 1) string
    protectedRoots string = strings(0, 1)
end

if strlength(strtrim(candidateRoot)) == 0
    error([mfilename, ':InvalidCandidateRoot'], ...
        'candidateRoot must be nonempty.');
end
candidateRoot = canonicalPath(candidateRoot);
candidateParent = fileparts(candidateRoot);
if isempty(candidateParent) || ~isfolder(candidateParent)
    error([mfilename, ':MissingCandidateParent'], ...
        'Candidate parent directory does not exist: "%s".', candidateParent);
end

protectedRoots = protectedRoots(:);
for idx = 1:numel(protectedRoots)
    if strlength(strtrim(protectedRoots(idx))) == 0
        error([mfilename, ':InvalidProtectedRoot'], ...
            'Protected roots must be nonempty.');
    end
    protectedRoot = canonicalPath(protectedRoots(idx));
    if pathsOverlap(candidateRoot, protectedRoot)
        error([mfilename, ':ProtectedCandidateRoot'], ...
            'Candidate root overlaps protected root "%s".', protectedRoot);
    end
end

stageRoot = [candidateRoot, '.partial'];
if isfile(candidateRoot) || isfolder(candidateRoot)
    error([mfilename, ':CandidateExists'], ...
        'Refusing to overwrite candidate root: "%s".', candidateRoot);
end
if isfile(stageRoot) || isfolder(stageRoot)
    error([mfilename, ':StageExists'], ...
        ['Staging path already exists. Inspect and remove it before ' ...
        'retrying: "%s".'], stageRoot);
end

[status, message] = mkdir(stageRoot);
if ~status
    error([mfilename, ':CreateStageFailed'], ...
        'Failed to create staging directory "%s": %s', stageRoot, message);
end
end

function value = canonicalPath(value)
value = char(java.io.File(char(value)).getCanonicalPath());
end

function tf = pathsOverlap(left, right)
tf = pathIsWithin(left, right) || pathIsWithin(right, left);
end

function tf = pathIsWithin(path, root)
path = lower(canonicalPath(path));
root = lower(canonicalPath(root));
if strcmp(path, root)
    tf = true;
    return;
end
if ~endsWith(root, filesep)
    root = [root, filesep];
end
tf = startsWith(path, root);
end
