function commit = preprocess_check_ptmdecoder(ptmdecoderRoot)
%PREPROCESS_CHECK_PTMDECODER Require a tracked-clean checkout and return its commit.

arguments
    ptmdecoderRoot (1, 1) string
end

ptmdecoderRoot = char(java.io.File(char(ptmdecoderRoot)).getCanonicalPath());
if ~isfolder(ptmdecoderRoot)
    error([mfilename, ':MissingCheckout'], ...
        'PTMdecoder checkout does not exist: "%s".', ptmdecoderRoot);
end

commit = gitOutput(ptmdecoderRoot, 'rev-parse HEAD');
status = gitOutput(ptmdecoderRoot, ...
    'status --porcelain=v1 --untracked-files=no');
if ~isempty(strtrim(status))
    error([mfilename, ':DirtyCheckout'], ...
        'Tracked PTMdecoder changes found in "%s".', ptmdecoderRoot);
end
end

function output = gitOutput(root, arguments)
if contains(root, '"')
    error([mfilename, ':UnsupportedGitPath'], ...
        'Git paths containing quotes are unsupported.');
end
command = sprintf('git -C "%s" %s 2>&1', root, arguments);
[status, output] = system(command);
output = strtrim(output);
if status ~= 0
    error([mfilename, ':GitCommandFailed'], ...
        'Git command failed in "%s": %s', root, output);
end
end
