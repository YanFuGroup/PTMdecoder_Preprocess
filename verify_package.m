function results = verify_package
%VERIFY_PACKAGE Run every fictional preprocessing test without skips.
repoRoot = fileparts(mfilename('fullpath'));
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
results = runtests(fullfile(repoRoot, 'tests'), 'IncludeSubfolders', true);
assertSuccess(results);
if any([results.Incomplete])
    error('verify_package:IncompleteTests', 'Tests were skipped or incomplete.');
end
fprintf('PACKAGE_TESTS: %d Passed, %d Failed, %d Incomplete\n', ...
    sum([results.Passed]), sum([results.Failed]), sum([results.Incomplete]));
end
