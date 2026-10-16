function tests = test_preprocess_helpers
%TEST_PREPROCESS_HELPERS Fictional tests for shared publication helpers.
tests = functiontests(localfunctions);
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'common'));
end

function setup(testCase)
sandbox = tempname;
mkdir(sandbox);
testCase.TestData.sandbox = sandbox;
end

function teardown(testCase)
if isfolder(testCase.TestData.sandbox)
    rmdir(testCase.TestData.sandbox, 's');
end
end

function testPtmdecoderUsesRepositorySubmoduleWithoutConfig(testCase)
repoRoot = testCase.TestData.sandbox;
expected = fullfile(testCase.TestData.sandbox, 'external', 'PTMdecoder');

actual = preprocess_read_local_path(repoRoot);

verifyEqual(testCase, actual, canonicalPath(expected));
end

function testReadQuotedDataRootWithSpaces(testCase)
repoRoot = testCase.TestData.sandbox;
expected = fullfile(testCase.TestData.sandbox, 'input data');
mkdir(expected);
writeText(fullfile(repoRoot, 'local_paths.yaml'), ...
    sprintf('data_root: "%s"\n', expected));
originalDirectory = pwd;
directoryCleanup = onCleanup(@() cd(originalDirectory));
unrelatedDirectory = fullfile(repoRoot, 'unrelated cwd');
mkdir(unrelatedDirectory);
cd(unrelatedDirectory);
verifyEqual(testCase, preprocess_read_local_path(repoRoot, "data_root"), ...
    canonicalPath(expected));
verifyEqual(testCase, preprocess_read_local_path(repoRoot), ...
    canonicalPath(fullfile(repoRoot, 'external', 'PTMdecoder')));
end

function testMissingAndAmbiguousDataConfigurationIsRejected(testCase)
repoRoot = testCase.TestData.sandbox;
legacyDirectory = fullfile(repoRoot, 'analysis_code');
mkdir(legacyDirectory);
writeText(fullfile(legacyDirectory, 'local_paths.yaml'), ...
    sprintf('data_root: "%s"\n', repoRoot));
verifyError(testCase, @() preprocess_read_local_path(repoRoot, "data_root"), ...
    'preprocess_read_local_path:MissingLocalPaths');
for content = {sprintf('data_root: ""\n'), sprintf('other: value\n'), ...
        sprintf('data_root: first\ndata_root: second\n')}
    writeText(fullfile(repoRoot, 'local_paths.yaml'), content{1});
    verifyError(testCase, @() preprocess_read_local_path(repoRoot, "data_root"), ...
        'preprocess_read_local_path:InvalidLocalPath');
end
verifyError(testCase, @() preprocess_read_local_path(repoRoot, "other"), ...
    'preprocess_read_local_path:InvalidKey');
end

function testProtectedCandidateOverlapIsRejected(testCase)
protected = fullfile(testCase.TestData.sandbox, 'protected');
mkdir(protected);
candidate = fullfile(protected, 'candidate');

verifyError(testCase, ...
    @() preprocess_begin_candidate(candidate, string(protected)), ...
    'preprocess_begin_candidate:ProtectedCandidateRoot');
verifyFalse(testCase, isfolder([candidate, '.partial']));
end

function testExistingCandidateAndStageAreRejected(testCase)
candidate = fullfile(testCase.TestData.sandbox, 'candidate');
mkdir(candidate);
verifyError(testCase, ...
    @() preprocess_begin_candidate(candidate, strings(0, 1)), ...
    'preprocess_begin_candidate:CandidateExists');
rmdir(candidate);

mkdir([candidate, '.partial']);
verifyError(testCase, ...
    @() preprocess_begin_candidate(candidate, strings(0, 1)), ...
    'preprocess_begin_candidate:StageExists');
end

function testExactTreePublishesWithSameParentRename(testCase)
candidate = fullfile(testCase.TestData.sandbox, 'candidate');
[candidate, stage] = preprocess_begin_candidate( ...
    candidate, strings(0, 1));
expected = ["run1/filtered_result_mascot.txt"; ...
    "run1/pepSpecFile.txt"];
writeFixtureTree(stage, expected);

preprocess_publish_candidate(stage, candidate, expected);

verifyTrue(testCase, isfolder(candidate));
verifyFalse(testCase, isfolder(stage));
verifyTrue(testCase, isfile(fullfile( ...
    candidate, 'run1', 'filtered_result_mascot.txt')));
end

function testMissingFileLeavesPartialForInspection(testCase)
candidate = fullfile(testCase.TestData.sandbox, 'candidate');
[candidate, stage] = preprocess_begin_candidate( ...
    candidate, strings(0, 1));
expected = ["run1/filtered_result_mascot.txt"; ...
    "run1/pepSpecFile.txt"];
writeFixtureTree(stage, expected(1));

verifyError(testCase, ...
    @() preprocess_publish_candidate(stage, candidate, expected), ...
    'preprocess_publish_candidate:UnexpectedTree');
verifyTrue(testCase, isfolder(stage));
verifyFalse(testCase, isfolder(candidate));
end

function testExtraFileIsRejected(testCase)
candidate = fullfile(testCase.TestData.sandbox, 'candidate');
[candidate, stage] = preprocess_begin_candidate( ...
    candidate, strings(0, 1));
expected = "run1/pepSpecFile.txt";
writeFixtureTree(stage, [expected; "run1/unexpected.txt"]);

verifyError(testCase, ...
    @() preprocess_publish_candidate(stage, candidate, expected), ...
    'preprocess_publish_candidate:UnexpectedTree');
end

function testEmptyFileIsRejected(testCase)
candidate = fullfile(testCase.TestData.sandbox, 'candidate');
[candidate, stage] = preprocess_begin_candidate( ...
    candidate, strings(0, 1));
expected = "run1/pepSpecFile.txt";
mkdir(fullfile(stage, 'run1'));
writeBytes(fullfile(stage, 'run1', 'pepSpecFile.txt'), uint8([]));

verifyError(testCase, ...
    @() preprocess_publish_candidate(stage, candidate, expected), ...
    'preprocess_publish_candidate:EmptyFile');
end

function testUnexpectedEmptyDirectoryIsRejected(testCase)
candidate = fullfile(testCase.TestData.sandbox, 'candidate');
[candidate, stage] = preprocess_begin_candidate( ...
    candidate, strings(0, 1));
expected = "run1/pepSpecFile.txt";
writeFixtureTree(stage, expected);
mkdir(fullfile(stage, 'unexpected-empty-run'));

verifyError(testCase, ...
    @() preprocess_publish_candidate(stage, candidate, expected), ...
    'preprocess_publish_candidate:UnexpectedTree');
end

function testCleanPtmdecoderReturnsCommitWithoutWarning(testCase)
checkout = fullfile(testCase.TestData.sandbox, 'ptmdecoder');
expectedCommit = makeGitCheckout(checkout);

verifyWarningFree(testCase, @() preprocess_check_ptmdecoder(checkout));
actualCommit = preprocess_check_ptmdecoder(checkout);

verifyEqual(testCase, actualCommit, expectedCommit);
end

function testDirtyPtmdecoderIsRejected(testCase)
checkout = fullfile(testCase.TestData.sandbox, 'ptmdecoder');
makeGitCheckout(checkout);
writeText(fullfile(checkout, 'tracked.txt'), 'changed');

verifyError(testCase, ...
    @() preprocess_check_ptmdecoder(checkout), ...
    'preprocess_check_ptmdecoder:DirtyCheckout');
end

function commit = makeGitCheckout(checkout)
mkdir(checkout);
runGit(checkout, 'init --quiet');
runGit(checkout, 'config user.email fixture@example.invalid');
runGit(checkout, 'config user.name "Fictional Fixture"');
writeText(fullfile(checkout, 'tracked.txt'), 'fictional');
runGit(checkout, 'add tracked.txt');
runGit(checkout, 'commit --quiet -m fixture');
commit = runGit(checkout, 'rev-parse HEAD');
end

function output = runGit(checkout, arguments)
command = sprintf('git -C "%s" %s 2>&1', checkout, arguments);
[status, output] = system(command);
assert(status == 0, 'Fixture git command failed: %s', output);
output = strtrim(output);
end

function writeFixtureTree(root, relativePaths)
for idx = 1:numel(relativePaths)
    path = fullfile(root, char(relativePaths(idx)));
    parent = fileparts(path);
    if ~isfolder(parent)
        mkdir(parent);
    end
    writeBytes(path, uint8('fictional fixture'));
end
end

function writeText(path, value)
writeBytes(path, unicode2native(value, 'UTF-8'));
end

function writeBytes(path, bytes)
fid = fopen(path, 'wb');
assert(fid >= 0, 'Unable to create fictional fixture: %s', path);
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, bytes, 'uint8');
end

function value = canonicalPath(value)
value = char(java.io.File(char(value)).getCanonicalPath());
end
