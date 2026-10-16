function tests = test_run_syn_preprocess
%TEST_RUN_SYN_PREPROCESS Fictional DAT smoke tests for both public runners.
tests = functiontests(localfunctions);
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'synthetic'));
end

function testBothRunnersPublishExactFictionalTrees(testCase)
ptmdecoderRoot = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'external', 'PTMdecoder');

sandbox = [tempname ' with spaces'];
mkdir(sandbox);
buildSyntheticRepo(sandbox, ptmdecoderRoot);
runnerDir = fullfile(sandbox, 'synthetic');
originalPath = path;
originalDirectory = pwd;
cleanupFixture = onCleanup(@() restoreFixture( ...
    sandbox, originalPath, originalDirectory));
addpath(runnerDir);
cd(runnerDir);

candidateParent = fullfile(sandbox, 'candidates');
mkdir(candidateParent);
indRoot = fullfile(candidateParent, 'syn_ind');
mixRoot = fullfile(candidateParent, 'syn_mix');
indSummary = run_syn_ind_preprocess(indRoot);
mixSummary = run_syn_mix_preprocess(mixRoot);

verifyEqual(testCase, indSummary.run_count, 16);
verifyEqual(testCase, indSummary.file_count, 32);
verifyEqual(testCase, mixSummary.run_count, 18);
verifyEqual(testCase, mixSummary.file_count, 36);
verifyTree(testCase, indRoot, 'ind', 16);
verifyTree(testCase, mixRoot, 'mix', 18);
verifyError(testCase, @() run_syn_ind_preprocess(indRoot), ...
    'preprocess_begin_candidate:CandidateExists');
stageCandidate = fullfile(candidateParent, 'existing_stage');
mkdir([stageCandidate, '.partial']);
verifyError(testCase, @() run_syn_mix_preprocess(stageCandidate), ...
    'preprocess_begin_candidate:StageExists');
for dataset = {'syn_ind', 'syn_mix'}
    runner = strrep(dataset{1}, 'syn_', 'run_syn_');
    runner = [runner, '_preprocess'];
    verifyError(testCase, @() feval(runner, ...
        fullfile(sandbox, dataset{1}, 'mascot_dat', 'candidate')), ...
        'preprocess_begin_candidate:ProtectedCandidateRoot');
    verifyError(testCase, @() feval(runner, ...
        fullfile(sandbox, dataset{1}, 'preprocess_result')), ...
        'preprocess_begin_candidate:ProtectedCandidateRoot');
end
dependency = fullfile(sandbox, 'external', 'PTMdecoder', ...
    '@CFdrFilteredResultIO', 'write.m');
movefile(dependency, [dependency, '.fixture_backup']);
callPath = path;
missingCandidate = fullfile(candidateParent, 'missing_dependency');
verifyError(testCase, @() run_syn_ind_preprocess(missingCandidate), ...
    'run_syn_ind_preprocess:MissingDependency');
verifyError(testCase, @() run_syn_mix_preprocess(missingCandidate), ...
    'run_syn_mix_preprocess:MissingDependency');
verifyFalse(testCase, isfolder([missingCandidate, '.partial']));
verifyEqual(testCase, path, callPath);
clear cleanupFixture;
end

function verifyTree(testCase, root, prefix, runCount)
files = dir(fullfile(root, '**', '*'));
files = files(~[files.isdir]);
verifyEqual(testCase, numel(files), 2 * runCount);
verifyFalse(testCase, any([files.bytes] == 0));
for idx = 1:runCount
    runName = sprintf('%s%d', prefix, idx);
    runDir = fullfile(root, runName);
    runFiles = dir(fullfile(runDir, '*'));
    runFiles = runFiles(~[runFiles.isdir]);
    verifyEqual(testCase, sort({runFiles.name}), ...
        {'filtered_result_mascot.txt', 'pepSpecFile.txt'});
    expected = sprintf(['GKGGKGLGKGGAKR\n' ...
        'SYN_%s.mgf\tSYN_SCAN:101\n' ...
        'SYN_%s.mgf\tSYN_SCAN:102\n' ...
        'SYN_%s.mgf\tSYN_SCAN:103\n'], ...
        runName, runName, runName);
    verifyEqual(testCase, ...
        fileread(fullfile(runDir, 'pepSpecFile.txt')), expected);
end
end

function buildSyntheticRepo(sandbox, ptmdecoderRoot)
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
sourceRunner = fullfile(repoRoot, 'synthetic');
sourceHelpers = fullfile(repoRoot, 'common');
targetRunner = fullfile(sandbox, 'synthetic');
targetHelpers = fullfile(sandbox, 'common');
mkdir(targetRunner);
mkdir(targetHelpers);

runnerFiles = {'run_syn_ind_preprocess.m', 'run_syn_mix_preprocess.m', ...
    'syn_pepspec_dataset_config.m', 'select_syn_pepspec_entries.m'};
for idx = 1:numel(runnerFiles)
    copyfile(fullfile(sourceRunner, runnerFiles{idx}), targetRunner);
end
helperFiles = {'preprocess_read_local_path.m', ...
    'preprocess_check_ptmdecoder.m', ...
    'preprocess_begin_candidate.m', ...
    'preprocess_publish_candidate.m'};
for idx = 1:numel(helperFiles)
    copyfile(fullfile(sourceHelpers, helperFiles{idx}), targetHelpers);
end
attachPtmdecoder(sandbox, ptmdecoderRoot);
writeText(fullfile(sandbox, 'local_paths.yaml'), ...
    sprintf('data_root: %s\n', sandbox));

writeDatasetDat(sandbox, 'syn_ind', 'ind', 16);
writeDatasetDat(sandbox, 'syn_mix', 'mix', 18);
end

function writeDatasetDat(sandbox, dataset, prefix, runCount)
for idx = 1:runCount
    runName = sprintf('%s%d', prefix, idx);
    datDir = fullfile(sandbox, dataset, 'mascot_dat', runName);
    mkdir(datDir);
    writeDat(fullfile(datDir, 'fictional.dat'), ...
        ['SYN_', runName, '.mgf']);
end
end

function writeDat(filename, mgfName)
peptide = 'GKGGKGLGKGGAKR';
lines = {'synthetic-dataset-boundary', 'MODS=', ...
    ['FILE=', mgfName], 'IT_MODS=Acetyl (K)'};
for idx = 1:3
    lines{end + 1} = sprintf( ...
        ['q%d_p1=0,1535.878357,0.000000,8,%s,8,' ...
        '0000000000000000,%d;"sp|P00000|TEST_HUMAN":0:1:7:2'], ...
        idx, peptide, 51 - idx); %#ok<AGROW>
end
lines{end + 1} = 'Content-Type: application/x-Mascot; name="proteins"';
for idx = 1:3
    lines{end + 1} = sprintf( ...
        'Content-Type: application/x-Mascot; name="query%d"', idx); %#ok<AGROW>
    lines{end + 1} = 'QUERY=0,0,0,0,0,0,0,0,0'; %#ok<AGROW>
    lines{end + 1} = sprintf('TITLE=SYN_SCAN:%d', 100 + idx); %#ok<AGROW>
    lines{end + 1} = 'charge=2+'; %#ok<AGROW>
end
writeText(filename, sprintf('%s\n', lines{:}));
end

function writeText(filename, value)
fid = fopen(filename, 'w');
assert(fid >= 0, 'Unable to create fictional fixture: %s', filename);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', value);
end

function cleanupDirectory(directory)
if isfolder(directory)
    rmdir(directory, 's');
end
end

function restoreFixture(sandbox, originalPath, originalDirectory)
cd(originalDirectory);
path(originalPath);
cleanupDirectory(sandbox);
end

function attachPtmdecoder(sandbox, sourceRoot)
checkout = fullfile(sandbox, 'external', 'PTMdecoder');
if isfolder(checkout)
    return;
end
mkdir(fileparts(checkout));
[status, output] = system(sprintf('git clone --quiet --shared "%s" "%s" 2>&1', sourceRoot, checkout));
assert(status == 0, 'Fictional dependency clone failed: %s', output);
[status, output] = system(sprintf('git -C "%s" checkout --quiet %s 2>&1', ...
    checkout, gitOutput(sourceRoot, 'rev-parse HEAD')));
assert(status == 0, 'Fictional dependency checkout failed: %s', output);
end

function output = gitOutput(root, arguments)
[status, output] = system(sprintf('git -C "%s" %s 2>&1', root, arguments));
assert(status == 0, 'Git fixture command failed: %s', output);
output = strtrim(output);
end
