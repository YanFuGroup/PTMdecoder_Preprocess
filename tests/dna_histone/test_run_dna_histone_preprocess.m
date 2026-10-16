function tests = test_run_dna_histone_preprocess
%TEST_RUN_DNA_HISTONE_PREPROCESS Fictional DAT smoke test for the public runner.
tests = functiontests(localfunctions);
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'dna_histone'));
end

function testRunnerPublishesExactFictionalTree(testCase)
ptmdecoderRoot = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'external', 'PTMdecoder');

sandbox = tempname;
mkdir(sandbox);
runnerDir = fullfile(sandbox, 'dna_histone');
originalPath = path;
originalDirectory = pwd;
cleanupFixture = onCleanup(@() restoreFixture( ...
    sandbox, originalPath, originalDirectory));
buildSyntheticRepo(sandbox, ptmdecoderRoot);
addpath(runnerDir);

candidateRoot = fullfile(sandbox, 'candidates', 'dna_histone');
mkdir(fileparts(candidateRoot));
summary = run_dna_histone_preprocess(candidateRoot);

verifyEqual(testCase, summary.run_count, 6);
verifyEqual(testCase, summary.file_count, 18);
verifyEqual(testCase, summary.ptmdecoder_commit, ...
    gitOutput(ptmdecoderRoot, 'rev-parse HEAD'));
verifyExactOutputTree(testCase, candidateRoot, ...
    {'wt1', 'wt2', 'wt3', 'ir1', 'ir2', 'ir3'});
verifyFalse(testCase, isfolder(fullfile( ...
    sandbox, 'DNA_histone', 'results_auto')));
verifyError(testCase, ...
    @() run_dna_histone_preprocess(candidateRoot), ...
    'preprocess_begin_candidate:CandidateExists');
verifyError(testCase, @() run_dna_histone_preprocess( ...
    fullfile(sandbox, 'DNA_histone', 'mascot_dat', 'candidate')), ...
    'preprocess_begin_candidate:ProtectedCandidateRoot');
verifyError(testCase, @() run_dna_histone_preprocess( ...
    fullfile(sandbox, 'DNA_histone', 'results_auto')), ...
    'preprocess_begin_candidate:ProtectedCandidateRoot');
clear cleanupFixture;
end

function testRunnerRejectsConflictingHelperBeforeReadingDat(testCase)
sourceRoot = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'external', 'PTMdecoder');
for location = {'root', 'fdr'}
    sandbox = tempname;
    mkdir(sandbox);
    originalPath = path;
    originalDirectory = pwd;
    cleanupFixture = onCleanup(@() restoreFixture( ...
        sandbox, originalPath, originalDirectory));
    checkout = buildShadowedPtmdecoder(sandbox, sourceRoot, location{1});
    buildSyntheticRepo(sandbox, checkout);
    addpath(fullfile(sandbox, 'dna_histone'));
    candidate = fullfile(sandbox, 'candidates', 'dna_histone');
    mkdir(fileparts(candidate));
    % Reject checkout contents even when local helpers have path precedence.
    cd(fullfile(sandbox, 'common'));
    callPath = path;
    verifyError(testCase, @() run_dna_histone_preprocess(candidate), ...
        'run_dna_histone_preprocess:ConflictingHelper');
    verifyFalse(testCase, isfolder(candidate));
    verifyFalse(testCase, isfolder([candidate, '.partial']));
    verifyEqual(testCase, path, callPath);
    clear cleanupFixture;
end
end

function verifyExactOutputTree(testCase, candidateRoot, runNames)
files = dir(fullfile(candidateRoot, '**', '*'));
files = files(~[files.isdir]);
verifyEqual(testCase, numel(files), 3 * numel(runNames));
verifyFalse(testCase, any([files.bytes] == 0));
requiredNames = sort({ ...
    'filtered_result_mascot.txt', ...
    'group_result_mascot.txt', ...
    'pepSpecFile.txt'});
for idx = 1:numel(runNames)
    runFiles = dir(fullfile(candidateRoot, runNames{idx}, '*'));
    runFiles = runFiles(~[runFiles.isdir]);
    verifyEqual(testCase, sort({runFiles.name}), requiredNames);
end
end

function buildSyntheticRepo(sandbox, ptmdecoderRoot)
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
sourceDna = fullfile(repoRoot, 'dna_histone');
sourcePreprocess = fullfile(repoRoot, 'common');
targetDna = fullfile(sandbox, 'dna_histone');
targetPreprocess = fullfile(sandbox, 'common');
mkdir(targetDna);
mkdir(targetPreprocess);

copyfile(fullfile(sourceDna, 'run_dna_histone_preprocess.m'), targetDna);
copyfile(fullfile(sourceDna, 'select_dna_histone_pepspec_entries.m'), targetDna);
helperNames = { ...
    'preprocess_read_local_path.m', ...
    'preprocess_check_ptmdecoder.m', ...
    'preprocess_begin_candidate.m', ...
    'preprocess_publish_candidate.m', ...
    'select_unique_dataset_scan.m'};
for idx = 1:numel(helperNames)
    copyfile(fullfile(sourcePreprocess, helperNames{idx}), targetPreprocess);
end
copyfile(fullfile(sourcePreprocess, 'FilterWithChemPrior.m'), ...
    targetPreprocess);
attachPtmdecoder(sandbox, ptmdecoderRoot);
writeText(fullfile(sandbox, 'local_paths.yaml'), ...
    sprintf('data_root: %s\n', sandbox));

runNames = {'wt1', 'wt2', 'wt3', 'ir1', 'ir2', 'ir3'};
for idx = 1:numel(runNames)
    runName = runNames{idx};
    datDir = fullfile(sandbox, 'DNA_histone', 'mascot_dat', runName);
    mkdir(datDir);
    writeDat(fullfile(datDir, 'fictional.dat'), ...
        ['SYN_', runName, '.mgf']);
end
end

function checkout = buildShadowedPtmdecoder(sandbox, sourceRoot, location)
% Untracked migrated helpers must also be rejected in a tracked-clean checkout.
checkout = fullfile(sandbox, 'external', 'PTMdecoder');
fdrDir = fullfile(checkout, 'FDR_control_generate_pep_spec_list');
mkdir(checkout);
copyfile(fullfile(sourceRoot, '*.m'), checkout);
copyfile(fullfile(sourceRoot, 'FDR_control_generate_pep_spec_list'), fdrDir);
directory = checkout;
if strcmp(location, 'fdr')
    directory = fdrDir;
end
writeText(fullfile(fdrDir, 'ReadDatResultFolder.m'), sprintf( ...
    ['function result = ReadDatResultFolder(varargin)\n' ...
    'error(''fixture:UnexpectedDatRead'', ''DAT reader must not run.'');\nend\n']));
gitOutput(checkout, 'init');
gitOutput(checkout, 'add .');
gitOutput(checkout, ['-c user.name="Fictional Fixture" ' ...
    '-c user.email=fictional@example.invalid -c commit.gpgsign=false ' ...
    'commit -m "Add shadowed preprocess fixture"']);
writeText(fullfile(directory, 'FilterWithChemPrior.m'), sprintf( ...
    ['function result = FilterWithChemPrior(result)\n' ...
    'error(''fixture:ExternalHelperInvoked'', ''External helper was selected.'');\nend\n']));
end

function writeDat(filename, mgfName)
peptides = {'SYNPEPTIDEA', 'SYNPEPTIDEA', 'SYNPEPTIDEB'};
scans = {'101', '102', '103'};
modificationCodes = {'00010', '00010', '00000'};
scores = {'50', '49', '48'};
lines = {'synthetic-dataset-boundary', 'MODS=', ...
    ['FILE=', mgfName], 'IT_MODS=Acetyl (K),Methyl (K)'};
for idx = 1:numel(peptides)
    lines{end + 1} = sprintf( ...
        ['q%d_p1=0,1000.000000,0.000000,5,%s,8,%s,%s;' ...
        '"sp|P00000|TEST_HUMAN":0:1:7:2'], ...
        idx, peptides{idx}, modificationCodes{idx}, scores{idx}); %#ok<AGROW>
end
lines{end + 1} = 'Content-Type: application/x-Mascot; name="proteins"';
for idx = 1:numel(peptides)
    lines{end + 1} = sprintf( ...
        'Content-Type: application/x-Mascot; name="query%d"', idx); %#ok<AGROW>
    lines{end + 1} = 'QUERY=0,0,0,0,0,0,0,0,0'; %#ok<AGROW>
    lines{end + 1} = sprintf('TITLE=SYN_SCAN:%s', scans{idx}); %#ok<AGROW>
    lines{end + 1} = 'charge=2+'; %#ok<AGROW>
end
writeText(filename, sprintf('%s\n', lines{:}));
end

function output = gitOutput(root, arguments)
[status, output] = system(sprintf('git -C "%s" %s 2>&1', root, arguments));
assert(status == 0, 'Git fixture command failed: %s', output);
output = strtrim(output);
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
