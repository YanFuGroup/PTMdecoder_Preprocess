function tests = test_mcf7_preprocess
%TEST_MCF7_PREPROCESS Fictional tests for MCF7 auto and Top 1 runners.
tests = functiontests(localfunctions);
addpath(repoPath('mcf7'));
end

function testWriterGoldenFormatPreserved(testCase)
preprocessDir = repoPath('mcf7', 'top1');
addpath(preprocessDir);
cleanupPath = onCleanup(@() rmpath(preprocessDir));

result = makePsmStruct( ...
    {'SYN_A.mgf', 'SYN_A.mgf'}, ...
    {'SYN_SCAN:101', 'SYN_SCAN:102'}, ...
    {'SYNPEPTIDEA', 'SYNPEPTIDEA'}, ...
    {'Acetyl (K)', 'Acetyl (K)'}, ...
    {'3', '3'});
outputDir = temporaryDirectory();
cleanupOutput = onCleanup(@() cleanupDirectory(outputDir));
outputFile = fullfile(outputDir, 'report_msms_top1.txt');
write_report_msms_top1(result, outputFile);

expectedLines = { ...
    'P\tSYNPEPTIDEA', ...
    'S\tSYN_A.mgf\tSYN_SCAN:101', ...
    '_SYN{Acetyl}PEPTIDEA_\t1.000000', ...
    'S\tSYN_A.mgf\tSYN_SCAN:102', ...
    '_SYN{Acetyl}PEPTIDEA_\t1.000000' ...
    };
expectedLines = cellfun(@(value) strrep(value, '\t', char(9)), ...
    expectedLines, 'UniformOutput', false);
expected = [strjoin(expectedLines, newline), newline];
verifyEqual(testCase, normalizeNewlines(fileread(outputFile)), char(expected));
clear cleanupOutput;
end

function testExplicitSelectorRules(testCase)
result = makePsmStruct( ...
    {'A.mgf', 'A.mgf', 'A.mgf', 'A.mgf'}, ...
    {'1', '2', '3', '4'}, ...
    {'PEPTIDEA', 'PEPTIDEB', 'PEPTIDEC', 'PEPTIDED'}, ...
    {'Acetyl (K)', '-', 'Acetyl (K),Oxidation (M)', 'Methyl (K)'}, ...
    {'3', '-', '2', '3'});

jib04Mask = select_mcf7_modified_peptides(result, 'MCF7_JIB04');
sahaMask = select_mcf7_modified_peptides(result, 'MCF7_SAHA');

verifyEqual(testCase, jib04Mask, [true, false, false, true]);
verifyEqual(testCase, sahaMask, [true, false, false, true]);
verifyError(testCase, ...
    @() select_mcf7_modified_peptides(result, 'MCF7_OTHER'), ...
    'select_mcf7_modified_peptides:UnknownExperiment');
end

function testTableWriterPreservesSignedZero(testCase)
ptmdecoderRoot = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'external', 'PTMdecoder');
originalPath = path;
cleanupPath = onCleanup(@() path(originalPath));
addpath(fullfile(ptmdecoderRoot, 'FDR_control_generate_pep_spec_list'));
addpath(ptmdecoderRoot);
outputDir = temporaryDirectory();
cleanupOutput = onCleanup(@() cleanupDirectory(outputDir));
entry = struct('Site', 'fictional.dat', 'DatasetName', 'fictional.mgf', ...
    'Scan', '1', 'Spectrum', 'fictional_scan', 'Charge', '2', ...
    'Calc_neutral_pepmass', 1000, 'precursor_neutral_mass', 1000, ...
    'massdiff', str2double('-0.000000'), 'num_match_ions', '8', ...
    'peptide', 'SYNPEPTIDEA', 'protein', 'TEST_HUMAN', ...
    'modification', 'Acetyl (K)', 'modificationlocation', '3', 'Score', 50);
entries = repmat(entry, 1, 3);
entries(2).massdiff = 0;
entries(3).massdiff = -0.000123;
original = entries;
outputFile = fullfile(outputDir, 'filtered_result_mascot.txt');
write_mcf7_mascot_result_table(entries, outputFile);
lines = readlines(outputFile);
negative = split(lines(2), char(9));
positive = split(lines(3), char(9));
nonzero = split(lines(4), char(9));
verifyEqual(testCase, negative(8), "-0.000000");
verifyEqual(testCase, positive(8), "0.000000");
verifyEqual(testCase, nonzero(8), "-0.000123");
verifyEqual(testCase, entries, original);
end

function testFourRunnersExecuteFictionalDatFullChain(testCase)
ptmdecoderRoot = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'external', 'PTMdecoder');

sandbox = temporaryDirectory();
originalPath = path;
originalDirectory = pwd;
cleanupFixture = onCleanup(@() restoreFixture( ...
    sandbox, originalPath, originalDirectory));
buildSyntheticRepo(sandbox, ptmdecoderRoot);
addpath(fullfile(sandbox, 'mcf7'));
addpath(fullfile(sandbox, 'mcf7', 'top1'));

candidateParent = fullfile(sandbox, 'candidates');
mkdir(candidateParent);
jib04Candidate = fullfile(candidateParent, 'jib04');
sahaCandidate = fullfile(candidateParent, 'saha');
jib04Auto = fullfile(candidateParent, 'jib04_auto');
sahaAuto = fullfile(candidateParent, 'saha_auto');
jib04Summary = run_mcf7_jib04_top1_preprocess(jib04Candidate);
sahaSummary = run_mcf7_saha_top1_preprocess(sahaCandidate);
jib04AutoSummary = run_mcf7_jib04_auto_preprocess(jib04Auto);
sahaAutoSummary = run_mcf7_saha_auto_preprocess(sahaAuto);

verifyEqual(testCase, jib04Summary.run_count, 4);
verifyEqual(testCase, jib04Summary.file_count, 16);
verifyEqual(testCase, sahaSummary.run_count, 4);
verifyEqual(testCase, sahaSummary.file_count, 16);
verifyEqual(testCase, jib04AutoSummary.file_count, 12);
verifyEqual(testCase, sahaAutoSummary.file_count, 12);
verifyEqual(testCase, jib04Summary.ptmdecoder_commit, ...
    gitOutput(ptmdecoderRoot, 'rev-parse HEAD'));
verifyExactOutputTree(testCase, jib04Candidate, ...
    {'DMSO_1', 'DMSO_2', 'JIB04_1', 'JIB04_2'}, true);
verifyExactOutputTree(testCase, sahaCandidate, ...
    {'DMSO_R1', 'DMSO_R2', 'SAHA_R1', 'SAHA_R2'}, true);
verifyExactOutputTree(testCase, jib04Auto, ...
    {'DMSO_1', 'DMSO_2', 'JIB04_1', 'JIB04_2'}, false);
verifyExactOutputTree(testCase, sahaAuto, ...
    {'DMSO_R1', 'DMSO_R2', 'SAHA_R1', 'SAHA_R2'}, false);
verifyCommonFilesEqual(testCase, jib04Auto, jib04Candidate, ...
    {'DMSO_1', 'DMSO_2', 'JIB04_1', 'JIB04_2'});
verifyCommonFilesEqual(testCase, sahaAuto, sahaCandidate, ...
    {'DMSO_R1', 'DMSO_R2', 'SAHA_R1', 'SAHA_R2'});

% An explicit DAT root must be read as-is; Site records that actual path.
explicitInputs = fullfile(sandbox, 'external', 'MCF7_JIB04', 'mascot_dat');
mkdir(fileparts(explicitInputs));
copyfile(fullfile(sandbox, 'MCF7_JIB04', 'mascot_dat'), explicitInputs);
explicitCandidate = fullfile(candidateParent, 'explicit_inputs');
explicitSummary = run_mcf7_jib04_auto_preprocess(explicitCandidate, explicitInputs);
verifyEqual(testCase, explicitSummary.input_root, string(explicitInputs));
explicitTable = fileread(fullfile(explicitCandidate, 'DMSO_1', 'filtered_result_mascot.txt'));
verifyTrue(testCase, contains(explicitTable, ...
    fullfile(explicitInputs, 'DMSO_1', 'fictional.dat')));
verifyError(testCase, ...
    @() run_mcf7_jib04_auto_preprocess(fullfile(explicitInputs, 'outputs'), explicitInputs), ...
    'preprocess_begin_candidate:ProtectedCandidateRoot');
verifyError(testCase, @() run_mcf7_jib04_auto_preprocess( ...
    fullfile(sandbox, 'MCF7_JIB04', 'results_auto'), explicitInputs), ...
    'preprocess_begin_candidate:ProtectedCandidateRoot');
verifyError(testCase, @() run_mcf7_jib04_auto_preprocess( ...
    fullfile(fileparts(explicitInputs), 'results_auto'), explicitInputs), ...
    'preprocess_begin_candidate:ProtectedCandidateRoot');

expectedLines = { ...
    'P\tSYNPEPTIDEA', ...
    'S\tSYN_DMSO_1.mgf\tSYN_SCAN:101', ...
    '_SYN{Acetyl}PEPTIDEA_\t1.000000', ...
    'S\tSYN_DMSO_1.mgf\tSYN_SCAN:102', ...
    '_SYN{Acetyl}PEPTIDEA_\t1.000000' ...
    };
expectedLines = cellfun(@(value) strrep(value, '\t', char(9)), ...
    expectedLines, 'UniformOutput', false);
expected = [strjoin(expectedLines, newline), newline];
verifyEqual(testCase, normalizeNewlines(fileread(fullfile( ...
    jib04Candidate, 'DMSO_1', 'report_msms_top1.txt'))), char(expected));

verifyFalse(testCase, isfolder(fullfile( ...
    sandbox, 'MCF7_JIB04', 'results_top1_auto')));
verifyFalse(testCase, isfolder(fullfile( ...
    sandbox, 'MCF7_SAHA', 'results_top1_auto')));
verifyError(testCase, ...
    @() run_mcf7_jib04_top1_preprocess(jib04Candidate), ...
    'preprocess_begin_candidate:CandidateExists');
verifyError(testCase, ...
    @() run_mcf7_jib04_auto_preprocess(fullfile( ...
        sandbox, 'MCF7_JIB04', 'results_auto')), ...
    'preprocess_begin_candidate:ProtectedCandidateRoot');
verifyError(testCase, ...
    @() run_mcf7_jib04_top1_preprocess(fullfile( ...
        sandbox, 'MCF7_JIB04', 'results_top1_auto')), ...
    'preprocess_begin_candidate:ProtectedCandidateRoot');
verifyError(testCase, ...
    @() run_mcf7_saha_auto_preprocess(fullfile( ...
        sandbox, 'MCF7_SAHA', 'results_auto')), ...
    'preprocess_begin_candidate:ProtectedCandidateRoot');
verifyError(testCase, ...
    @() run_mcf7_saha_top1_preprocess(fullfile( ...
        sandbox, 'MCF7_SAHA', 'results_top1_auto')), ...
    'preprocess_begin_candidate:ProtectedCandidateRoot');
clear cleanupFixture;
end

function testRunnersRejectConflictingHelpersBeforeReadingDat(testCase)
sourceRoot = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'external', 'PTMdecoder');
runnerNames = {'run_mcf7_jib04_top1_preprocess', ...
    'run_mcf7_saha_top1_preprocess', 'run_mcf7_jib04_auto_preprocess', ...
    'run_mcf7_saha_auto_preprocess'};
for helperName = {'FilterWithChemPrior', 'write_report_msms_top1'}
    for location = {'root', 'fdr'}
        sandbox = temporaryDirectory();
        originalPath = path;
        originalDirectory = pwd;
        cleanupFixture = onCleanup(@() restoreFixture( ...
            sandbox, originalPath, originalDirectory));
        checkout = buildShadowedPtmdecoder( ...
            sandbox, sourceRoot, helperName{1}, location{1});
        buildSyntheticRepo(sandbox, checkout);
        addpath(fullfile(sandbox, 'mcf7'));
        addpath(fullfile(sandbox, 'mcf7', 'top1'));
        candidateParent = fullfile(sandbox, 'candidates');
        mkdir(candidateParent);
        % Reject checkout contents even when local helpers have path precedence.
        cd(fullfile(sandbox, 'common'));
        callPath = path;
        runnerCount = numel(runnerNames);
        if strcmp(helperName{1}, 'write_report_msms_top1')
            runnerCount = 2;
        end
        for idx = 1:runnerCount
            candidate = fullfile(candidateParent, runnerNames{idx});
            verifyError(testCase, @() feval(runnerNames{idx}, candidate), ...
                [runnerNames{idx}, ':ConflictingHelper']);
            verifyFalse(testCase, isfolder(candidate));
            verifyFalse(testCase, isfolder([candidate, '.partial']));
            verifyEqual(testCase, path, callPath);
        end
        clear cleanupFixture;
    end
end
end

function verifyExactOutputTree(testCase, candidateRoot, runNames, includeTop1)
files = dir(fullfile(candidateRoot, '**', '*'));
files = files(~[files.isdir]);
verifyEqual(testCase, numel(files), (3 + includeTop1) * numel(runNames));
verifyFalse(testCase, any([files.bytes] == 0));
requiredNames = sort({ ...
    'filtered_result_mascot.txt', ...
    'group_result_mascot.txt', ...
    'pepSpecFile.txt'});
if includeTop1
    requiredNames = sort([requiredNames, {'report_msms_top1.txt'}]);
end
for idx = 1:numel(runNames)
    runFiles = dir(fullfile(candidateRoot, runNames{idx}, '*'));
    runFiles = runFiles(~[runFiles.isdir]);
    verifyEqual(testCase, sort({runFiles.name}), requiredNames);
end
end

function verifyCommonFilesEqual(testCase, autoRoot, top1Root, runNames)
commonNames = {'group_result_mascot.txt', ...
    'filtered_result_mascot.txt', 'pepSpecFile.txt'};
for idxRun = 1:numel(runNames)
    for idxFile = 1:numel(commonNames)
        autoBytes = readBytes(fullfile(autoRoot, runNames{idxRun}, commonNames{idxFile}));
        top1Bytes = readBytes(fullfile(top1Root, runNames{idxRun}, commonNames{idxFile}));
        verifyEqual(testCase, autoBytes, top1Bytes);
    end
end
end

function bytes = readBytes(filename)
fid = fopen(filename, 'rb');
assert(fid >= 0, 'Unable to read fictional output: %s', filename);
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
end

function buildSyntheticRepo(sandbox, ptmdecoderRoot)
repoRoot = repoPath();
sourceMcf7 = fullfile(repoRoot, 'mcf7');
sourcePreprocess = fullfile(repoRoot, 'common');
targetMcf7 = fullfile(sandbox, 'mcf7');
targetTop1 = fullfile(targetMcf7, 'top1');
targetPreprocess = fullfile(sandbox, 'common');
mkdir(targetMcf7);
mkdir(targetTop1);
mkdir(targetPreprocess);

copyfile(fullfile(sourceMcf7, '*.m'), targetMcf7);
copyfile(fullfile(sourceMcf7, 'top1', '*.m'), targetTop1);
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

writeExperimentDat(sandbox, 'MCF7_JIB04', ...
    {'DMSO_1', 'DMSO_2', 'JIB04_1', 'JIB04_2'});
writeExperimentDat(sandbox, 'MCF7_SAHA', ...
    {'DMSO_R1', 'DMSO_R2', 'SAHA_R1', 'SAHA_R2'});
end

function checkout = buildShadowedPtmdecoder(sandbox, sourceRoot, helperName, location)
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
writeText(fullfile(directory, [helperName, '.m']), sprintf( ...
    ['function varargout = %s(varargin)\n' ...
    'error(''fixture:ExternalHelperInvoked'', ''External helper was selected.'');\nend\n'], ...
    helperName));
end

function writeExperimentDat(sandbox, experiment, runNames)
for idx = 1:numel(runNames)
    runName = runNames{idx};
    datDir = fullfile(sandbox, experiment, 'mascot_dat', runName);
    mkdir(datDir);
    writeDat(fullfile(datDir, 'fictional.dat'), ...
        ['SYN_', runName, '.mgf'], ...
        {'SYNPEPTIDEA', 'SYNPEPTIDEA', 'SYNPEPTIDEB'}, ...
        {'101', '102', '103'}, ...
        {'00010', '00010', '00000'}, ...
        {'50', '49', '48'});
end
end

function writeDat(filename, mgfName, peptides, scans, modificationCodes, scores)
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

function result = makePsmStruct( ...
        datasetNames, scans, peptides, modifications, locations)
template = struct( ...
    'DatasetName', '', 'Scan', '', 'Spectrum', '', ...
    'peptide', '', 'modification', '', 'modificationlocation', '');
result = repmat(template, size(datasetNames));
for idx = 1:numel(datasetNames)
    result(idx).DatasetName = datasetNames{idx};
    result(idx).Scan = scans{idx};
    result(idx).Spectrum = scans{idx};
    result(idx).peptide = peptides{idx};
    result(idx).modification = modifications{idx};
    result(idx).modificationlocation = locations{idx};
end
end

function rootDir = repoPath(varargin)
rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
rootDir = fullfile(rootDir, varargin{:});
end

function directory = temporaryDirectory()
scratchRoot = tempdir;
if ~isfolder(scratchRoot)
    mkdir(scratchRoot);
end
directory = tempname(scratchRoot);
mkdir(directory);
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

function value = normalizeNewlines(value)
value = regexprep(char(value), '\r\n', '\n');
end

function output = gitOutput(root, arguments)
[status, output] = system(sprintf('git -C "%s" %s 2>&1', root, arguments));
assert(status == 0, 'Git fixture command failed: %s', output);
output = strtrim(output);
end

function writeText(path, value)
fid = fopen(path, 'w');
assert(fid >= 0, 'Unable to create fictional fixture: %s', path);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', value);
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
