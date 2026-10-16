function summary = run_mcf7_saha_top1_preprocess(candidateRoot, inputRoot)
%RUN_MCF7_SAHA_TOP1_PREPROCESS Rebuild SAHA top-1 preprocess inputs.

arguments
    candidateRoot (1, 1) string
    inputRoot (1, 1) string = ""
end

experimentName = 'MCF7_SAHA';
runNames = {'DMSO_R1', 'DMSO_R2', 'SAHA_R1', 'SAHA_R2'};
scriptDir = fileparts(mfilename('fullpath'));
mcf7Dir = fileparts(scriptDir);
repoRoot = fileparts(mcf7Dir);

commonDir = fullfile(repoRoot, 'common');
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(commonDir);

dataRoot = preprocess_read_local_path(repoRoot, "data_root");
experimentRoot = fullfile(dataRoot, experimentName);
if strlength(inputRoot) == 0
    inputRoot = string(fullfile(experimentRoot, 'mascot_dat'));
else
    inputRoot = string(java.io.File(char(inputRoot)).getCanonicalPath());
end
inputExperimentRoot = fileparts(inputRoot);
protectedRoots = string({ ...
    inputRoot
    fullfile(experimentRoot, 'results')
    fullfile(experimentRoot, 'results_auto')
    fullfile(experimentRoot, 'results_top1')
    fullfile(experimentRoot, 'results_top1_auto')
    fullfile(experimentRoot, 'results_top1_candidate_auto')
    fullfile(inputExperimentRoot, 'results')
    fullfile(inputExperimentRoot, 'results_auto')
    fullfile(inputExperimentRoot, 'results_top1')
    fullfile(inputExperimentRoot, 'results_top1_auto')
    fullfile(inputExperimentRoot, 'results_top1_candidate_auto')
    });

addpath(mcf7Dir);
addpath(scriptDir);
validateInputs(inputRoot, runNames);
ptmdecoderRoot = preprocess_read_local_path(repoRoot);
validateDependencies(ptmdecoderRoot, mcf7Dir, ...
    commonDir, scriptDir);
ptmdecoderCommit = preprocess_check_ptmdecoder(ptmdecoderRoot);

addpath(ptmdecoderRoot);
addpath(fullfile(ptmdecoderRoot, 'FDR_control_generate_pep_spec_list'));

[candidateRoot, stageRoot] = preprocess_begin_candidate( ...
    candidateRoot, protectedRoots);
expectedFiles = expectedOutputFiles(runNames);
runSummaries = repmat(struct( ...
    'run', '', 'group_count', 0, 'filtered_count', 0, ...
    'selected_count', 0), numel(runNames), 1);
for idx = 1:numel(runNames)
    runName = runNames{idx};
    fprintf('[%d/%d] Rebuilding %s/%s\n', ...
        idx, numel(runNames), experimentName, runName);
    runSummaries(idx) = processRun( ...
        fullfile(inputRoot, runName), stageRoot, runName, experimentName);
end
preprocess_publish_candidate(stageRoot, candidateRoot, expectedFiles);

summary = struct( ...
    'candidate_root', candidateRoot, ...
    'input_root', inputRoot, ...
    'ptmdecoder_commit', ptmdecoderCommit, ...
    'run_count', numel(runNames), ...
    'file_count', numel(expectedFiles), ...
    'runs', {runSummaries});
clear pathCleanup;
end

function runSummary = processRun( ...
        inputDir, stageRoot, runName, experimentName)
result = ReadDatResultFolder(inputDir);
if isempty(result)
    error([mfilename, ':EmptyDatResult'], ...
        'No Mascot PSMs were parsed from "%s".', inputDir);
end
result = FilterWithChemPrior(result);
if isempty(result)
    error([mfilename, ':EmptyChemPriorResult'], ...
        'Chemical-prior filtering removed every PSM from "%s".', inputDir);
end

[decoyType, groupType, ~, scores, numResult, indices] = JudgeGroup( ...
    result, 'Protein', 'REVERSE_', {'_HUMAN'});
[~, identified] = ComputeFDR(decoyType, groupType, scores, ...
    numResult, indices, 0.01);
groupResult = result(indices(~groupType));
fdrFilteredResult = result(identified.SF);

selectedMask = select_unique_dataset_scan(fdrFilteredResult);
uniqueResult = fdrFilteredResult(selectedMask);
selectorMask = select_mcf7_modified_peptides( ...
    uniqueResult, experimentName);
finalResult = uniqueResult(selectorMask);
[~, order] = sort(string({finalResult.peptide}), 'ascend');
outputResult = finalResult(order);

outputDir = fullfile(stageRoot, runName);
createDirectory(outputDir);
write_mcf7_mascot_result_table(groupResult, ...
    fullfile(outputDir, 'group_result_mascot.txt'));
write_mcf7_mascot_result_table(fdrFilteredResult, ...
    fullfile(outputDir, 'filtered_result_mascot.txt'));
write_peptide_spectra_list_file(outputResult, ...
    fullfile(outputDir, 'pepSpecFile.txt'));
write_report_msms_top1(outputResult, ...
    fullfile(outputDir, 'report_msms_top1.txt'));

runSummary = struct( ...
    'run', runName, ...
    'group_count', numel(groupResult), ...
    'filtered_count', numel(fdrFilteredResult), ...
    'selected_count', numel(outputResult));
end

function validateInputs(inputRoot, runNames)
if ~isfolder(inputRoot)
    error([mfilename, ':MissingDatRoot'], ...
        'Mascot DAT root does not exist: "%s".', inputRoot);
end
for idx = 1:numel(runNames)
    inputDir = fullfile(inputRoot, runNames{idx});
    if ~isfolder(inputDir)
        error([mfilename, ':MissingDatDirectory'], ...
            'Mascot DAT run directory does not exist: "%s".', inputDir);
    end
    if isempty(dir(fullfile(inputDir, '*.dat')))
        error([mfilename, ':MissingDatFiles'], ...
            'No .dat files found in: "%s".', inputDir);
    end
end
end

function validateDependencies(ptmdecoderRoot, mcf7Dir, ...
        commonDir, top1Dir)
fdrDir = fullfile(ptmdecoderRoot, 'FDR_control_generate_pep_spec_list');
requiredFiles = {
    fullfile(fdrDir, 'ReadDatResultFolder.m')
    fullfile(fdrDir, 'JudgeGroup.m')
    fullfile(fdrDir, 'ComputeFDR.m')
    fullfile(fdrDir, 'write_mascot_result_table.m')
    fullfile(fdrDir, 'write_peptide_spectra_list_file.m')
    fullfile(commonDir, 'FilterWithChemPrior.m')
    fullfile(top1Dir, 'write_report_msms_top1.m')
    fullfile(commonDir, 'select_unique_dataset_scan.m')
    fullfile(mcf7Dir, 'select_mcf7_modified_peptides.m')
    fullfile(mcf7Dir, 'write_mcf7_mascot_result_table.m')
    };
for idx = 1:numel(requiredFiles)
    if ~isfile(requiredFiles{idx})
        error([mfilename, ':MissingDependency'], ...
            'Missing required dependency: "%s".', requiredFiles{idx});
    end
end
conflictingFiles = {
    fullfile(ptmdecoderRoot, 'FilterWithChemPrior.m')
    fullfile(fdrDir, 'FilterWithChemPrior.m')
    fullfile(ptmdecoderRoot, 'write_report_msms_top1.m')
    fullfile(fdrDir, 'write_report_msms_top1.m')
    };
for idx = 1:numel(conflictingFiles)
    if isfile(conflictingFiles{idx})
        error([mfilename, ':ConflictingHelper'], ...
            ['PTMdecoder checkout still contains migrated helper "%s". ' ...
            'Use a checkout with the dataset-specific helpers removed.'], ...
            conflictingFiles{idx});
    end
end
end

function expectedFiles = expectedOutputFiles(runNames)
filenames = {'group_result_mascot.txt', ...
    'filtered_result_mascot.txt', 'pepSpecFile.txt', ...
    'report_msms_top1.txt'};
expectedFiles = strings(numel(runNames) * numel(filenames), 1);
idxExpected = 0;
for idxRun = 1:numel(runNames)
    for idxFile = 1:numel(filenames)
        idxExpected = idxExpected + 1;
        expectedFiles(idxExpected) = fullfile( ...
            runNames{idxRun}, filenames{idxFile});
    end
end
end

function createDirectory(directory)
[status, message] = mkdir(directory);
if ~status
    error([mfilename, ':CreateRunOutputFailed'], ...
        'Failed to create output directory "%s": %s', directory, message);
end
end
