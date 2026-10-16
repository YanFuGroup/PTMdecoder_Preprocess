function summary = run_dna_histone_preprocess(candidateRoot)
%RUN_DNA_HISTONE_PREPROCESS Rebuild all DNA histone preprocess inputs.

arguments
    candidateRoot (1, 1) string
end

runNames = {'wt1', 'wt2', 'wt3', 'ir1', 'ir2', 'ir3'};
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
commonDir = fullfile(repoRoot, 'common');
addpath(commonDir);
dataRoot = preprocess_read_local_path(repoRoot, "data_root");
inputRoot = fullfile(dataRoot, 'DNA_histone', 'mascot_dat');
protectedRoots = string({ ...
    inputRoot
    fullfile(dataRoot, 'DNA_histone', 'results_backup')
    fullfile(dataRoot, 'DNA_histone', 'results_auto')
    });

validateInputs(inputRoot, runNames);
ptmdecoderRoot = preprocess_read_local_path(repoRoot);
validateDependencies(ptmdecoderRoot, scriptDir, ...
    commonDir);
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
    fprintf('[%d/%d] Rebuilding DNA_histone/%s\n', ...
        idx, numel(runNames), runName);
    runSummaries(idx) = processRun( ...
        fullfile(inputRoot, runName), stageRoot, runName);
end
preprocess_publish_candidate(stageRoot, candidateRoot, expectedFiles);

summary = struct( ...
    'candidate_root', candidateRoot, ...
    'ptmdecoder_commit', ptmdecoderCommit, ...
    'run_count', numel(runNames), ...
    'file_count', numel(expectedFiles), ...
    'runs', {runSummaries});
clear pathCleanup;
end

function runSummary = processRun(inputDir, stageRoot, runName)
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
selectedEntries = select_dna_histone_pepspec_entries(fdrFilteredResult);

outputDir = fullfile(stageRoot, runName);
createDirectory(outputDir);
write_mascot_result_table(groupResult, ...
    fullfile(outputDir, 'group_result_mascot.txt'));
write_mascot_result_table(fdrFilteredResult, ...
    fullfile(outputDir, 'filtered_result_mascot.txt'));
write_peptide_spectra_list_file(selectedEntries, ...
    fullfile(outputDir, 'pepSpecFile.txt'));

runSummary = struct( ...
    'run', runName, ...
    'group_count', numel(groupResult), ...
    'filtered_count', numel(fdrFilteredResult), ...
    'selected_count', numel(selectedEntries));
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

function validateDependencies(ptmdecoderRoot, scriptDir, ...
        commonDir)
fdrDir = fullfile(ptmdecoderRoot, 'FDR_control_generate_pep_spec_list');
requiredFiles = {
    fullfile(fdrDir, 'ReadDatResultFolder.m')
    fullfile(fdrDir, 'JudgeGroup.m')
    fullfile(fdrDir, 'ComputeFDR.m')
    fullfile(fdrDir, 'write_mascot_result_table.m')
    fullfile(fdrDir, 'write_peptide_spectra_list_file.m')
    fullfile(commonDir, 'FilterWithChemPrior.m')
    fullfile(commonDir, 'select_unique_dataset_scan.m')
    fullfile(scriptDir, 'select_dna_histone_pepspec_entries.m')
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
    'filtered_result_mascot.txt', 'pepSpecFile.txt'};
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
