function summary = run_syn_ind_preprocess(candidateRoot)
%RUN_SYN_IND_PREPROCESS Rebuild all synthetic individual-run inputs.

arguments
    candidateRoot (1, 1) string
end

dataset = 'syn_ind';
runNames = arrayfun(@(idx) sprintf('ind%d', idx), 1:16, ...
    'UniformOutput', false);
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(scriptDir);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
commonDir = fullfile(repoRoot, 'common');
addpath(commonDir);
dataRoot = preprocess_read_local_path(repoRoot, "data_root");
inputRoot = fullfile(dataRoot, dataset, 'mascot_dat');
protectedRoots = string({ ...
    inputRoot
    fullfile(dataRoot, dataset, 'preprocess_result')
    });

validateInputs(inputRoot, runNames);
ptmdecoderRoot = preprocess_read_local_path(repoRoot);
validateDependencies(ptmdecoderRoot, scriptDir);
ptmdecoderCommit = preprocess_check_ptmdecoder(ptmdecoderRoot);

addpath(ptmdecoderRoot);
addpath(fullfile(ptmdecoderRoot, 'FDR_control_generate_pep_spec_list'));

[candidateRoot, stageRoot] = preprocess_begin_candidate( ...
    candidateRoot, protectedRoots);
expectedFiles = expectedOutputFiles(runNames);
runSummaries = repmat(struct( ...
    'run', '', 'filtered_count', 0, 'selected_count', 0), ...
    numel(runNames), 1);
for idx = 1:numel(runNames)
    runName = runNames{idx};
    fprintf('[%d/%d] Rebuilding %s/%s\n', ...
        idx, numel(runNames), dataset, runName);
    runSummaries(idx) = processRun( ...
        fullfile(inputRoot, runName), stageRoot, runName, dataset);
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

function runSummary = processRun(inputDir, stageRoot, runName, dataset)
result = ReadDatResultFolder(inputDir);
if isempty(result)
    error([mfilename, ':EmptyDatResult'], ...
        'No Mascot PSMs were parsed from "%s".', inputDir);
end

[decoyType, groupType, ~, scores, numResult, indices] = JudgeGroup( ...
    result, 'Protein', 'REVERSE_', {'_HUMAN'});
[~, identified] = ComputeFDR(decoyType, groupType, scores, ...
    numResult, indices, 0.01);
fdrFilteredResult = result(identified.SF);

outputDir = fullfile(stageRoot, runName);
createDirectory(outputDir);
filteredPath = fullfile(outputDir, 'filtered_result_mascot.txt');
write_mascot_result_table(fdrFilteredResult, filteredPath);
[~, criteria] = syn_pepspec_dataset_config(dataset);
selectedEntries = select_syn_pepspec_entries( ...
    fdrFilteredResult, criteria);
write_peptide_spectra_list_file(selectedEntries, ...
    fullfile(outputDir, 'pepSpecFile.txt'));

runSummary = struct( ...
    'run', runName, ...
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

function validateDependencies(ptmdecoderRoot, scriptDir)
fdrDir = fullfile(ptmdecoderRoot, 'FDR_control_generate_pep_spec_list');
requiredFiles = {
    fullfile(fdrDir, 'ReadDatResultFolder.m')
    fullfile(fdrDir, 'JudgeGroup.m')
    fullfile(fdrDir, 'ComputeFDR.m')
    fullfile(fdrDir, 'write_mascot_result_table.m')
    fullfile(fdrDir, 'write_peptide_spectra_list_file.m')
    fullfile(ptmdecoderRoot, '@CFdrFilteredResultIO', 'write.m')
    fullfile(scriptDir, 'syn_pepspec_dataset_config.m')
    fullfile(scriptDir, 'select_syn_pepspec_entries.m')
    };
for idx = 1:numel(requiredFiles)
    if ~isfile(requiredFiles{idx})
        error([mfilename, ':MissingDependency'], ...
            'Missing required dependency: "%s".', requiredFiles{idx});
    end
end
end

function expectedFiles = expectedOutputFiles(runNames)
expectedFiles = strings(2 * numel(runNames), 1);
for idx = 1:numel(runNames)
    expectedFiles(2 * idx - 1) = fullfile( ...
        runNames{idx}, 'filtered_result_mascot.txt');
    expectedFiles(2 * idx) = fullfile(runNames{idx}, 'pepSpecFile.txt');
end
end

function createDirectory(directory)
[status, message] = mkdir(directory);
if ~status
    error([mfilename, ':CreateRunOutputFailed'], ...
        'Failed to create output directory "%s": %s', directory, message);
end
end
