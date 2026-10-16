function selectedMask = select_unique_dataset_scan(result)
%SELECT_UNIQUE_DATASET_SCAN Legacy-compatible uniqueness mask.
%
%   selectedMask = select_unique_dataset_scan(result) returns a
%   logical mask for an FDR-filtered PSM struct array. The key is the exact
%   text produced by the legacy expression:
%
%       strcat({result.DatasetName}', ' ', {result.Scan}')
%
%   MATLAB strcat removes trailing whitespace from character-array inputs;
%   in particular, the single-space argument becomes empty. The resulting key
%   is therefore normally DatasetName text directly followed by Scan text,
%   not a true (DatasetName, Scan) tuple and not a numeric Scan comparison.
%
%   A PSM is selected only when its legacy text key occurs exactly once in
%   the complete input. Every PSM sharing a key that occurs more than once is
%   rejected. The helper never chooses a representative PSM and never
%   compares peptide, score, mass, or modification fields.
%
%   INPUT SCHEMA
%   result must be a 1-by-N or N-by-1 struct array. Every element must contain
%   DatasetName and Scan as char row vectors. Other fields may be present and
%   are ignored. Values are compared using the legacy strcat construction,
%   without trimming, numeric conversion, case folding, or other content
%   normalization. Scan does not need to be a decimal integer.
%
%   OUTPUT
%   selectedMask is a logical array with the same shape as result. The mask
%   follows input order; this function does not sort, filter the struct array
%   itself, or modify any input field.
%
%   ERRORS
%   select_unique_dataset_scan:InvalidResult       Invalid type/shape/fields
%   select_unique_dataset_scan:InvalidDatasetName  Invalid DatasetName value
%   select_unique_dataset_scan:InvalidScan         Invalid Scan value
%
%   FICTIONAL EXAMPLE
%   DatasetName values {'A1', 'A'} with Scan values {'23', '123'} both create
%   the legacy key 'A123', so both PSMs are dropped. Under DatasetName 'A',
%   Scan values {'001'} and {'1'} create distinct text keys and both remain.
%
%   See also README.md (Preprocessing rules) for dataset selection policies.

resultSize = size(result);
% size returns every nonsingleton trailing dimension, so reject N-D inputs
% before testing for the supported 1-by-N or N-by-1 shapes.
isSupportedShape = numel(resultSize) == 2 && ...
    (resultSize(1) == 1 || resultSize(2) == 1);
if ~isstruct(result) || ~isSupportedShape || ...
        ~isfield(result, 'DatasetName') || ~isfield(result, 'Scan')
    error('select_unique_dataset_scan:InvalidResult', ...
        'Result must be a 1-by-N or N-by-1 struct array with DatasetName and Scan fields.');
end

recordCount = numel(result);
if recordCount == 0
    selectedMask = false(resultSize);
    return;
end

datasetNames = reshape({result.DatasetName}, recordCount, 1);
validDatasetName = cellfun(@ischar, datasetNames) & ...
    cellfun(@isrow, datasetNames);
if ~all(validDatasetName)
    error('select_unique_dataset_scan:InvalidDatasetName', ...
        'DatasetName must be a nonempty char row vector for every PSM.');
end
blankDatasetName = cellfun(@isempty, strtrim(datasetNames));
if any(blankDatasetName)
    error('select_unique_dataset_scan:InvalidDatasetName', ...
        'DatasetName must contain a non-whitespace character for every PSM.');
end

scanText = reshape({result.Scan}, recordCount, 1);
validScanType = cellfun(@ischar, scanText) & cellfun(@isrow, scanText);
if ~all(validScanType) || any(cellfun(@isempty, scanText))
    error('select_unique_dataset_scan:InvalidScan', ...
        'Scan must be a nonempty char row vector for every PSM.');
end

% Reproduce the historical key verbatim. Do not replace this expression with
% DatasetName/Scan concatenation that differs from MATLAB strcat semantics.
legacyKeys = strcat(datasetNames, ' ', scanText);
[~, ~, groupIndex] = unique(legacyKeys, 'stable');
groupCounts = accumarray(groupIndex, 1);
selectedMask = reshape(groupCounts(groupIndex) == 1, resultSize);
end
