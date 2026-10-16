function selectedEntries = select_dna_histone_pepspec_entries(fdrFilteredResult)
%SELECT_DNA_HISTONE_PEPSPEC_ENTRIES Apply the legacy post-FDR policy.
%
%   selectedEntries = select_dna_histone_pepspec_entries(fdrFilteredResult)
%   keeps only legacy-unique dataset/scan keys, excludes unmodified,
%   Oxidation-containing, and label-only PSMs, then sorts by peptide text.
%   The input struct array is not modified.

selectedMask = select_unique_dataset_scan(fdrFilteredResult);
uniqueEntries = fdrFilteredResult(selectedMask);

modifications = {uniqueEntries.modification}';
isNotUnmodified = ~strcmp(modifications, '-');
isNotOxidation = ~contains(modifications, 'Oxidation');
candidateMask = isNotUnmodified & isNotOxidation;

labelOnlyMask = false(size(modifications));
labelOnlyModifications = { ...
    'Label:13C(6) (K)', ...
    'Label:13C(6)15N(4) (R)'};
labelOnlyMask(candidateMask) = cellfun( ...
    @(text) all(ismember(strsplit(text, ','), labelOnlyModifications)), ...
    modifications(candidateMask));

selectedEntries = uniqueEntries(candidateMask & ~labelOnlyMask);
[~, order] = sort({selectedEntries.peptide});
selectedEntries = selectedEntries(order);
end
