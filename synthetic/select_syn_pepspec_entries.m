function selectedEntries = select_syn_pepspec_entries(entries, criteria)
%SELECT_SYN_PEPSPEC_ENTRIES Apply the paper-specific synthetic PSM rules.

validateCriteria(criteria);
if isempty(entries)
    selectedEntries = entries;
    return;
end

datasetNames = {entries.DatasetName}';
scanNames = {entries.Scan}';
groupIndex = findgroups(datasetNames, scanNames);
groupCounts = accumarray(groupIndex, 1);
isUniqueSpectrum = groupCounts(groupIndex) == 1;

peptides = {entries.peptide}';
charges = numericColumn({entries.Charge}');
neutralMasses = numericColumn({entries.Calc_neutral_pepmass}');
isSelected = isUniqueSpectrum ...
    & strcmp(peptides, criteria.targetPeptide) ...
    & charges == criteria.charge ...
    & abs(neutralMasses - criteria.neutralMass) < criteria.massTolerance;

if ~isempty(criteria.minimumMatchedIons)
    matchedIonCounts = numericColumn({entries.num_match_ions}');
    isSelected = isSelected & matchedIonCounts >= criteria.minimumMatchedIons;
end

selectedEntries = entries(isSelected);
[~, sortedIndex] = sort({selectedEntries.peptide});
selectedEntries = selectedEntries(sortedIndex);
end

function validateCriteria(criteria)
requiredFields = {'targetPeptide', 'neutralMass', 'massTolerance', ...
    'charge', 'minimumMatchedIons'};
if ~isstruct(criteria) || ~all(isfield(criteria, requiredFields))
    error('select_syn_pepspec_entries:InvalidCriteria', ...
        'Criteria must contain target peptide, mass, tolerance, charge, and matched-ion policy.');
end
end

function values = numericColumn(rawValues)
values = cellfun(@toDouble, rawValues);
end

function value = toDouble(rawValue)
if isnumeric(rawValue)
    value = double(rawValue);
else
    value = str2double(char(rawValue));
end
end
