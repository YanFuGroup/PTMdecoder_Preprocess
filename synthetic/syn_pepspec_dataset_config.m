function [config, criteria] = syn_pepspec_dataset_config(dataset)
%SYN_PEPSPEC_DATASET_CONFIG Return paper-specific synthetic dataset rules.

dataset = char(dataset);
switch dataset
    case 'syn_ind'
        config = struct('runPrefix', 'ind', 'runCount', 16);
    case 'syn_mix'
        config = struct('runPrefix', 'mix', 'runCount', 18);
    otherwise
        error('prepare_syn_pepspec:InvalidDataset', ...
            'Dataset must be syn_ind or syn_mix. Got: %s', dataset);
end

criteria = struct( ...
    'targetPeptide', 'GKGGKGLGKGGAKR', ...
    'neutralMass', 1535.878357, ...
    'massTolerance', 0.0001, ...
    'charge', 2, ...
    'minimumMatchedIons', []);
if strcmp(dataset, 'syn_ind')
    criteria.minimumMatchedIons = 8;
end
end
