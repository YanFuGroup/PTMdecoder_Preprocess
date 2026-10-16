function selectedMask = select_mcf7_modified_peptides(result, experimentName)
%SELECT_MCF7_MODIFIED_PEPTIDES Apply the explicit MCF7 selector.
%
%   selectedMask = select_mcf7_modified_peptides(result, experimentName)
%   returns a mask for the post-FDR and uniqueness-filtered result. The
%   current JIB04 and SAHA rules are intentionally spelled out as separate
%   branches: both keep modified peptides and reject unmodified entries or
%   entries whose modification text contains Oxidation. Keeping the branches
%   explicit avoids hiding a future dataset-specific change behind generic
%   configuration.
%
%   INPUT
%       result          1-by-N or N-by-1 PSM struct array.
%       experimentName  'MCF7_JIB04' or 'MCF7_SAHA'.
%
%   OUTPUT
%       selectedMask    Logical array with the same shape as result.

if ~ischar(experimentName) && ~isstring(experimentName)
    error('select_mcf7_modified_peptides:InvalidExperiment', ...
        'Experiment name must be char or string.');
end


resultSize = size(result);
isSupportedShape = numel(resultSize) == 2 && ...
    (resultSize(1) == 1 || resultSize(2) == 1);
if ~isSupportedShape || ~isstruct(result) || ...
        ~isfield(result, 'modification')
    error('select_mcf7_modified_peptides:InvalidResult', ...
        'Result must be a 1-by-N or N-by-1 struct array with modification.');
end

recordCount = numel(result);
modifications = reshape({result.modification}, recordCount, 1);
if ~all(cellfun(@ischar, modifications) & cellfun(@isrow, modifications))
    error('select_mcf7_modified_peptides:InvalidModification', ...
        'Every PSM modification must be a char row vector.');
end

% Keep the two experiment cases explicit. They are currently the same
% legacy-compatible scientific rule, but a future rule change must be made
% and reviewed for one dataset without silently changing the other.
switch lower(char(experimentName))
    case 'mcf7_jib04'
        selected = ~strcmp(modifications, '-') & ...
            ~contains(modifications, 'Oxidation');
    case 'mcf7_saha'
        selected = ~strcmp(modifications, '-') & ...
            ~contains(modifications, 'Oxidation');
    otherwise
        error('select_mcf7_modified_peptides:UnknownExperiment', ...
            'Unsupported MCF7 experiment: "%s".', char(experimentName));
end

selectedMask = reshape(selected, resultSize);
end
