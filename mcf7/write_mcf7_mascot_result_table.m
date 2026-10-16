function write_mcf7_mascot_result_table(result, filename)
%WRITE_MCF7_MASCOT_RESULT_TABLE Preserve the legacy signed-zero spelling.
%
% The PTMdecoder table writer uses num2str, which drops a negative zero's
% sign. Format only zero massdiff values with the original fprintf precision
% before delegating. MATLAB passes the struct by value, so the PSM values
% used for selection remain unchanged.

for idx = 1:numel(result)
    value = result(idx).massdiff;
    if isnumeric(value) && isscalar(value) && value == 0
        result(idx).massdiff = sprintf('%.6f', value);
    end
end
write_mascot_result_table(result, filename);
end
