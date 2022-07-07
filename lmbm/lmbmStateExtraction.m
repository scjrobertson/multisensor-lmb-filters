function [cardinalityEstimate, extractionIndices] = lmbmStateExtraction(hypotheses, useEapOnLmbm)
% LMBMSTATEEXTRACTION -- Heuristically determine the number of objects
% present, and their respective parameters' indices.
%    [cardinalityEstimate, extractionIndices] = lmbmStateExtraction(hypotheses, useEapOnLmbm)
%
%   This function computes an approximate cardinality estimate for the LMBM filter.
%
%   See also runLmbmFilter
%
%   Inputs
%       hypotheses - struct. A struct containing posterior LMBM
%           hypotheses, but with unnormalised hypothesis weights.
%       useEapOnLmbm - bool. A 'true' if we want to use an approximate EAP for 
%           state extraction, a 'false' we want to use an heuristic MAP.
%
%   Output
%       cardinalityEstimate - integer. An estimate of the number of objects present.
%       extractionIndices - array. An array of indices of the objects' with
%           the greatest existence probabilities.

rTotal = sum([hypotheses.w] .* [hypotheses.r], 2);
if (useEapOnLmbm)
    %% Heuristic EAP estimate
    % Determine expected cardinality
    cardinalityEstimate = floor(sum(rTotal));
    [~,  extractionIndices] = maxk(hypotheses(1).r, cardinalityEstimate);
else
    %% Very heuristic MAP estimate
    % Determine MAP cardinality estimate of the most likely hypothesis
    [cardinalityEstimate, extractionIndices] = lmbMapCardinalityEstimate([rTotal]);
end

end