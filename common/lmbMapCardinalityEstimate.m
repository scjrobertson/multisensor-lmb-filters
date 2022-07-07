function [nMap, mapIndices] = lmbMapCardinalityEstimate(r)
% LMBMAPCARDINALITYESTIMATE -- Determine approximate LMB MAP estimate
%   [nMap, mapIndices] = lmbMapCardinalityEstimate(r)
%
%   This function computes an approximate MAP estimate for the LMB filter.
%
%   See also runLmbFilter.
%
%   Inputs
%       r - array. Each object's posterior existence probability.
%
%   Output
%       nMap - integer. The MAP estimate for the LMB
%           cardinality estimate.
%       mapIndices - array. The indices of the nMap greatest indices
%           of r

% Use Mahler's aglorithm to determine the LMB cardinality distribution
r = r - 1e-6; % Does not work with unit existence probabilities
rho = prod(1 - r)*esf(r./(1-r));
% Determine the MAP estimate of the distribution
[~, maxCardinalityIndex] = max(rho);
% The MAP estimate cannot be larger than the number of objects
nMap = min(maxCardinalityIndex - 1, length(r));
% Sort r in descending order
[~, sortedIndices] = sort(-r);
% Choose the nMap largest indices of r
mapIndices = sortedIndices(1:nMap);
end