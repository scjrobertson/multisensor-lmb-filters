function [r, W, V] = lmbMurtysAlgorithm(associationMatrices, numberOfAssignments)
% LMBMURTYSALGORITHM -- Determine posterior existence probabilities and association weights using Murty's algorithm
%   [r, W, V] = lmbMurtysAlgorithm(associationMatrices, numberOfAssignments)
%
%   This function determines each object's posterior existence and marginal
%   association probabilities using Murty's algorithm. This function uses
%   Vo et al.'s code.
%
%   See also runLmbFilter, generateLmbAssociationMatrices,
%   computePosteriorLmbSpatialDistributions, lmbGibbsSampling,
%   loopyBeliefPropagation
%
%   Inputs
%       associationMatrices - struct. A struct whose fields are the arrays required 
%           by the various data association algorithms.
%       numberOfAssignments - double. The number of association events we
%           want Murty's algorithm to generate
%
%   Output
%       r - array. Each object's posterior existence probability.
%       W - array. An array of marginal association probabilities, where
%           each row is an object's marginal association probabilities.
%       V - array. An array of association events, where each row is
%           assignment of objects to measurements.

[n, m] = size(associationMatrices.C);
%% Determine most likely assignments using Vo et al.'s implementation of Murty's algorithm
V = murtysAlgorithmWrapper(associationMatrices.C, numberOfAssignments);
%% Determine marginal distributions
W = repmat(V, 1, 1, m+1) == reshape(0:m, 1, 1, m+1);
J = reshape(associationMatrices.L(n * V + (1:n)), size(V, 1), n);
L = permute(sum(prod(J, 2) .* W, 1), [2 1 3]);
Sigma = reshape(L, n, m+1);
% Normalise
Tau = (Sigma .* associationMatrices.R) ./ sum(Sigma, 2);
%% Determine existence probabilities
r = sum(Tau, 2);
%% Determine marginal association probabilities
W =  Tau ./ r;
end