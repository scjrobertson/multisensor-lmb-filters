function [r, W, V] = lmbGibbsSampling(associationMatrices, numberOfSamples)
% LMBGIBBSSAMPLING -- Determine posterior existence probabilities and association weights using a Gibbs sampler
%  [r, W, V] = lmbGibbsSampling(associationMatrices, numberOfSamples)
%
%   This function determines each object's posterior existence and marginal
%   association probabilities using Gibbs sampling. This function is a bit
%   more optimised for Matlab.
%
%   See also runLmbFilter, generateLmbAssociationMatrices,
%       computePosteriorLmbSpatialDistributions, lmbMurtysAlgorithm, 
%       loopyBeliefPropagation, lmbGibbsFrequencySampling
%
%   Inputs
%       associationMatrices - struct. A struct whose fields are the arrays required 
%           by the various data association algorithms.
%       numberOfSamples - double. The number of Gibbs samples we want to
%           generate.
%
%   Output
%       r - array. Each object's posterior existence probability.
%       W - array. An array of marginal association probabilities, where
%           each row is an object's marginal association probabilities.
%       V - array. An array of association events, where each row is
%           assignment of objects to measurements.

%% Declare variables
[n, m] = size(associationMatrices.P);
% Association vectors
[v, w] = initialiseGibbsAssociationVectors(associationMatrices.C);
V = zeros(numberOfSamples, n);
%% Gibbs sampling
for i = 1:numberOfSamples
    %% Generate a new Gibbs sample
    [v, w] = generateGibbsSample(associationMatrices.P, v, w);
    %% Store Gibbs sample
    V(i, :) = v;
end
%% Determine marginal distributions
% Distinct samples
V = unique(V, 'rows');
% Marginal distributions
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