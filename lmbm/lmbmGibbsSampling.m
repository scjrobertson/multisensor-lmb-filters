function V = lmbmGibbsSampling(P, C, numberOfSamples)
% LMBMGIBBSSAMPLING -- Generate association events using a Gibbs sampler.
%  [r, W] = lmbmGibbsSampling(gibbsParameters, numberOfSamples)
%
%   This function generates a set of posterior hypotheses for a given prior
%   hypothesis' input matrix using Gibbs sampling. 
%
%   See also runLmbmFilter, generateLmbmAssociationMatrices, lmbGibbsSampling
%
%   Inputs
%       P - (n, m) array. An array of sampling probabilites for the Gibbs
%           sampler. See also generateLbmmGibbsMatrices.
%       C - (n, m) array. The cost matrix declared for Murty's algorithm in 
%           generateLmbAssociationMatrices.
%       numberOfSamples - double. The number of Gibbs samples we want to
%           generate.
%
%   Output
%       V - array. An array of distinct association events, where each row of the
%           array is an association event.

%% Declare variables
n = size(P, 1);
% Association vectors
[v, w] = initialiseGibbsAssociationVectors(C);
V = zeros(numberOfSamples, n);
%% Gibbs sampling
for i = 1:numberOfSamples
    %% Generate a new Gibbs sample
    [v, w] = generateGibbsSample(P, v, w);
    %% Store Gibbs sample
    V(i, :) = v;
end
%% Keep only the distinct samples
V = unique(V, 'rows');
end