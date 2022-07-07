function A = multisensorLmbmGibbsSampling(L, numberOfSamples)
% MULTISENSORLMBMGIBBSSAMPLING -- Generate association events using a multi-sensor Gibbs sampler.
%   A = multisensorLmbmGibbsSampling(L, numberOfSamples)
%
%   This function generates a set of posterior hypotheses for a given prior
%   hypothesis' input matrix using Gibbs sampling. 
%
%   See also runMultisensorLmbmFilter, generateMultisensorLmbmAssociationMatrices, generateMultisensorAssociationEvent
%
%   Inputs
%       L - (m1 + 1, m2 + 1, ..., ms + 1, n) array. A log likelihood
%           matrix used for Gibbs sampling and evaluating the probability
%           density of an association event.
%       numberOfSamples - double. The number of Gibbs samples we want to
%           generate.
%
%   Output
%       A - array. An array of distinct association events, where each row of the
%           array is an association event. Each assoication event is a
%           matrix, but it is flattened here.

%% Declare variables
% Sizes
ell = size(L);
m = ell(1:end-1) - 1;
n = ell(end);
numberOfSensors = length(m);
% Association vectors
V = zeros(n, numberOfSensors);
W = zeros(max(m), numberOfSensors);
A = zeros(numberOfSamples, n * numberOfSensors);
%% Gibbs sampling
for i = 1:numberOfSamples
    %% Generate a new Gibbs sample
    [V, W] = generateMultisensorAssociationEvent(L, V, W);
    %% Store Gibbs sample
    A(i, :) = reshape(V, 1, n * numberOfSensors);
end
%% Keep only the distinct samples
A = unique(A, 'rows');
end