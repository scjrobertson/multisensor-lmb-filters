function [v, w] = initialiseGibbsAssociationVectors(C)
% INITIALISEGIBBSASSOCIATIONVECTORS -- Use Murty's algorithm to initialise
% the Gibbs association vectors
%   [v, w] = initialiseGibbsAssociationVectors(C)
%
%   This function determines the most likely association event, and uses it
%   to intialise the Gibbs samplers association vectors.
%
%   See also runLmbFilter, generateLmbAssociationMatrices,  lmbGibbsSampling
%
%   Inputs
%       C - (n, m) array. The cost matrix declared for Murty's algorithm in 
%           generateLmbAssociationMatrices.
%
%   Output
%       v - (n, 1) array. The most likely object-to-measurement association
%           vector.
%       w - (1, m) array. The corresponding most likly
%           measurement-to-object association vector.


[n, m] = size(C);
%% Determine most likely assignments using Vo et al.'s implementation of Murty's algorithm
v = murtysAlgorithmWrapper(C, 1)';
%% Determine w
w = zeros(1, m);
ell = v > 0;
p = 1:n;
w(v(ell)) = p(ell);
end