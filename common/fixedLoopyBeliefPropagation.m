function [r, W] = fixedLoopyBeliefPropagation(associationMatrices, maximumNumberOfLbpIterations)
% FIXEDLOOPYBELIEFPROPAGATION -- Determine posterior existence probabilities and association weights using loopy belief propagation
%   [r, W] = fixedLoopyBeliefPropagation(associationMatrices, epsilon)
%
%   This function determines each object's posterior existence and marginal
%   association probabilities using loopy belief propagation (LBP).
%   This algorithm uses a fixed number of iteratrions, and it is only used
%   to verify the LMB filter's asymptotic computational complexity.
%
%   See also runLmbFilter, generateLmbAssociationMatrices,
%   computePosteriorLmbSpatialDistributions, lmbGibbsSampling,
%   lmbMurtysAlgorithm, lmbFilterTimeTrials, loopyBeliefPropagation
%
%   Inputs
%       associationMatrices - struct. A struct whose fields are the arrays required 
%           by the various data association algorithms.
%       maximumNumberOfLbpIterations - integer. Maximum allowable number
%           of LBP iterations.
%
%   Output
%       r - array. Each object's posterior existence probability.
%       W - array. An array of marginal association probabilities, where
%           each row is an object's marginal association probabilities.

%% Declare variables
SigmaMT = ones(size(associationMatrices.Psi));
%% Loopy belief propagation
for i = 1:maximumNumberOfLbpIterations
    % Pass messages from the object to the measurement clusters
    B = associationMatrices.Psi .* SigmaMT;
    SigmaTM = associationMatrices.Psi ./ (-B + sum(B, 2) + 1);
    % Pass messages from the measurement to the object clusters
    SigmaMT = 1./ (-SigmaTM + sum(SigmaTM, 1) + 1);
end
Gamma = [associationMatrices.phi B .* associationMatrices.eta];
q = sum(Gamma, 2);
%% Determine association probabilities
W = Gamma ./ q;
%% Determine existence probabilities
r = q ./ (associationMatrices.eta + q - associationMatrices.phi);
end