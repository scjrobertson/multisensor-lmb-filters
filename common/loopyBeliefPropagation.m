function [r, W] = loopyBeliefPropagation(associationMatrices, epsilon, maximumNumberOfLbpIterations)
% LOOPYBELIEFPROPAGATION -- Determine posterior existence probabilities and association weights using loopy belief propagation
%   [r, W] = loopyBeliefPropagation(associationMatrices, epsilon, maximumNumberOfLbpIterations)
%
%   This function determines each object's posterior existence and marginal
%   association probabilities using loopy belief propagation (LBP).
%
%   See also runLmbFilter, generateLmbAssociationMatrices,
%   computePosteriorLmbSpatialDistributions, lmbGibbsSampling,
%   lmbMurtysAlgorithm
%
%   Inputs
%       associationMatrices - struct. A struct whose fields are the arrays required 
%           by the various data association algorithms.
%       epsilon - double. The convergence tolerance for the LBP aglorithm.
%       maximumNumberOfLbpIterations - integer. Maximum allowable number
%           of LBP iterations.
%
%   Output
%       r - (n, 1) array. Each object's posterior existence probability.
%       W - (n, m) array. An array of marginal association probabilities, where
%           each row is an object's marginal association probabilities.

%% Declare variables
SigmaMT = ones(size(associationMatrices.Psi));
notConverged = true;
counter = 0;
%% Loopy belief propagation
while notConverged
    % Cache previous iteration's messages
    SigmaMTOld = SigmaMT;
    % Pass messages from the object to the measurement clusters
    B = associationMatrices.Psi .* SigmaMT;
    SigmaTM = associationMatrices.Psi ./ (-B + sum(B, 2) + 1);
    % Pass messages from the measurement to the object clusters
    SigmaMT = 1./ (-SigmaTM + sum(SigmaTM, 1) + 1);
    % Check for convergence
    counter = counter + 1;
    delta = abs(SigmaMT - SigmaMTOld);
    notConverged = (max(delta(:)) > epsilon) && (counter < maximumNumberOfLbpIterations);
end
Gamma = [associationMatrices.phi B .* associationMatrices.eta];
q = sum(Gamma, 2);
%% Determine association probabilities
W = Gamma ./ q;
%% Determine existence probabilities
r = q ./ (associationMatrices.eta + q - associationMatrices.phi);
end