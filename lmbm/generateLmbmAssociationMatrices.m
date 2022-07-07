function [associationMatrices, posteriorParameters] = generateLmbmAssociationMatrices(hypothesis, z, model)
% GENERATELMBMASSOCIATIONMATRICES -- Compute the association matrices required by the data association algorithms.
%   [gibbsParameters, posteriorParameters] = generateLmbGibbsMatrices(objects, z, model)
%
%   This function computes the association matrices required by the LBP.
%   Gibbs sampler, and Murty's algorithms. It also determines the measurement-updated components that
%   are used to determine each object's posterior spatial distribution.
%
%   See also runLmbmFilter, generateModel, lmbmGibbsSampling,
%
%   Inputs
%       hypothesis - struct. A struct containing the prior LMBM hypothesis' Bernoulli components.
%       z - cell array. A cell array of measurements for the
%           current time-step.
%       model - struct. A struct with the fields declared in generateModel.
%
%   Output
%       associationMatrices - struct. A struct whose fields are the arrays required 
%           by the various data association algorithms.
%       posteriorParameters - struct. A struct whose fields are a
%           hypothesis posterior distribution parameters.

%% Declare output structs
numberOfObjects = numel(hypothesis.r);
numberOfMeasurements = numel(z);
% Log likelihood matrix
R = zeros(numberOfObjects, numberOfMeasurements);
% Auxiliary variables
phi = (1 - model.detectionProbability) * hypothesis.r;
eta = 1 - model.detectionProbability * hypothesis.r;
% Updated components for the objects' posterior spatial distributions
posteriorParameters.r = phi ./ eta;
posteriorParameters.mu = cell(numberOfObjects, numberOfMeasurements + 1);
posteriorParameters.Sigma = hypothesis.Sigma;
%% Populate the Gibbs sampler's likelihood array, and compute posterior components
for i = 1:numberOfObjects
    % Missed detection event
    posteriorParameters.mu{i, 1} = hypothesis.mu{i};
    % Determine posterior parameters for each object's spatial distribution
    muZ = model.C * hypothesis.mu{i};
    Z = model.C * hypothesis.Sigma{i} * model.C' + model.Q;
    logGaussianNormalisingConstant = - (0.5 * model.zDimension) * log(2 * pi) - 0.5 * log(det(Z));
    logLikelihoodRatioTerms = log(hypothesis.r(i)) + log(model.detectionProbability) - log(model.clutterPerUnitVolume);
    ZInv = inv(Z);
    K = hypothesis.Sigma{i} * model.C' * ZInv;
    posteriorParameters.Sigma{i} = (eye(model.xDimension) - K * model.C) * hypothesis.Sigma{i};
    % Determine total marginal likelihood, and determine posterior components
    for j = 1:numberOfMeasurements
        % Determine marginal likelihood ratio
        nu = z{j} - muZ;
        gaussianLogLikelihood = logGaussianNormalisingConstant - 0.5 * nu' * ZInv * nu;
        R(i, j) = logLikelihoodRatioTerms + gaussianLogLikelihood;
        % Determine posterior mean for each measurement
        posteriorParameters.mu{i, j+1} = hypothesis.mu{i} + K * nu;
    end
end
%% Determine Gibbs sampler parameters
RLinear = exp(R);
% Gibbs sampler association matrices
associationMatrices.P = RLinear ./ (RLinear + eta);
associationMatrices.L = [log(eta) R];
% Murty's algorithm cost matrix
associationMatrices.C = -R;
end