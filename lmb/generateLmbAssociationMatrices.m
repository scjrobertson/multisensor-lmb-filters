function [associationMatrices, posteriorParameters] = generateLmbAssociationMatrices(objects, z, model)
% GENERATELMBASSOCIATIONMATRICES -- Compute the association matrices required for the data association algorithms
%   [associationMatrices, posteriorParameters] = generateLmbAssociationMatrices(objects, z, model)
%
%   This function computes the association matrices required by the LBP,
%   Gibbs sampler, and Murty's algorithms. It also determines the measurement-updated components that
%   are used to determine each object's posterior spatial distribution.
%
%   See also runLmbFilter, generateModel, loopyBeliefPropagation, lmbGibbsSampling, lmbMurtysAlgorithm
%
%   Inputs
%       objects - struct. A struct containing the prior LMB's Bernoulli components.
%       z - cell array. A cell array of measurements for the
%           current time-step.
%       model - struct. A struct with the fields declared in generateModel.
%
%   Output
%       associationMatrices - struct. A struct whose fields are the arrays required 
%           by the various data association algorithms.
%       posteriorParameters - struct. A struct whose fields are an object's
%           posterior spatial distribution parameters.
%

%% Declare output structs
numberOfObjects = numel(objects);
numberOfMeasurements = numel(z);
% Auxillary matrices
L = zeros(numberOfObjects, numberOfMeasurements);
phi = zeros(numberOfObjects, 1);
eta = zeros(numberOfObjects, 1);
% Updated components for the objects' posterior spatial distributions
posteriorParameters.w = [];
posteriorParameters.mu = {};
posteriorParameters.Sigma = {};
posteriorParameters = repmat(posteriorParameters, 1, numberOfObjects);
%% Populate the LBP arrays, and compute posterior components
for i = 1:numberOfObjects
    % Predeclare the object's posterior components, and include missed detection event
    posteriorParameters(i).w = repmat(log(objects(i).w * (1 - model.detectionProbability)), numberOfMeasurements + 1, 1);
    posteriorParameters(i).mu  = repmat(objects(i).mu, numberOfMeasurements + 1, 1);
    posteriorParameters(i).Sigma = repmat(objects(i).Sigma, numberOfMeasurements + 1, 1);
    % Populate auxiliary LBP parameters
    phi(i) = (1 -  model.detectionProbability) * objects(i).r;
    eta(i) = 1 - model.detectionProbability * objects(i).r;
    %% Determine marginal likelihood ratio of the object generating each measurement
    for j = 1:objects(i).numberOfGmComponents
        % Update components for a mixture component
        muZ = model.C * objects(i).mu{j};
        Z = model.C * objects(i).Sigma{j} * model.C' + model.Q;
        logGaussianNormalisingConstant = - (0.5 * model.zDimension) * log(2 * pi) - 0.5 * log(det(Z));
        logLikelihoodRatioTerms = log(objects(i).r) + log(model.detectionProbability) + log(objects(i).w(j)) - log(model.clutterPerUnitVolume);
        ZInv = inv(Z);
        K = objects(i).Sigma{j} * model.C' * ZInv;
        SigmaUpdated = (eye(model.xDimension) - K * model.C) * objects(i).Sigma{j};
        % Determine total marginal likelihood, and determine posterior components
        for k = 1:numberOfMeasurements
            % Determine marginal likelihood ratio
            nu = z{k} - muZ;
            gaussianLogLikelihood = logGaussianNormalisingConstant - 0.5 * nu' * ZInv * nu;
            L(i, k) = L(i, k) + exp(logLikelihoodRatioTerms + gaussianLogLikelihood);
            % Determine updated mean and covariance for each mixture component
            posteriorParameters(i).w(k+1, j) = log(objects(i).w(j)) + gaussianLogLikelihood + log(model.detectionProbability) - log(model.clutterPerUnitVolume);
            posteriorParameters(i).mu{k+1, j} = objects(i).mu{j} + K * nu;
            posteriorParameters(i).Sigma{k+1, j} = SigmaUpdated;
        end
    end
    % Normalise weights
    maximumWeights = max(posteriorParameters(i).w, [], 2);
    offsetWeights = posteriorParameters(i).w - maximumWeights;
    posteriorParameters(i).w = exp(offsetWeights) ./ sum(exp(offsetWeights), 2);
end
%% Output association matrices
% LBP association matrices
associationMatrices.Psi = L ./ eta;
associationMatrices.phi = phi;
associationMatrices.eta = eta;
% Gibbs sampler association matrices
associationMatrices.P = L./ (L + eta);
associationMatrices.L = [eta L];
associationMatrices.R = [(phi ./ eta) ones(numberOfObjects, numberOfMeasurements)];
% Murty's algorithm association matrices
associationMatrices.C = -log(L);
end