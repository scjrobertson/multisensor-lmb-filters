function model = generateSimplifiedModel(numberOfObjects, detectionProbability, clutterRate)
% GENERATESIMPLIFIEDMODEL -- Generates a structure containing all simulation info.
%   model = generateSimplifiedModel(numberOfObjects, detectionProbability, clutterRate);
%
%   Declares all simulation information.
%
%   See also generateAssociationMatrices.
%
%   Inputs
%       numberOfObjects - integer. The number of objects for a simulation.
%       clutterRate - integer. The number of expected clutter per
%           time-step. Clutter is assumed to generated by a Poisson point
%           process.
%       detectionProbability - double. The probability of the sensor
%           detecting a target-generated measurement.
%
%   Output
%       model - struct. A struct with the fields declared in this function.

%% State space
model.xDim = 2;
model.numberOfObjects = numberOfObjects;
model.r = rand(model.numberOfObjects, 1); %0.95 * ones(numberOfObjects, 1);
model.stateSpaceLimits = 2 * [0 1; 0 1];
model.mu = generateRandomMeans(model.xDim, model.numberOfObjects, model.stateSpaceLimits);
model.Sigma = generateRandomCovarianceMatrices(model.xDim, model.numberOfObjects);
%% Observation space
% Measurement model
model.zDim = 2;
model.detectionProbability = detectionProbability;
model.C = [eye(model.zDim) zeros(model.zDim, max(model.xDim - model.zDim, 0))];
model.Q = (2^2) * eye(model.zDim);
% Clutter
model.observationSpaceLimits = 2 * [0 1; 0 1];
expectedNumberOfClutterReturns = clutterRate + (clutterRate == 0);
model.observationSpaceVolume = prod(model.observationSpaceLimits(:, 2) - model.observationSpaceLimits(:, 1));
model.clutterDensity = expectedNumberOfClutterReturns / model.observationSpaceVolume;
model.expectedNumberOfClutterReturns = clutterRate;
%% Gibbs sampling 
model.numberOfSamples = 1e6;
%% Loopy belief propagation
model.maximumNumberOfLbpIterations = 1e4;
model.lbpCovergenceTolerance = 1e-18;
%% Murty's algorithm
model.numberOfAssignments = 1e4;
%% Generate random mean
function mu = generateRandomMeans(d, numberOfGmComponents, spatialLimits)
mu = cell(1, numberOfGmComponents);
for i = 1:numberOfGmComponents
    mu{i} = spatialLimits(:, 1) + 2 * spatialLimits(:, 2) .* rand(1, d)';
end
end
%% Generate random covariance
function Sigma = generateRandomCovarianceMatrices(d, numberOfGmComponents)
Sigma = cell(1, numberOfGmComponents);
for i = 1:numberOfGmComponents
    Q = randn(d, d);
    S = diag(abs( 5^2 + randn(d,1) ) );
    Sigma{i} = Q' * S * Q;
end
end
end