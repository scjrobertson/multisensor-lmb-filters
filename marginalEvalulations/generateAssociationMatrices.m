function associationMatrices = generateAssociationMatrices(model)
% GENERATEASSOCIATIONMATRICES -- Compute the association matrices required for the data association algorithms
%   [associationMatrices, posteriorParameters] = generateAssociationMatrices(objects, z, model)
%
%   This function computes the association matrices required by the LBP,
%   Gibbs sampler, and Murty's algorithms. It simulates a basic, but
%   difficult, data association problem, and then determines the relevant
%   association matrices.
%
%   See also generateSimplifiedModel
%
%   Inputs
%       model - struct. A struct with the fields declared in generateModel.
%
%   Output
%       associationMatrices - struct. A struct whose fields are the arrays required 
%           by the various data association algorithms.

%% Generate measurements
objectsWhichGeneratedMeasurements = rand(1, model.numberOfObjects) < 1; %model.detectionProbability;
numberOfObjectGeneratedMeasurements = sum(objectsWhichGeneratedMeasurements);
numberOfClutterMeasurements = 0;
if (model.expectedNumberOfClutterReturns > 0)
    numberOfClutterMeasurements = poissrnd(model.expectedNumberOfClutterReturns);
end
numberOfMeasurements = numberOfObjectGeneratedMeasurements + numberOfClutterMeasurements;
measurements = cell(1, numberOfMeasurements);
% Object generated measurements
QChol = chol(model.Q, 'lower');
counter = 0;
for i = 1:model.numberOfObjects
    if (objectsWhichGeneratedMeasurements(i))
        counter = counter + 1;
        measurements{counter} = model.C * model.mu{i} + QChol * randn(1, model.zDim)';
    end
end
% Clutter measurements
for i = 1:numberOfClutterMeasurements
    measurements{counter + i} = model.observationSpaceLimits(:, 1) + 2 * model.observationSpaceLimits(:, 2) .* rand(1, model.zDim)';
end
% Shuffle
measurements = measurements(randperm(numberOfMeasurements));
%% Generate Association matrices
L = zeros(model.numberOfObjects, numberOfMeasurements);
for i = 1:model.numberOfObjects
   % Determine marginal distribution
   mu = model.C * model.mu{i};
   K = inv(model.C * model.Sigma{i} * model.C' + model.Q);
   normalisingConstant = sqrt( det(K / (2 * pi) ) );
   % Determine likelihood
   for j = 1:numberOfMeasurements
       nu = measurements{j} - mu;
       L(i, j) = model.r(i) * model.detectionProbability *  normalisingConstant * exp( -0.5 * ( nu' * K * nu ) ) / model.clutterDensity;
   end
end
% LBP association matrices
associationMatrices.eta = 1 - model.detectionProbability * model.r;
associationMatrices.phi = (1 - model.detectionProbability) * model.r;
associationMatrices.Psi = L ./ associationMatrices.eta;
% Gibbs matrices
associationMatrices.P = L ./ (associationMatrices.eta + L);
associationMatrices.L = [associationMatrices.eta L];
associationMatrices.R = [associationMatrices.phi ./ associationMatrices.eta ones(model.numberOfObjects, numberOfMeasurements)];
% Murty's algorithm cost matrix
associationMatrices.C = -log(L);
end