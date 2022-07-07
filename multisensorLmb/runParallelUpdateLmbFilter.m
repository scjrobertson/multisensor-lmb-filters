function stateEstimates = runParallelUpdateLmbFilter(model, measurements)
% RUNPARALLELUPDATELMBFILTER -- Run a multi-sensor LMB filter that uses a parallel measurment update.
%   stateEstimates = runParallelUpdateLmbFilter(model, measurements)
%
%   Run a multi-sensor LMB filter that uses a parallel measurment update.
%   Measurement update variants include arithmetic average (AA), geometric
%   average (GA), and parallel update (PU) update.
%
%   See also generateMultisensorModel, generateMultisensorGroundTruth, lmbPredictionStep,
%   generateLmbSensorAssociationMatrices, loopyBeliefPropagation, lmbGibbsSampling, 
%   lmbMurtysAlgorithm, computePosteriorLmbSpatialDistributions,
%   lmbMapCardinalityEstimate, aaLmbTrackMerging, gaLmbTrackMerging,
%   puLmbTrackMerging.
%
%   Inputs
%       model - struct. A struct with the fields declared in generateModel.
%       measurements - cell array. An array containing the measurements for
%           each time-step of the simulation. See also generateModel.
%
%   Output
%       stateEstimates - struct. A struct containing the LMB filter's
%           approximate MAP estimate for each time-step of the simulation, as
%           well as the objects' trajectories.

%% Initialise variables
simulationLength = length(measurements);
% Struct containing objects' Bernoulli parameters and metadata
objects = model.object;
% Output struct
stateEstimates.labels = cell(simulationLength, 1);
stateEstimates.mu = cell(simulationLength, 1);
stateEstimates.Sigma = cell(simulationLength, 1);
stateEstimates.objects = objects;
%% Run the LMB filter
for t = 1:simulationLength
    %% Prediction
    objects = lmbPredictionStep(objects, model, t);
    %% Measurement update
    measurementUpdatedDistributions = cell(1, model.numberOfSensors);
    for s = 1:model.numberOfSensors
        if (numel(measurements{s, t}))
            % Populate the association matrices required by the data association algorithms
            [associationMatrices, posteriorParameters] = generateLmbSensorAssociationMatrices(objects, measurements{s, t}, model, s);
            if (strcmp(model.dataAssociationMethod, 'LBP'))
                % Data association by way of loopy belief propagation
                [r, W] = loopyBeliefPropagation(associationMatrices, model.lbpConvergenceTolerance, model.maximumNumberOfLbpIterations);
            elseif(strcmp(model.dataAssociationMethod, 'Gibbs'))
                % Data association by way of Gibbs sampling
                [r, W] = lmbGibbsSampling(associationMatrices, model.numberOfSamples);
            else
                % Data association by way of Murty's algorithm
                [r, W] = lmbMurtysAlgorithm(associationMatrices, model.numberOfAssignments);
            end
            % Compute measurement-updated spatial distributions
            measurementUpdatedDistributions{s} = computePosteriorLmbSpatialDistributions(objects, r, W, posteriorParameters, model);
        else
            % No measurements collected
            measurementUpdatedDistributions{s} = objects;
            for i = 1:numel(objects)
                measurementUpdatedDistributions{s}(i).r = (measurementUpdatedDistributions{s}(i).r * (1 - model.detectionProbability(s))) / (1 - measurementUpdatedDistributions{s}(i).r * model.detectionProbability(s));
            end
        end
    end
    %% Track merging
    if (strcmp(model.lmbParallelUpdateMode, 'AA'))
        objects = aaLmbTrackMerging(measurementUpdatedDistributions, model);
    elseif (strcmp(model.lmbParallelUpdateMode, 'GA'))
        objects = gaLmbTrackMerging(measurementUpdatedDistributions, model);
    else
        objects = puLmbTrackMerging(measurementUpdatedDistributions, objects, model);
    end
    %% Gate tracks
    % Determine which objects have high existence probabilities
    objectsLikelyToExist = [objects.r] > model.existenceThreshold;
    % Objects with low existence probabilities and long trajectories are worth exporting
    discardedObjects = objects(~objectsLikelyToExist & ([objects.trajectoryLength] > model.minimumTrajectoryLength));
    stateEstimates.objects(end+1:end+numel(discardedObjects)) =  discardedObjects;
    % Keep objects with high existence probabilities
    objects = objects(objectsLikelyToExist);
    %% MAP cardinality extraction
    % Determine approximate MAP estimate of the posterior LMB
    [nMap, mapIndices] = lmbMapCardinalityEstimate([objects.r]);
    % Extract RFS state estimate
    stateEstimates.labels{t} = zeros(2, nMap);
    stateEstimates.mu{t} = cell(1, nMap);
    stateEstimates.Sigma{t} = cell(1, nMap);
    for i = 1:nMap
        j = mapIndices(i);
        % Gaussians in the posterior GM are sorted according to weight
        stateEstimates.labels{t}(:, i) = [objects(j).birthTime; objects(j).birthLocation];
        stateEstimates.mu{t}{i} = objects(j).mu{1};
        stateEstimates.Sigma{t}{i} = objects(j).Sigma{1};
    end
    %% Update each object's trajectory
    for i = 1:numel(objects)
        j = objects(i).trajectoryLength;
        objects(i).trajectoryLength = j + 1;
        objects(i).trajectory(:, j+1) = objects(i).mu{1};
        objects(i).timestamps(j+1) = t;
    end 
end
%% Get any long trajectories that weren't extracted
discardedObjects = objects(([objects.trajectoryLength] > model.minimumTrajectoryLength));
numberOfDiscardedObjects = numel(discardedObjects);
stateEstimates.objects(end+1:end+numberOfDiscardedObjects) =  discardedObjects;
end