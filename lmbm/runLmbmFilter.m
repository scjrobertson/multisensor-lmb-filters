function stateEstimates = runLmbmFilter(model, measurements)
% RUNLMBMFILTER -- Run the LMBM filter for a given simulated scenario.
%   stateEstimates = runLmbmFilter(model, measurements)
%
%   Determine the objects' state estimates using the LMBM filter.
%
%   See also generateModel, generateGroundTruth
%
%   Inputs
%       model - struct. A struct with the fields declared in generateModel.
%       measurements - cell array. An array containing the measurements for
%           each time-step of the simulation. See also generateGroundTruth.
%
%   Output
%       stateEstimates - struct. A struct containing the LMB filter's
%           approximate MAP estimate for each time-step of the simulation, as
%           well as the objects' trajectories.

%% Initialise variables
simulationLength = length(measurements);
% Struct containing objects' Bernoulli parameters and metadata
hypotheses = model.hypotheses;
objects = model.trajectory;
% Output struct
stateEstimates.labels = cell(simulationLength, 1);
stateEstimates.mu = cell(simulationLength, 1);
stateEstimates.Sigma = cell(simulationLength, 1);
stateEstimates.objects = objects;
%% Run the LMBM filter
for t = 1:simulationLength
    %% Add in new trajectory structs
    [model.birthTrajectory.birthTime] = deal(t);
    objects(end+1:end+model.numberOfBirthLocations) = model.birthTrajectory; 
    %% Preallocate posterior hypotheses
    posteriorHypotheses = repmat(model.hypotheses, 0, 1);
    %% Generate posterior hypotheses for each prior hyptohesis
    for i = 1:numel(hypotheses)
        %% Prediction step
        priorHypothesis = lmbmPredictionStep(hypotheses(i), model, t);
        %% Measurement update
        if (numel(measurements{t}))
            % Generate the Gibbs sampler matrices, and determine the posterior spatial distributions' parameters
            [associationMatrices, posteriorParameters] = generateLmbmAssociationMatrices(priorHypothesis, measurements{t}, model);
            % Generate posterior hypotheses using Gibbs sampling
            if(strcmp(model.dataAssociationMethod, 'Murty'))
                V = murtysAlgorithmWrapper(associationMatrices.C, model.numberOfAssignments);
            else
                V = lmbmGibbsSampling(associationMatrices.P, associationMatrices.C, model.numberOfSamples);
            end
            % Determine each posterior hypothesis' parameters
            newHypotheses = determinePosteriorHypothesisParameters(V, associationMatrices.L, posteriorParameters, priorHypothesis);
            % Add posterior hypotheses to the pile
            posteriorHypotheses(end+1:end+numel(newHypotheses)) = newHypotheses;
        else
            priorHypothesis.r = ((1-model.detectionProbability) * priorHypothesis.r) ./ (1 - model.detectionProbability * priorHypothesis.r);
            posteriorHypotheses(end+1) = priorHypothesis; 
        end
    end
    %% Normalise posterior hypothesis weights and discard unlikely hypothesis
    [hypotheses, objectsLikelyToExist] = lmbmNormalisationAndGating(posteriorHypotheses, model);
    %% Gate trajectories
    % Objects with low existence probabilities and long trajectories are worth exporting
    discardedObjects = objects(~objectsLikelyToExist' & ([objects.trajectoryLength] > model.minimumTrajectoryLength));
    stateEstimates.objects(end+1:end+numel(discardedObjects)) =  discardedObjects;
    % Keep objects with high existence probabilities
    objects = objects(objectsLikelyToExist);
    %% State extraction
    [cardinalityEstimate, extractionIndices] = lmbmStateExtraction(hypotheses, false);
    % Extract RFS state estimate
    stateEstimates.labels{t} = zeros(2, cardinalityEstimate);
    stateEstimates.mu{t} = cell(1, cardinalityEstimate);
    stateEstimates.Sigma{t} = cell(1, cardinalityEstimate);
    for i = 1:cardinalityEstimate
        j = extractionIndices(i);
        % Hypotheses in the posterior LMBM are sorted according to weight
        stateEstimates.labels{t}(:, i) = [hypotheses(1).birthTime(j); hypotheses(1).birthLocation(j)];
        stateEstimates.mu{t}{i} = hypotheses(1).mu{j};
        stateEstimates.Sigma{t}{i} = hypotheses(1).Sigma{j};
    end
    %% Update each object's trajectory
    for i = 1:numel(objects)
        j = objects(i).trajectoryLength;
        objects(i).trajectoryLength = j + 1;
        objects(i).trajectory(:, j+1) = hypotheses(1).mu{i};
        objects(i).timestamps(j+1) = t;
    end 
end
%% Get any long trajectories that weren't extracted
discardedObjects = objects(([objects.trajectoryLength] > model.minimumTrajectoryLength));
numberOfDiscardedObjects = numel(discardedObjects);
stateEstimates.objects(end+1:end+numberOfDiscardedObjects) =  discardedObjects;
end