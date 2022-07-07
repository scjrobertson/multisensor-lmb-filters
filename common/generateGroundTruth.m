function [groundTruth, measurements, groundTruthRfs] = generateGroundTruth(model, varargin)
% GENERATEGROUNDTRUTH -- Generates the measurements and groundtruth for a simple hard-coded example.
%   [groundTruth, measurements, groundTruthRfs] = generateGroundTruth(model)
%
%   Generates the objects' groundtruths for a simple scenario, and also their measurements.
%
%   See also generateModel, plotResults
%
%   Inputs
%       model - struct. A struct with the fields declared in generateModel.
%       numberOfObjects - integer. The number of objects to be simulated
%           for a 'Random' scenario.
%
%   Output
%       groundTruth - cell array. An array of each object's groundtruth
%           trajectory.
%       measurements - cell array. An array containing the measurements for
%           each time-step of the simulation.
%       groundTruthRfs - struct. Groundtruth RFS and optimal Kalman filter
%           output in RFS form.


%% Simple, hard-coded scenario
if (strcmp(model.scenarioType, 'Fixed'))
    numberOfObjects = 10;
    % Object birth times
    simulationLength = 100;
    objectBirthTimes = [1, 1, 20, 20, 40, 40, 60, 60, 60, 60];
    objectDeathTimes = [70, 70, 80, 80, 90, 90, 100, 100, 100, 100];
    % Object birth states
    birthLocationIndex = [1 2 3 4 1 4 1 2 3 4];
    priorLocations = [-80.0 -20.0 0.75 1.5;
        -20.0 80.0 -1.0 -2.0;
        0.0 0.0 -0.5 -1.0;
        40.0 -60.0 -0.25 -0.5;
        -80.0 -20.0 1.0 1.0;
        40.0 -60.0 -1.0 2.0;
        -80.0 -20.0 1.0 -0.5;
        -20.0 80.0 1.0 -1.0;
        0.0 0.0 1.0 -1.0;
        40.0 -60.0 -1.0 0.5]';
elseif (strcmp(model.scenarioType, 'Random'))
    % Number of objects
    if (nargin == 2)
        numberOfObjects = varargin{1};
    else
        error('You must specify the number of objects for a Random scenario');
    end
     % Object birth times
    simulationLength = 100;
    objectBirthTimes = ones(1, numberOfObjects); %randi([1 simulationLength], 1, numberOfObjects);
    objectDeathTimes = simulationLength * ones(1, numberOfObjects); %randi([1 simulationLength], 1, numberOfObjects) + objectBirthTimes + model.minimumTrajectoryLength;
    objectDeathTimes(objectDeathTimes > simulationLength) = simulationLength;
    % Object birth states
    birthLocationIndex = 1:model.numberOfBirthLocations; %randi([1 model.numberOfBirthLocations],1 , numberOfObjects);
    priorLocations = [model.muB{birthLocationIndex}];
    priorLocations(3:4, :) = 3 * randn(numberOfObjects, 2)';
end
%% Allocate output
measurements = repmat({}, simulationLength, 1);
groundTruth = repmat({}, numberOfObjects, 1);
groundTruthRfs.x = repmat({{}}, 1, simulationLength);
groundTruthRfs.mu = groundTruthRfs.x;
groundTruthRfs.Sigma = groundTruthRfs.x;
groundTruthRfs.cardinality = zeros(1, simulationLength);
%% Add in clutter measurements
for i = 1:simulationLength
    numberOfClutterMeasurements = poissrnd(model.clutterRate);
    measurements{i} = cell(numberOfClutterMeasurements, 1);
    for j = 1:numberOfClutterMeasurements
        measurements{i}{j} = model.observationSpaceLimits(:, 1) + 2 * model.observationSpaceLimits(:, 2) .* rand(model.zDimension, 1);
    end
end
%% Add in each object's measurements
QChol = chol(model.Q, 'lower');
for i = 1:numberOfObjects
    % Initialise the object's trajectory
    trajectoryLength = objectDeathTimes(i) - objectBirthTimes(i) + 1;
    t = objectBirthTimes(i);
    x =  priorLocations(:, i);
    mu = model.muB{birthLocationIndex(i)};
    Sigma = model.SigmaB{birthLocationIndex(i)};
    groundTruth{i} = repmat([t; x], 1, trajectoryLength);
    % Simulate the object's trajectory
    for j = 1:trajectoryLength
        % Predict the object's state
        if (j > 1)
            % Point estimate
            x = model.A * x + model.u;
            t = t + 1;
            groundTruth{i}(:, j) = [t; x];
            % Kalman filter prediction
            mu = model.A * mu + model.u;
            Sigma = model.A * Sigma * model.A' + model.R;
        end
        % Determine if object missed detection
        if (rand < model.detectionProbability)
            z = model.C * x + QChol * randn(1, model.zDimension)';
            measurements{t}{end+1} = z;
            % Kalman filter measurement update
            K = Sigma * model.C'/(model.C * Sigma * model.C' + model.Q);
            mu = mu + K *(z - model.C * mu);
            Sigma = (eye(model.xDimension) - K * model.C) * Sigma;
        end
        % Add to RFS
        groundTruthRfs.x{t}{end+1} = x;  
        groundTruthRfs.mu{t}{end+1} = mu;
        groundTruthRfs.Sigma{t}{end+1} = Sigma;
        groundTruthRfs.cardinality(t) = groundTruthRfs.cardinality(t) + 1;
    end
end
end