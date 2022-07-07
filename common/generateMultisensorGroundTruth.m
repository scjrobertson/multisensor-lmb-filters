function [groundTruth, measurements, groundTruthRfs] = generateMultisensorGroundTruth(model, varargin)
% GENERATEMULTISENSORGROUNDTRUTH -- Generates the measurements and groundtruth for a simple hard-coded example.
%   [groundTruth, measurements, groundTruthRfs] = generateGroundTruth(model)
%
%   Generates the objects' groundtruths for a simple scenario, and also their measurements.
%
%   See also generateMultisensorModel, plotMultisensorResults
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
measurements = repmat({}, model.numberOfSensors, simulationLength);
groundTruth = repmat({}, numberOfObjects, 1);
groundTruthRfs.x = repmat({{}}, 1, simulationLength);
groundTruthRfs.mu = groundTruthRfs.x;
groundTruthRfs.Sigma = groundTruthRfs.x;
groundTruthRfs.cardinality = zeros(1, simulationLength);
%% Add in clutter measurements
for i = 1:simulationLength
    for s = 1:model.numberOfSensors
        numberOfClutterMeasurements = poissrnd(model.clutterRate(s));
        measurements{s, i} = cell(numberOfClutterMeasurements, 1);
        for j = 1:numberOfClutterMeasurements
            measurements{s, i}{j} = model.observationSpaceLimits(:, 1) + 2 * model.observationSpaceLimits(:, 2) .* rand(model.zDimension, 1);
        end
    end
end
%% Add in each object's measurements
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
        % Preallocate parameters
        generatedMeasurement = rand(model.numberOfSensors, 1) < model.detectionProbability;
        numberOfAssignments = sum(generatedMeasurement);
        if (numberOfAssignments)
            counter = 0;
            z = zeros(model.zDimension * numberOfAssignments, 1);
            C = zeros(model.zDimension * numberOfAssignments, model.xDimension);
            % Generate measurements
            for s = 1:model.numberOfSensors
                % Determine if object missed detection
                if (generatedMeasurement(s))
                    y = model.C{s} * x + chol(model.Q{s}, 'lower') * randn(1, model.zDimension)';
                    measurements{s, t}{end+1} = y;
                    % Stack measurements and matrices
                    start = model.zDimension * counter + 1;
                    finish = start + model.zDimension - 1;
                    z(start:finish) = y;
                    C(start:finish, :) = model.C{s};
                    counter = counter + 1;
                end
            end
            % Multisensor Kalman filter update
            Q = blkdiag(model.Q{generatedMeasurement});
            K = Sigma * C' / ( C * Sigma * C' + Q);
            mu = mu + K * (z - C * mu);
            Sigma = (eye(model.xDimension) - K * C) * Sigma;
        end
        % Add to RFS
        groundTruthRfs.x{t}{end+1} = x;  
        groundTruthRfs.mu{t}{end+1} = mu;
        groundTruthRfs.Sigma{t}{end+1} = Sigma;
        groundTruthRfs.cardinality(t) = groundTruthRfs.cardinality(t) + 1;
    end
end
end