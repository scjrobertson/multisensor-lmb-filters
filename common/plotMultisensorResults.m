function plotMultisensorResults(model, measurements, groundTruth, stateEstimates, groundTruthRfs)
% PLOTMULTISENSORRESULTS -- Plots the ground truth and filter state estimates
%   plotMultisensorResults(model, measurements, groundTruth, stateEstimates, groundTruthRfs)
%
%   Plots the ground truth and filter state estimates. Plots include
%   instantaneous RFS state estimates, as well as coherent trajectories.
%   Coherent trajectories have a state estimate for each time-step of an
%   object's trajectory. Instantaneous RFS state estimates may have
%   "gaps" in an object's trajectory, as the state extraction algorithm may
%   not extract a state estimate for an object at a particular time-step.
%
%   See also generateModel, generateGroundtruth, runLmbFilter
%
%   Inputs
%       model - struct. A struct with the fields declared in generateModel.
%       measurements - cell array. An array containing the measurements for
%           each time-step of the simulation.
%       groundTruth - cell array. An array of each object's groundtruth
%           trajectory. This produced by generateGroundtruth.
%       stateEstimates - struct. A struct containing a filter's instantaneous
%           state estimates, and the object trajectories.
%       groundTruthRfs - struct. Groundtruth RFS and optimal Kalman filter
%           output in RFS form.

%% Line specifiers
LINE_COLOURS = {'r', 'g', 'c', 'm', 'y', 'b'};
LINE_MARKERS = {'+', 'o', '*', '.', 'x', 's', 'd'};
%% Unique line specifiers
numberOfColours = size(LINE_COLOURS, 2);
numberOfMarkers = size(LINE_MARKERS, 2);
colourIndices = repmat(1:numberOfColours, [1 numberOfMarkers]);
markerIndices = reshape(repmat(1:numberOfMarkers, [numberOfColours 1]),[1 numberOfColours*numberOfMarkers]);
lineSpecifiers = [colourIndices; markerIndices];
lineSpecifiers = lineSpecifiers(:, randperm(numberOfColours*numberOfMarkers));
numberOfUniqueLines = size(lineSpecifiers, 2);
%% Determine labels' line specifiers
% Determine all unique labels
instantaneousLabels = [stateEstimates.labels{:}];
trajectoryLabels = [stateEstimates.objects.birthTime; stateEstimates.objects.birthLocation];
labels = [instantaneousLabels trajectoryLabels];
distinctLabels = unique(labels', 'rows')';
numberOfDistinctLabels = size(distinctLabels, 2);
% Determine their line specifiers
lineSpecifierIndices = repmat(1:numberOfUniqueLines, 1, ceil(numberOfDistinctLabels / numberOfUniqueLines));
lineSpecifierIndices = lineSpecifierIndices(1:numberOfDistinctLabels);
%% Populate time stamps
simulationLength = numel(stateEstimates.mu);
timestamps = cell(1, simulationLength);
cardinalityEstimate = size(1, simulationLength);
for i = 1:simulationLength
    cardinalityEstimate(i) = numel(stateEstimates.mu{i});
    timestamps{i} = model.T * (i-1) * ones(1,  cardinalityEstimate(i));
end
packedTimestamps = [timestamps{:}];
%% Stack measurements
measurementTimestamps = repmat({zeros(1, 0)}, 1, model.numberOfSensors);
z = repmat({zeros(model.zDimension, 0)}, 1, model.numberOfSensors);
for s = 1:model.numberOfSensors
    for i = 1:simulationLength
        repeatedTimestamps = model.T * (i-1) * ones(1, numel(measurements{s, i}));
        measurementTimestamps{s} = horzcat(measurementTimestamps{s}, repeatedTimestamps);
        y = measurements{s, i};
        if (numel(y) > 0)
            z{s} = horzcat(z{s}, [y{:}]);
        end
    end
end
%% RFS estimate plot
figure(); hold on; grid on; box on;
% Ground truth
% Plot ground truth
for i = 1:length(groundTruth)
    projectedGroundTruth = groundTruth{i}(2:end, :);
    birthPoint = line(projectedGroundTruth(1, 1), projectedGroundTruth(2, 1), 'LineStyle','none','Marker','o', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
    truthLine = line(projectedGroundTruth(1, :), projectedGroundTruth(2, :), 'LineStyle','-','Marker','none','LineWidth', 3,'Color',0*ones(1,3));
    deathPoint = line(projectedGroundTruth(1, end), projectedGroundTruth(2, end), 'LineStyle','none','Marker','^', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
end
% Plot measurements
for s = 1:model.numberOfSensors
    plot(z{s}(1, :), z{s}(2, :), 'LineStyle', 'none', 'Marker', LINE_MARKERS{s}, 'MarkerSize', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Sensor %d measurements', s));
end
% State estimates
mu = [stateEstimates.mu{:}];
rfsTrajectories = struct('lineIndex', 0, 'displayName', [], 'trajectory', [], 'timestamps', []);
rfsTrajectories = repmat(rfsTrajectories, 1, numberOfDistinctLabels);
numberOfRfsTrajectories = 0;
for i = 1:numberOfDistinctLabels
    % Pack trajectory
    labelIndicator = logical(prod(instantaneousLabels == distinctLabels(:, i), 1));
    estimates = mu(labelIndicator);
    trajectory = [estimates{:}];
    trajectoryTimestamps = packedTimestamps(labelIndicator);
    if (~isempty(trajectory))
        % Plot trajectory
        j = lineSpecifierIndices(i);
        displayName = sprintf('Label = (%d, %d)\n', distinctLabels(1, i), distinctLabels(2, i));
        plot(trajectory(1, :), trajectory(2, :), 'LineStyle','none', 'Marker', LINE_MARKERS{lineSpecifiers(2, j)}, 'Markersize', 8, 'Color', LINE_COLOURS{lineSpecifiers(1, j)}, 'DisplayName', displayName);
        % Save result
        numberOfRfsTrajectories = numberOfRfsTrajectories + 1;
        rfsTrajectories(numberOfRfsTrajectories).lineIndex = j;
        rfsTrajectories(numberOfRfsTrajectories).displayName = displayName;
        rfsTrajectories(numberOfRfsTrajectories).trajectory = trajectory;
        rfsTrajectories(numberOfRfsTrajectories).timestamps = trajectoryTimestamps;
    end
end
rfsTrajectories = rfsTrajectories(1:numberOfRfsTrajectories);
% Plot control
xlabel('x (m)');
ylabel('y (m)');
title('Instantaneous state estimates');
xlim(model.observationSpaceLimits(1, :));
ylim(model.observationSpaceLimits(2, :));
% matlab2tikz([model.lmbParallelUpdateMode 'LmbXyPlot.tex']);
%% Plot results versus time
figure; hold on; 
%% x-happenings vs. time
subplot(211); box on; hold on; grid on;
% Plot ground truth
for i = 1:length(groundTruth)
    t = groundTruth{i}(1, :);
    projectedGroundTruth = groundTruth{i}(2:end, :);
    birthPoint = line(t(1, 1), projectedGroundTruth(1, 1), 'LineStyle','none','Marker','o', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
    truthLine = line(t(1, :), projectedGroundTruth(1, :), 'LineStyle','-','Marker','none','LineWidth', 3,'Color',0*ones(1,3));
    deathPoint = line(t(1, end), projectedGroundTruth(1, end), 'LineStyle','none','Marker','^', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
end
% Plot measurements
for s = 1:model.numberOfSensors
    plot(measurementTimestamps{s}, z{s}(1, :), 'LineStyle', 'none', 'Marker', LINE_MARKERS{s}, 'MarkerSize', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Sensor %d measurements', s));
end
% Plot simulation
for i = 1:numberOfRfsTrajectories
    % Get saved data
    j = rfsTrajectories(i).lineIndex;
    displayName = rfsTrajectories(i).displayName;
    trajectory = rfsTrajectories(i).trajectory(1, :);
    timestamps = rfsTrajectories(i).timestamps;
    % Plot
    plot(timestamps, trajectory, 'LineStyle','none', 'Marker', LINE_MARKERS{lineSpecifiers(2, j)}, 'Markersize', 8, 'Color', LINE_COLOURS{lineSpecifiers(1, j)}, 'DisplayName', displayName);
end
% Plot control
xlim([0 simulationLength]); ylim(model.observationSpaceLimits(1, :));
xlabel('t (s)');
ylabel('x (m)');
title('Instantaneous state estimates');
%% y-happenings vs. time
subplot(212); box on; hold on; grid on;
% Plot ground truth
for i = 1:length(groundTruth)
    t = groundTruth{i}(1, :);
    projectedGroundTruth = groundTruth{i}(2:end, :);
    birthPoint = line(t(1, 1), projectedGroundTruth(2, 1), 'LineStyle','none','Marker','o', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
    truthLine = line(t(1, :), projectedGroundTruth(2, :), 'LineStyle','-','Marker','none','LineWidth', 3,'Color',0*ones(1,3));
    deathPoint = line(t(1, end), projectedGroundTruth(2, end), 'LineStyle','none','Marker','^', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
end
% Plot measurements
for s = 1:model.numberOfSensors
    plot(measurementTimestamps{s}, z{s}(2, :), 'LineStyle', 'none', 'Marker', LINE_MARKERS{s}, 'MarkerSize', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Sensor %d measurements', s));
end
% Plot simulation
for i = 1:numberOfRfsTrajectories
    % Get saved data
    j = rfsTrajectories(i).lineIndex;
    displayName = rfsTrajectories(i).displayName;
    trajectory = rfsTrajectories(i).trajectory(2, :);
    timestamps = rfsTrajectories(i).timestamps;
    % Plot
    plot(timestamps, trajectory, 'LineStyle','none', 'Marker', LINE_MARKERS{lineSpecifiers(2, j)}, 'Markersize', 8, 'Color', LINE_COLOURS{lineSpecifiers(1, j)}, 'DisplayName', displayName);
end
% Plot simulation
xlim([0 simulationLength]); ylim(model.observationSpaceLimits(2, :));
xlabel('t (s)');
ylabel('y (m)');
% matlab2tikz([model.lmbParallelUpdateMode 'LmbXyVTime.tex']);
%% Built up trajectory plot
figure(); hold on; grid on; box on;
% Ground truth
% Plot ground truth
for i = 1:length(groundTruth)
    projectedGroundTruth = groundTruth{i}(2:end, :);
    birthPoint = line(projectedGroundTruth(1, 1), projectedGroundTruth(2, 1), 'LineStyle','none','Marker','o', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
    truthLine = line(projectedGroundTruth(1, :), projectedGroundTruth(2, :), 'LineStyle','-','Marker','none','LineWidth', 3,'Color',0*ones(1,3));
    deathPoint = line(projectedGroundTruth(1, end), projectedGroundTruth(2, end), 'LineStyle','none','Marker','^', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
end
% Plot measurements
for s = 1:model.numberOfSensors
    plot(z{s}(1, :), z{s}(2, :), 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Sensor %d measurements', s));
end
% State estimates
for i = 1:numel(stateEstimates.objects)
    object = stateEstimates.objects(i);
    objectLabel =  [object.birthTime; object.birthLocation];
    trajectoryLength = object.trajectoryLength;
    % Plot trajectory
    j = lineSpecifierIndices(logical(prod(distinctLabels == objectLabel, 1)));
    displayName = sprintf('Label = (%d, %d)\n', objectLabel(1, 1), objectLabel(2, 1));
    plot(stateEstimates.objects(i).trajectory(1, 1:trajectoryLength), stateEstimates.objects(i).trajectory(2, 1:trajectoryLength), 'LineStyle','none', 'Marker', LINE_MARKERS{lineSpecifiers(2, j)}, 'Markersize', 8, 'Color', LINE_COLOURS{lineSpecifiers(1, j)}, 'DisplayName', displayName);
end
% Plot control
xlabel('x (m)');
ylabel('y (m)');
title('Coherent object trajectories');
xlim(model.observationSpaceLimits(1, :));
ylim(model.observationSpaceLimits(2, :));
%% Plot results versus time
figure; hold on; 
%% x-happenings vs. time
subplot(211); box on; hold on; grid on;
% Plot ground truth
for i = 1:length(groundTruth)
    t = groundTruth{i}(1, :);
    projectedGroundTruth = groundTruth{i}(2:end, :);
    birthPoint = line(t(1, 1), projectedGroundTruth(1, 1), 'LineStyle','none','Marker','o', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
    truthLine = line(t(1, :), projectedGroundTruth(1, :), 'LineStyle','-','Marker','none','LineWidth', 3,'Color',0*ones(1,3));
    deathPoint = line(t(1, end), projectedGroundTruth(1, end), 'LineStyle','none','Marker','^', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
end
% Plot measurements
for s = 1:model.numberOfSensors
    plot(measurementTimestamps{s}, z{s}(1, :), 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Sensor %d measurements', s));
end
% State estimates
for i = 1:numel(stateEstimates.objects)
    object = stateEstimates.objects(i);
    objectLabel =  [object.birthTime; object.birthLocation];
    trajectoryLength = object.trajectoryLength;
    % Plot trajectory
    j = lineSpecifierIndices(logical(prod(distinctLabels == objectLabel, 1)));
    displayName = sprintf('Label = (%d, %d)\n', objectLabel(1, 1), objectLabel(2, 1));
    plot(stateEstimates.objects(i).timestamps(1, 1:trajectoryLength), stateEstimates.objects(i).trajectory(1, 1:trajectoryLength), 'LineStyle','none', 'Marker', LINE_MARKERS{lineSpecifiers(2, j)}, 'Markersize', 8, 'Color', LINE_COLOURS{lineSpecifiers(1, j)}, 'DisplayName', displayName);
end
% Plot control
xlim([0 simulationLength]); ylim(model.observationSpaceLimits(1, :));
xlabel('t (s)');
ylabel('x (m)');
title('Coherent state estimates');
%% y-happenings vs. time
subplot(212); box on; hold on; grid on;
% Plot ground truth
for i = 1:length(groundTruth)
    t = groundTruth{i}(1, :);
    projectedGroundTruth = groundTruth{i}(2:end, :);
    birthPoint = line(t(1, 1), projectedGroundTruth(2, 1), 'LineStyle','none','Marker','o', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
    truthLine = line(t(1, :), projectedGroundTruth(2, :), 'LineStyle','-','Marker','none','LineWidth', 3,'Color',0*ones(1,3));
    deathPoint = line(t(1, end), projectedGroundTruth(2, end), 'LineStyle','none','Marker','^', 'Markersize', 12,'LineWidth',1,'Color',0*ones(1,3));
end
% Plot measurements
for s = 1:model.numberOfSensors
    plot(measurementTimestamps{s}, z{s}(2, :), 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 2, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Sensor %d measurements', s));
end
% State estimates
for i = 1:numel(stateEstimates.objects)
    object = stateEstimates.objects(i);
    objectLabel =  [object.birthTime; object.birthLocation];
    trajectoryLength = object.trajectoryLength;
    % Plot trajectory
    j = lineSpecifierIndices(logical(prod(distinctLabels == objectLabel, 1)));
    displayName = sprintf('Label = (%d, %d)\n', objectLabel(1, 1), objectLabel(2, 1));
    plot(stateEstimates.objects(i).timestamps(1, 1:trajectoryLength), stateEstimates.objects(i).trajectory(2, 1:trajectoryLength), 'LineStyle','none', 'Marker', LINE_MARKERS{lineSpecifiers(2, j)}, 'Markersize', 8, 'Color', LINE_COLOURS{lineSpecifiers(1, j)}, 'DisplayName', displayName);
end
% Plot simulation
xlim([0 simulationLength]); ylim(model.observationSpaceLimits(2, :));
xlabel('t (s)');
ylabel('y (m)');
%% Plot error metrics
% Determine errors
t = model.T * (0:simulationLength-1);
[eOspa, hOspa] = computeSimulationOspa(model, groundTruthRfs, stateEstimates);
% Plot
figure();
% E-OSPA
subplot(311); box on; hold on; grid on;
plot(t, eOspa, 'Color', 'r', 'LineWidth', 1.5);
xlim([t(1) t(end)]);
ylim([0 model.ospaParameters.eC + 0.1])
xlabel('t (s)');
ylabel('E-OSPA');
% H-OSPA
subplot(312); box on; hold on; grid on;
plot(t, hOspa, 'Color', 'r', 'LineWidth', 1.5);
xlim([t(1) t(end)]);
ylim([0 model.ospaParameters.hC + 0.1])
xlabel('t (s)');
ylabel('H-OSPA');
% Cardinality
subplot(313); box on; hold on; grid on;
stairs(t, groundTruthRfs.cardinality, 'Color', 'b', 'LineWidth', 2);
stairs(t, cardinalityEstimate, 'Color', 'r', 'LineWidth', 1.5);
yMaxLimit = max([max(groundTruthRfs.cardinality), max(cardinalityEstimate)]);
xlim([t(1) t(end)]);
ylim([0 yMaxLimit + 1])
xlabel('t (s)');
ylabel('Cardinality');
end