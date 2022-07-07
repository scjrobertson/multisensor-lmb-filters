% SINGLESENSORACCURACYTRIAL -- Compare our single-sensors filters for a
% basic scenario.
%

%% Admin
clc; close all;
setPath;
%% Open existing file
filename = 'blimp.mat';
fileExists = isfile(fullfile(cd, filename));
%% Load old simulation, if it exists
if (fileExists)
    load(filename);
else
    % Number of trials
    numberOfTrials = 100;
    model = generateModel(2, 0.95, 'LBP'); % rBLbmbm = 0.06
    simulationLength = 100;
    trialIndex = 0;
    % Data assocation methods
    lmbDataAssociationMethods = {'LBP', 'Gibbs', 'Murty'};
    numberOfLmbAssociationMethods = numel(lmbDataAssociationMethods);
    lmbmDataAssociationMethods = {'Gibbs', 'Murty'};
    numberOfLmbmAssociationMethods = numel(lmbmDataAssociationMethods);
    % E-OSPA error metrics
    eOspaLmb = repmat({zeros(numberOfTrials, simulationLength)}, 1, numberOfLmbAssociationMethods);
    eOspaLmbm = repmat({zeros(numberOfTrials, simulationLength)}, 1, numberOfLmbmAssociationMethods);
    % H-OSPA error metrics
    hOspaLmb = repmat({zeros(numberOfTrials, simulationLength)}, 1, numberOfLmbAssociationMethods);
    hOspaLmbm = repmat({zeros(numberOfTrials, simulationLength)}, 1, numberOfLmbmAssociationMethods);    
    % Cardinalities
    lmbCardinality = repmat({zeros(numberOfTrials, simulationLength)}, 1, numberOfLmbAssociationMethods);
    lmbmCardinality = repmat({zeros(numberOfTrials, simulationLength)}, 1, numberOfLmbmAssociationMethods);
    % Plot colours
    lmbPlotColours = {'r', 'g', 'b'};
    lmbmPlotColours = {'c', 'm'};
    % Save all current variables
    save(filename);
end
%% Run the experiments
startingPoint = trialIndex;
for t = (startingPoint+1):numberOfTrials
    fprintf('Trial %d \n', t);
    %% Generate new measurements
    [groundTruth, measurements, groundTruthRfs] = generateGroundTruth(model);
    %% LMB filters
    for i = 1:numberOfLmbAssociationMethods
        fprintf(['LMB: ' lmbDataAssociationMethods{i} '\n']);
        model.dataAssociationMethod = lmbDataAssociationMethods{i};
        stateEstimates = runLmbFilter(model, measurements);
        [eOspaLmb{i}(t, :), hOspaLmb{i}(t, :), lmbCardinality{i}(t, :)] = computeSimulationOspa(model, groundTruthRfs, stateEstimates);
    end
    %% LMBM filters
    for i = 1:numberOfLmbmAssociationMethods
        fprintf(['LMBM: ' lmbmDataAssociationMethods{i} '\n'])
        model.dataAssociationMethod = lmbmDataAssociationMethods{i};
        stateEstimates = runLmbmFilter(model, measurements);
        [eOspaLmbm{i}(t, :), hOspaLmbm{i}(t, :), lmbmCardinality{i}(t, :)] = computeSimulationOspa(model, groundTruthRfs, stateEstimates);
    end
    %% Save the current results (loadshedding)
    trialIndex = t;
    save(filename, 'eOspaLmb', 'eOspaLmbm', 'hOspaLmb', 'hOspaLmbm', 'lmbCardinality', 'lmbmCardinality' , 'trialIndex', '-append')
end
%% Plot error metrics
t = model.T * (0:simulationLength-1);
% Plot
figure();
%% E-OSPA
subplot(311); box on; hold on; grid on;
% Lmb 
for i = 1:numberOfLmbAssociationMethods
    meanEOspaLmb = mean(eOspaLmb{i}, 1);
    plot(t, meanEOspaLmb, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Lmbm
for i = 1:numberOfLmbmAssociationMethods
    meanEOspaLmbm = mean(eOspaLmbm{i}, 1);
    plot(t, meanEOspaLmbm, 'LineWidth', 1.5, 'Color', lmbmPlotColours{i});
end
% Plot control
xlim([t(1) t(end)]);
ylim([0 model.ospaParameters.eC + 0.1])
xlabel('t (s)');
ylabel('E-OSPA');
%% H-OSPA
subplot(312); box on; hold on; grid on;
% Lmb 
for i = 1:numberOfLmbAssociationMethods
    meanHOspaLmb = mean(hOspaLmb{i}, 1);
    plot(t, meanHOspaLmb, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Lmbm
for i = 1:numberOfLmbmAssociationMethods
    meanHOspaLmbm = mean(hOspaLmbm{i}, 1);
    plot(t, meanHOspaLmbm, 'LineWidth', 1.5, 'Color', lmbmPlotColours{i});
end
% Plot control
xlim([t(1) t(end)]);
ylim([0 model.ospaParameters.hC + 0.1])
xlabel('t (s)');
ylabel('H-OSPA');
%% Cardinality
subplot(313); box on; hold on; grid on;
maxCardinality = -inf;
% Ground truth
stairs(t, groundTruthRfs.cardinality, 'Color', 'k', 'LineWidth', 2);
% Lmb 
for i = 1:numberOfLmbAssociationMethods
    meanLmbCardinality = round(mean(lmbCardinality{i}, 1));
    maxCardinality = max(max(meanLmbCardinality), maxCardinality);
    plot(t, meanLmbCardinality, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Lmbm
for i = 1:numberOfLmbmAssociationMethods
    meanLmbmCardinality = round(mean(lmbmCardinality{i}, 1));
    maxCardinality = max(max(meanLmbmCardinality), maxCardinality);
    plot(t, meanLmbmCardinality, 'LineWidth', 1.5, 'Color', lmbmPlotColours{i});
end
% Plot control
yMaxLimit = max(max(groundTruthRfs.cardinality), maxCardinality);
xlim([t(1) t(end)]);
ylim([0 yMaxLimit + 1])
xlabel('t (s)');
ylabel('Cardinality');
% matlab2tikz([filename(1:end-4) 'Ospa.tex']);