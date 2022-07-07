% MULTISENSORACCURACYTRIAL -- Compare our multi-sensors filters for a
% basic scenario.
%

%% Admin
clc; close all;
setPath;
%% Open existing file
filename = 'multisensorAccuracyTrials2.mat';
fileExists = isfile(fullfile(cd, filename));
%% Load old simulation, if it exists
if (fileExists)
    load(filename);
else
    % Number of trials
    numberOfLmbTrials = 1000;
    numberOfLmbmTrials = 100;
    numberOfSensors = 3;
    clutterRates = [5 5 5];
    detectionProbabilities = [0.67 0.7 0.73];
    q = [4 3 2];
    model = generateMultisensorModel(numberOfSensors, clutterRates, detectionProbabilities, q, 'PU', 'LBP', 'Fixed');
    simulationLength = 100;
    lmbTrialIndex = 0;
    lmbmTrialIndex = 0;
    % Data assocation methods
    multisensorLmbUpdateMethod = {'IC', 'PU', 'GA', 'AA'};
    numberOfLmbUpdateMethods = numel(multisensorLmbUpdateMethod);
    % E-OSPA error metrics
    eOspaLmb = repmat({zeros(numberOfLmbTrials, simulationLength)}, 1, numberOfLmbUpdateMethods);
    eOspaLmbm = zeros(numberOfLmbmTrials, simulationLength);
    % H-OSPA error metrics
    hOspaLmb = repmat({zeros(numberOfLmbTrials, simulationLength)}, 1, numberOfLmbUpdateMethods);
    hOspaLmbm = zeros(numberOfLmbmTrials, simulationLength);    
    % Cardinalities
    lmbCardinality = repmat({zeros(numberOfLmbTrials, simulationLength)}, 1, numberOfLmbUpdateMethods);
    lmbmCardinality = zeros(numberOfLmbmTrials, simulationLength);
    % Plot colours
    lmbPlotColours = {'c', 'r', 'g', 'b'};
    lmbmPlotColours = 'm';
    % Save all current variables
    save(filename);
end
%% Run the Lmb experiments
startingPoint = lmbTrialIndex;
for t = (startingPoint+1):numberOfLmbTrials
    fprintf('Lmb Trial %d \n', t);
    %% Generate new measurements
    [groundTruth, measurements, groundTruthRfs] = generateMultisensorGroundTruth(model);
    %% LMB filters
    for i = 1:numberOfLmbUpdateMethods
        fprintf([multisensorLmbUpdateMethod{i} '\n']);
        if (strcmp(multisensorLmbUpdateMethod{i}, 'IC'))
            stateEstimates = runIcLmbFilter(model, measurements);
        else
            model.lmbParallelUpdateMode = multisensorLmbUpdateMethod{i};
            stateEstimates = runParallelUpdateLmbFilter(model, measurements);
        end
        [eOspaLmb{i}(t, :), hOspaLmb{i}(t, :), lmbCardinality{i}(t, :)] = computeSimulationOspa(model, groundTruthRfs, stateEstimates);
    end
    %% Save the current results (loadshedding)
    lmbTrialIndex = t;
    save(filename, 'eOspaLmb', 'hOspaLmb', 'lmbCardinality', 'lmbTrialIndex', '-append')
end
%% Run the LMBM experiments
startingPoint = lmbmTrialIndex;
for t = (startingPoint+1):numberOfLmbmTrials
    fprintf('Lmbm Trial %d \n', t);
    %% Generate new measurements
    [groundTruth, measurements, groundTruthRfs] = generateMultisensorGroundTruth(model);
    %% LMBM filters
    stateEstimates = runMultisensorLmbmFilter(model, measurements);
    [eOspaLmbm(t, :), hOspaLmbm(t, :), lmbmCardinality(t, :)] = computeSimulationOspa(model, groundTruthRfs, stateEstimates);
    %% Save results
    lmbmTrialIndex = t;
    save(filename, 'eOspaLmbm', 'hOspaLmbm', 'lmbmCardinality' , 'lmbmTrialIndex', '-append')
end
%% Plot error metrics
t = model.T * (0:simulationLength-1);
% Plot
figure();
%% E-OSPA
subplot(311); box on; hold on; grid on;
% Lmb 
for i = 1:numberOfLmbUpdateMethods
    meanEOspaLmb = mean(eOspaLmb{i}, 1);
    plot(t, meanEOspaLmb, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Lmbm
meanEOspaLmbm = mean(eOspaLmbm, 1);
plot(t, meanEOspaLmbm, 'LineWidth', 1.5, 'Color', lmbmPlotColours);
% Plot control
xlim([t(1) t(end)]);
ylim([0 model.ospaParameters.eC + 0.1])
xlabel('t (s)');
ylabel('E-OSPA');
%% H-OSPA
subplot(312); box on; hold on; grid on;
% Lmb 
for i = 1:numberOfLmbUpdateMethods
    meanHOspaLmb = mean(hOspaLmb{i}, 1);
    plot(t, meanHOspaLmb, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Lmbm
meanHOspaLmbm = mean(hOspaLmbm, 1);
plot(t, meanHOspaLmbm, 'LineWidth', 1.5, 'Color', lmbmPlotColours);
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
for i = 1:numberOfLmbUpdateMethods
    meanLmbCardinality = round(mean(lmbCardinality{i}, 1));
    maxCardinality = max(max(meanLmbCardinality), maxCardinality);
    plot(t, meanLmbCardinality, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Lmbm
meanLmbmCardinality = round(mean(lmbmCardinality, 1));
maxCardinality = max(max(meanLmbmCardinality), maxCardinality);
plot(t, meanLmbmCardinality, 'LineWidth', 1.5, 'Color', lmbmPlotColours);
% Plot control
yMaxLimit = max(max(groundTruthRfs.cardinality), maxCardinality);
xlim([t(1) t(end)]);
ylim([0 yMaxLimit + 1])
xlabel('t (s)');
ylabel('Cardinality');
% matlab2tikz([filename(1:end-4) 'Ospa.tex']);