% MULTISENSORDETECTIONPROBABILITYTRIAL -- Compare our multi-sensors filters for
% varying detection probabilities
%

%% Admin
clc; close all;
setPath;
%% Open existing file
filename = 'multisensorDetectionProbabilityTrials.mat';
fileExists = isfile(fullfile(cd, filename));
%% Load old simulation, if it exists
if (fileExists)
    load(filename);
else
    % Number of trials
    numberOfTrials = 100;
    numberOfSensors = 3;
    detectionProbabilities = [0.5 0.6 0.7 0.8 0.9 0.999];
    numberOfExperimentsPerTrial = numel(detectionProbabilities);
    clutterRates = [5 5 5];
    q = [4 3 2];
    model = generateMultisensorModel(numberOfSensors, clutterRates, 0.5 * ones(1, numberOfSensors), q, 'PU', 'LBP', 'Fixed');
    simulationLength = 100;
    trialIndex = 0;
    % Data assocation methods
    multisensorLmbUpdateMethod = {'IC', 'PU', 'GA', 'AA'};
    numberOfLmbUpdateMethods = numel(multisensorLmbUpdateMethod);
    % E-OSPA error metrics
    eOspaLmb = repmat({zeros(numberOfTrials, numberOfExperimentsPerTrial)}, 1, numberOfLmbUpdateMethods);
    % H-OSPA error metrics
    hOspaLmb = repmat({zeros(numberOfTrials, numberOfExperimentsPerTrial)}, 1, numberOfLmbUpdateMethods);
    % Cardinalities
    lmbCardinality = repmat({zeros(numberOfTrials, numberOfExperimentsPerTrial)}, 1, numberOfLmbUpdateMethods);
    % Plot colours
    lmbPlotColours = {'c', 'r', 'g', 'b'};
    % Save all current variables
    save(filename);
end
%% Run the experiments
startingPoint = trialIndex;
for t = (startingPoint+1):numberOfTrials
    %% LMB filters
    for i = 1:numberOfExperimentsPerTrial
       fprintf('Trial %d, Detection probability %d \n', t, detectionProbabilities(i));
        %% Generate new model and measurements
        model = generateMultisensorModel(numberOfSensors, clutterRates, detectionProbabilities(i) * ones(1, numberOfSensors), q, 'PU', 'LBP', 'Fixed');
        [groundTruth, measurements, groundTruthRfs] = generateMultisensorGroundTruth(model);
        %% Lmb filters
        for j = 1:numberOfLmbUpdateMethods
            fprintf([multisensorLmbUpdateMethod{j} '\n']);
            if (strcmp(multisensorLmbUpdateMethod{j}, 'IC'))
                stateEstimates = runIcLmbFilter(model, measurements);
            else
                model.lmbParallelUpdateMode = multisensorLmbUpdateMethod{j};
                stateEstimates = runParallelUpdateLmbFilter(model, measurements);
            end
            [eOspa, hOspa] = computeSimulationOspa(model, groundTruthRfs, stateEstimates);
            eOspaLmb{j}(t, i) = mean(eOspa, 2);
            hOspaLmb{j}(t, i) = mean(hOspa, 2);
        end
    end
    %% Save the current results (loadshedding)
    trialIndex = t;
    save(filename, 'eOspaLmb', 'hOspaLmb', 'trialIndex', '-append')
end
%% Plot error metrics
figure();
% E-OSPA
subplot(211); hold on; grid on; box on;
for i = 1:numberOfLmbUpdateMethods
    meanEOspaError = mean(eOspaLmb{i}, 1);
    stdEOspaError = std(eOspaLmb{i}, 0, 1);
    errorbar(detectionProbabilities, meanEOspaError, stdEOspaError, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Plot control
xlim([min(detectionProbabilities) max(detectionProbabilities)]);
ylim([0 model.ospaParameters.eC + 0.5])
xlabel('Detection probability');
ylabel('E-OSPA');
% H-OSPA
subplot(212); hold on; grid on; box on;
for i = 1:numberOfLmbUpdateMethods
    meanHOspaError = mean(hOspaLmb{i}, 1);
    stdHOspaError = std(hOspaLmb{i}, 0, 1);
    errorbar(detectionProbabilities, meanHOspaError, stdHOspaError, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Plot control
%ax = gca;
%ax.XTick = numberOfClutterReturns;
xlim([min(detectionProbabilities) max(detectionProbabilities)]);
ylim([0 model.ospaParameters.hC + 0.1])
xlabel('Detection probability');
ylabel('H-OSPA');
% matlab2tikz([filename(1:end-4) 'OSPA.tex']);