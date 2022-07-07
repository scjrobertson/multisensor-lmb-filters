% SINGLESENSORCLUTTERTRIAL -- Compare our single-sensors filters for
% varying numbers of clutter returns
%

%% Admin
clc; close all;
setPath;
%% Open existing file
filename = 'highDetectionClutterTrials.mat';
fileExists = isfile(fullfile(cd, filename));
%% Load old simulation, if it exists
if (fileExists)
    load(filename);
else
    % Number of trials
    numberOfTrials = 100;
    numberOfClutterReturns = [10 20 30 40 50 60 70 80 90 100];
    numberOfExperimentsPerTrial = numel(numberOfClutterReturns);
    simulationLength = 100;
    trialIndex = 0;
    % Data assocation methods
    lmbDataAssociationMethods = {'LBP', 'Gibbs', 'Murty'};
    numberOfLmbAssociationMethods = numel(lmbDataAssociationMethods);
    lmbmDataAssociationMethods = {'Gibbs', 'Murty'};
    numberOfLmbmAssociationMethods = numel(lmbmDataAssociationMethods);
    % E-OSPA error metrics
    eOspaLmb = repmat({zeros(numberOfTrials, numberOfExperimentsPerTrial)}, 1, numberOfLmbAssociationMethods);
    eOspaLmbm = repmat({zeros(numberOfTrials, numberOfExperimentsPerTrial)}, 1, numberOfLmbmAssociationMethods);
    % H-OSPA error metrics
    hOspaLmb = repmat({zeros(numberOfTrials, numberOfExperimentsPerTrial)}, 1, numberOfLmbAssociationMethods);
    hOspaLmbm = repmat({zeros(numberOfTrials, numberOfExperimentsPerTrial)}, 1, numberOfLmbmAssociationMethods);    
    % Plot colours
    lmbPlotColours = {'r', 'g', 'b'};
    lmbmPlotColours = {'c', 'm'};
    % Save all current variables
    save(filename);
end
%% Run the experiments
startingPoint = trialIndex;
for t = (startingPoint+1):numberOfTrials
    for i = 1:numberOfExperimentsPerTrial
       fprintf('Trial %d, Number of Clutter returns %d \n', t, numberOfClutterReturns(i));
        %% Generate new model and measurements
        model = generateModel(numberOfClutterReturns(i), 0.95, 'LBP');
        [groundTruth, measurements, groundTruthRfs] = generateGroundTruth(model);
        %% Lmb filters
        for j = 1:numberOfLmbAssociationMethods
            fprintf(['LMB: ' lmbDataAssociationMethods{j} '\n']);
            model.dataAssociationMethod = lmbDataAssociationMethods{j};
            stateEstimates = runLmbFilter(model, measurements);
            [eOspa, hOspa] = computeSimulationOspa(model, groundTruthRfs, stateEstimates);
            eOspaLmb{j}(t, i) = mean(eOspa, 2);
            hOspaLmb{j}(t, i) = mean(hOspa, 2);
        end
        %% LMBM filters
        for j = 1:numberOfLmbmAssociationMethods
            fprintf(['LMBM: ' lmbmDataAssociationMethods{j} '\n'])
            model.dataAssociationMethod = lmbmDataAssociationMethods{j};
            stateEstimates = runLmbmFilter(model, measurements);
            [eOspa, hOspa] = computeSimulationOspa(model, groundTruthRfs, stateEstimates);
            eOspaLmbm{j}(t, i) = mean(eOspa, 2);
            hOspaLmbm{j}(t, i) = mean(hOspa, 2);
        end

    end
    %% Save the current results (loadshedding)
    trialIndex = t;
    save(filename, 'eOspaLmb', 'eOspaLmbm', 'hOspaLmb', 'hOspaLmbm', 'trialIndex', '-append')
end
%% Plot error metrics
figure();
% E-OSPA
subplot(211); hold on; grid on; box on;
% Lmb 
for i = 1:numberOfLmbAssociationMethods
    meanEOspaError = mean(eOspaLmb{i}, 1);
    stdEOspaError = std(eOspaLmb{i}, 0, 1);
    errorbar(numberOfClutterReturns, meanEOspaError, stdEOspaError, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Lmbm
for i = 1:numberOfLmbmAssociationMethods
    meanEOspaError = mean(eOspaLmbm{i}, 1);
    stdEOspaError = std(eOspaLmbm{i}, 0, 1);
    errorbar(numberOfClutterReturns, meanEOspaError, stdEOspaError, 'LineWidth', 1.5, 'Color', lmbmPlotColours{i});
end
% Plot control
xlim([min(numberOfClutterReturns) max(numberOfClutterReturns)]);
ylim([0 model.ospaParameters.eC + 0.5])
xlabel('Number of clutter returns');
ylabel('E-OSPA');
% H-OSPA
subplot(212); hold on; grid on; box on;
% Lmb 
for i = 1:numberOfLmbAssociationMethods
    meanHOspaError = mean(hOspaLmb{i}, 1);
    stdHOspaError = std(hOspaLmb{i}, 0, 1);
    errorbar(numberOfClutterReturns, meanHOspaError, stdHOspaError, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Lmbm
for i = 1:numberOfLmbmAssociationMethods
    meanHOspaError = mean(hOspaLmbm{i}, 1);
    stdHOspaError = std(hOspaLmbm{i}, 0, 1);
    errorbar(numberOfClutterReturns, meanHOspaError, stdHOspaError, 'LineWidth', 1.5, 'Color', lmbmPlotColours{i});
end
subplot(212); hold on; grid on; box on;
% Plot control
%ax = gca;
%ax.XTick = numberOfClutterReturns;
xlim([min(numberOfClutterReturns) max(numberOfClutterReturns)]);
ylim([0 model.ospaParameters.hC + 0.1])
xlabel('Number of clutter returns');
ylabel('H-OSPA');
%matlab2tikz([filename(1:end-4) 'OSPA.tex']);
