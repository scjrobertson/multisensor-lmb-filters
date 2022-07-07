% LMBFILTERTIMETRIALS -- Determine the average run time for the LMB filter
% under various conditions.

%% NB Change model parameters to a "clean" model using the following:
% model.r0 = 0.1;
% model.q0 = 1;
% model.SigmaB = repmat( {diag( 2 * ones(model.xDimension, 1)  ).^2}, model.numberOfBirthLocations, 1);

%% Admin
clc; close all;
setPath;
%% Open existing file
filename = 'basicRuntimesLowBirthUncertainty.mat';
fileExists = isfile(fullfile(cd, filename));
%% Load old simulation, if it exists
if (fileExists)
    load(filename);
else
    % Number of trials
    numberOfTrials = 10;
    trialIndex = 0;
    % Clutter
    numberOfClutterReturns = [10:10:50 100:50:200];
    numberOfClutterExperiments = numel(numberOfClutterReturns);
    % Objects
    numberOfObjects = [10:10:50 100:50:150];
    numberOfObjectExperiments = numel(numberOfObjects);
    % Number of LBP iterations
    numberOfLbpIterations = 1e3:1e3:1e4;
    numberOfLbpIterationExperiments = numel(numberOfLbpIterations);
    % Number of LBP iterations
    numberOfGibbsSamples = [1000:1000:5000 1e4:1e4:1e5];
    numberOfGibbsSamplesExperiments = numel(numberOfGibbsSamples);
    % Number of LBP iterations
    numberOfMurtyAssignments = [10:10:50 100:50:250];
    numberOfMurtyAssignmentExperiments = numel(numberOfMurtyAssignments);
    % Simulation information
    simulationLength = 100;
    % Data assocation methods
    lmbDataAssociationMethods = {'LBP', 'Gibbs', 'Murty'};
    numberOfLmbAssociationMethods = numel(lmbDataAssociationMethods);
    % Run times
    clutterRuntimes = repmat({zeros(numberOfTrials, numberOfClutterExperiments)}, 1, numberOfLmbAssociationMethods);
    objectRuntimes = repmat({zeros(numberOfTrials, numberOfObjectExperiments)}, 1, numberOfLmbAssociationMethods);
    lbpIterationRuntimes = zeros(numberOfTrials, numberOfLbpIterationExperiments);
    gibbsSampleRuntimes = zeros(numberOfTrials, numberOfGibbsSamplesExperiments);
    murtyAssignmentRuntimes = zeros(numberOfTrials, numberOfMurtyAssignmentExperiments);
    % Plot colours
    lmbPlotColours = {'r', 'g', 'b'};
    % Save all current variables
    save(filename);
end
%% Run the experiments
startingPoint = trialIndex;
for t = (startingPoint+1):numberOfTrials
    %% Clutter trials
    for i = 1:numberOfClutterExperiments
        fprintf('Trial %d, Number of Clutter returns %d \n', t, numberOfClutterReturns(i));
        %% Generate new model and measurements
        model = generateModel(numberOfClutterReturns(i), 0.99, 'LBP', 'Fixed');
        [groundTruth, measurements] = generateGroundTruth(model);
        %% Data association methods
        for j = 1:numberOfLmbAssociationMethods
            fprintf([lmbDataAssociationMethods{j} '\n']);
            model.dataAssociationMethod = lmbDataAssociationMethods{j};
            tic;
            stateEstimates = runLmbFilter(model, measurements);
            runtime = toc;
            clutterRuntimes{j}(t, i) = runtime;
        end
    end
    %% Object trials
    for i = 1:numberOfObjectExperiments
        fprintf('Trial %d, Number of objects %d \n', t, numberOfObjects(i));
        %% Generate new model and measurements
        model = generateModel(5, 0.99, 'LBP', 'Random', numberOfObjects(i));
        [groundTruth, measurements] = generateGroundTruth(model, numberOfObjects(i));
        %% Data association methods
        for j = 1:numberOfLmbAssociationMethods
            fprintf([lmbDataAssociationMethods{j} '\n']);
            model.dataAssociationMethod = lmbDataAssociationMethods{j};
            tic;
            stateEstimates = runLmbFilter(model, measurements);
            runtime = toc;
            objectRuntimes{j}(t, i) = runtime;
        end
    end
    %% LBP iterations
    for i = 1:numberOfLbpIterationExperiments
        fprintf('Trial %d, Number of LBP iterations %d \n', t, numberOfLbpIterations(i));
        %% Generate new model and measurements
        model = generateModel(5, 0.99, 'LBPFixed');
        model.maximumNumberOfLbpIterations = numberOfLbpIterations(i);
        [groundTruth, measurements] = generateGroundTruth(model);
        tic;
        stateEstimates = runLmbFilter(model, measurements);
        runtime = toc;
        lbpIterationRuntimes(t, i) = runtime;
    end
    %% Gibbs samples
    for i = 1:numberOfGibbsSamplesExperiments
        fprintf('Trial %d, Number of Gibbs samples %d \n', t, numberOfGibbsSamples(i));
        %% Generate new model and measurements
        model = generateModel(5, 0.99, 'Gibbs');
        model.numberOfSamples = numberOfGibbsSamples(i);
        [groundTruth, measurements] = generateGroundTruth(model);
        tic;
        stateEstimates = runLmbFilter(model, measurements);
        runtime = toc;
        gibbsSampleRuntimes(t, i) = runtime;
    end
    %% Murty assignments
    for i = 1:numberOfMurtyAssignmentExperiments
        fprintf('Trial %d, Number of Murt assignments %d \n', t, numberOfMurtyAssignments(i));
        %% Generate new model and measurements
        model = generateModel(5, 0.95, 'Murty');
        model.numberOfAssignments = numberOfMurtyAssignments(i);
        [groundTruth, measurements] = generateGroundTruth(model);
        tic;
        stateEstimates = runLmbFilter(model, measurements);
        runtime = toc;
        murtyAssignmentRuntimes(t, i) = runtime;
    end
    %% Save the current results (loadshedding)
    trialIndex = t;
    save(filename, 'clutterRuntimes', 'objectRuntimes', 'lbpIterationRuntimes', 'gibbsSampleRuntimes', 'murtyAssignmentRuntimes', 'trialIndex', '-append');
end
%% Clutter plots
figure(); hold on; grid on; box on;
maxYValue = -inf;
% Lmb 
for i = 1:numberOfLmbAssociationMethods
    meanRuntime = mean(clutterRuntimes{i}, 1);
    stdRuntime = std(clutterRuntimes{i}, 0, 1);
    maxYValue = max(maxYValue, max(meanRuntime));
    errorbar(numberOfClutterReturns, meanRuntime, stdRuntime, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Plot control
% ax = gca;
% ax.XTick = numberOfClutterReturns;
xlim([min(numberOfClutterReturns) max(numberOfClutterReturns)]);
ylim([0 maxYValue + 0.5]);
xlabel('Number of clutter returns');
ylabel('Runtime (s)');
% matlab2tikz([filename(1:end-4) 'ClutterRuntimes.tex']);
%% Object plots
figure(); hold on; grid on; box on;
maxYValue = -inf;
% Lmb 
for i = 1:numberOfLmbAssociationMethods
    meanRuntime = mean(objectRuntimes{i}, 1);
    stdRuntime = std(objectRuntimes{i}, 0, 1);
    maxYValue = max(maxYValue, max(meanRuntime));
    errorbar(numberOfObjects, meanRuntime, stdRuntime, 'LineWidth', 1.5, 'Color', lmbPlotColours{i});
end
% Plot control
% ax = gca;
% ax.XTick = numberOfObjects;
xlim([min(numberOfObjects) max(numberOfObjects)]);
ylim([0 maxYValue + 0.5]);
xlabel('Number of objects');
ylabel('Runtime (s)');
% matlab2tikz([filename(1:end-4) 'ObjectsRuntimes.tex']);
%% LBP iterations
figure(); hold on; grid on; box on;
% Lmb 
meanRuntime = mean(lbpIterationRuntimes, 1);
stdRuntime = std(lbpIterationRuntimes, 0, 1);
errorbar(numberOfLbpIterations, meanRuntime, stdRuntime, 'LineWidth', 1.5, 'Color', lmbPlotColours{1});
% Plot control
%ax = gca;
%ax.XTick = numberOfLbpIterations;
xlim([min(numberOfLbpIterations) max(numberOfLbpIterations)]);
ylim([0 max(meanRuntime) + 0.5]);
xlabel('Number of LBP iterations');
ylabel('Runtime (s)');
% matlab2tikz([filename(1:end-4) 'LbpRuntimes.tex']);
%% Gibbs samples
figure(); hold on; grid on; box on;
% Lmb 
meanRuntime = mean(gibbsSampleRuntimes, 1);
stdRuntime = std(gibbsSampleRuntimes, 0, 1);
errorbar(numberOfGibbsSamples, meanRuntime, stdRuntime, 'LineWidth', 1.5, 'Color', lmbPlotColours{2});
% Plot control
%ax = gca;
%ax.XTick = numberOfGibbsSamples;
xlim([min(numberOfGibbsSamples) max(numberOfGibbsSamples)]);
ylim([0 max(meanRuntime) + 0.5]);
xlabel('Number of samples');
ylabel('Runtime (s)');
%% Murty assignments
figure(); hold on; grid on; box on;
% Lmb 
meanRuntime = mean(murtyAssignmentRuntimes, 1);
stdRuntime = std(murtyAssignmentRuntimes, 0, 1);
errorbar(numberOfMurtyAssignments, meanRuntime, stdRuntime, 'LineWidth', 1.5, 'Color', lmbPlotColours{3});
% Plot control
%ax = gca;
%ax.XTick = numberOfMurtyAssignments;
xlim([min(numberOfMurtyAssignments) max(numberOfMurtyAssignments)]);
ylim([0 max(meanRuntime) + 0.5]);
xlabel('Number of assignments');
ylabel('Runtime (s)');
