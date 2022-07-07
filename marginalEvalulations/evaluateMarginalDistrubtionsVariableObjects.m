% evaluateMarginalDistrubtionsVariableObjects -- This script compares the output the
% marginal distributions produced the various data association algorithms.
%

%% Admin
clc; close all;
setPath;
%% Open existing file
filename = 'generalExamplesMultipleObjects.mat';
fileExists = isfile(fullfile(cd, filename));
%% Load old simulation, if it exists
if (fileExists)
    load(filename);
else
    % Number of trials
    numberOfTrials = 500;
    numberOfGibbsSamples = [1 10:10:490 500:500:1e4 1e4:1e4:1e5 1e5:1e5:1e6];
    numberOfMurtyAssignments = [1 10:10:490 500:500:1e4];
    % Simulated scenario
    numberOfObjects = 10:10:30;
    numberOfObjectExperiments = numel(numberOfObjects);
    detectionProbability = 0.75;
    expectedNumberOfClutterReturns = 20;
    % Trial counter
    trialIndex = 0;
    % Average Error output
    m = numel(numberOfGibbsSamples);
    rKlG = repmat({zeros(numberOfTrials, m)}, 1, numberOfObjectExperiments);
    WKlG = repmat({zeros(numberOfTrials, m)}, 1, numberOfObjectExperiments);
    % Distinct samples
    numberOfDistinctSamples = repmat({zeros(numberOfTrials, m)}, 1, numberOfObjectExperiments);
    % Murty's algorithm
    n = numel(numberOfMurtyAssignments);
    rKlM = repmat({zeros(numberOfTrials, n)}, 1, numberOfObjectExperiments);
    WKlM = repmat({zeros(numberOfTrials, n)}, 1, numberOfObjectExperiments);
    % Save all current variables
    save(filename);
end
%% Run the experiments
startingPoint = trialIndex;
for t = (startingPoint+1):numberOfTrials
    for ell = 1:numberOfObjectExperiments
        %% Generate association matrices for new scenario
        model = generateSimplifiedModel(numberOfObjects(ell), detectionProbability, expectedNumberOfClutterReturns);
        associationMatrices = generateAssociationMatrices(model);
        %% Run LBP
        [rLbp, WLbp] = loopyBeliefPropagation(associationMatrices, model.lbpCovergenceTolerance, model.maximumNumberOfLbpIterations);
        %% Run the Gibbs sampler
        for i = 1:m
            fprintf('Trial %d, Gibbs experiment %d, Number of objects %d, Number of Gibbs samples %d\n', t, i, numberOfObjects(ell), numberOfGibbsSamples(i));
            [rGibbs, WGibbs, VGibbs] = lmbGibbsSampling(associationMatrices, numberOfGibbsSamples(i));
            % Kullback-Leibler
            rKlG{ell}(t, i) = averageKullbackLeiblerDivergence([1 - rGibbs rGibbs], [1 - rLbp rLbp]);
            WKlG{ell}(t, i) = averageKullbackLeiblerDivergence(WGibbs, WLbp);
            % Distinct samples
            numberOfDistinctSamples{ell}(t, i) = size(VGibbs, 1);
        end
        %% Run Murty's algorithm
        % It is (probably) more efficient to determine all assignments first, then work out the marginals for each subset of assignments
        for i = 1:n
            fprintf('Trial %d, Murtys experiment %d, Number of objects %d, Number of assignments %d\n', t, i, numberOfObjects(ell), numberOfMurtyAssignments(i));
            [rMurty, WMurty] = lmbMurtysAlgorithm(associationMatrices, numberOfMurtyAssignments(i));
            % Kullback-Keibler
            rKlM{ell}(t, i) = averageKullbackLeiblerDivergence([1 - rMurty rMurty], [1 - rLbp rLbp]);
            WKlM{ell}(t, i) = averageKullbackLeiblerDivergence(WMurty, WLbp);
        end
    end
    %% Save the current results (loadshedding)
    trialIndex = t;
    save(filename, 'rKlG', 'WKlG', 'rKlM', 'WKlM', 'numberOfDistinctSamples', 'trialIndex', '-append')
end
%% Plotting
PLOT_COLOURS = {'r', 'g', 'b'};
%% Plot Gibbs existence
figure(); hold on; grid on; box on;
for i = 1:numberOfObjectExperiments
    semilogy(numberOfGibbsSamples, mean(rKlG{i}, 1) + eps(), 'Color', PLOT_COLOURS{i}, 'LineWidth', 2);
end
set(gca, 'YScale', 'log');
xlabel('Number of samples');
ylabel('D_{KL} (p || q)');
matlab2tikz([filename(1:end-4) 'GibbsExistence.tex']);
%% Plot Gibbs association
figure(); hold on; grid on; box on;
for i = 1:numberOfObjectExperiments
    semilogy(numberOfGibbsSamples, mean(WKlG{i}, 1) + eps(), 'Color', PLOT_COLOURS{i}, 'LineWidth', 2);
end
set(gca, 'YScale', 'log');
xlabel('Number of samples');
ylabel('D_{KL} (p || q)');
matlab2tikz([filename(1:end-4) 'GibbsAssociation.tex']);
%% Plot Murty's existence
figure(); hold on; grid on; box on;
for i = 1:numberOfObjectExperiments
    semilogy(numberOfMurtyAssignments, mean(rKlM{i}, 1) + eps(), 'Color', PLOT_COLOURS{i}, 'LineWidth', 2);
end
set(gca, 'YScale', 'log');
xlabel('Number of assignments');
ylabel('D_{KL} (p || q)');
matlab2tikz([filename(1:end-4) 'MurtyExistence.tex']);
%% Plot Murty's association
figure(); hold on; grid on; box on;
for i = 1:numberOfObjectExperiments
    semilogy(numberOfMurtyAssignments, mean(WKlM{i}, 1) + eps(), 'Color', PLOT_COLOURS{i}, 'LineWidth', 2);
end
set(gca, 'YScale', 'log');
xlabel('Number of assignments');
ylabel('log(D_{KL} (p || q)');
matlab2tikz([filename(1:end-4) 'MurtyAssociation.tex']);
%% Distinct sample plot
figure(); hold on; grid on; box on;
set(gca, 'YScale', 'log');
for i = 1:numberOfObjectExperiments
    semilogy(numberOfGibbsSamples, mean(numberOfDistinctSamples{i}, 1) + eps(), 'Color', PLOT_COLOURS{i}, 'LineWidth', 2);
end
xlabel('Number of samples');
ylabel('Number of distinct samples');
matlab2tikz([filename(1:end-4) 'DistinctSamples.tex']);
%% Average KLD
function kl = averageKullbackLeiblerDivergence(p, q)
    logPQ = log(p ./ q);
    logPQ(isinf(logPQ)) = 0;
    kl = mean(sum(p .* logPQ, 2));
end