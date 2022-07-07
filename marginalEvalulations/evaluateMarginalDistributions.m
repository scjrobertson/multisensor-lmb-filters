% EVALUATEMARGINALDISTRIBUTIONS -- This script compares the output the
% marginal distributions produced the various data association algorithms.
%

%% Admin
clc; close all;
setPath;
%% Open existing file
filename = 'blimples.mat';
fileExists = isfile(fullfile(cd, filename));
%% Load old simulation, if it exists
if (fileExists)
    load(filename);
else
    % Number of trials
    numberOfTrials = 500;
    numberOfGibbsSamples =  [1 10:10:490 500:500:1e4 1e4:1e4:1e5 1e5:1e5:1e6];
    numberOfMurtyAssignments = [1 10:10:490 500:500:1e4];
    % Simulated scenario
    numberOfObjects = 10;
    detectionProbability = 0.75;
    expectedNumberOfClutterReturns = 20;
    % Trial counter
    trialIndex = 0;
    % Average Error output
    m = numel(numberOfGibbsSamples);
    rKlG = zeros(numberOfTrials, m);
    WKlG = zeros(numberOfTrials, m);
    rHG = zeros(numberOfTrials, m);
    WHG = zeros(numberOfTrials, m);
    % Distinct samples
    numberOfDistinctSamples = zeros(numberOfTrials, m);
    % Murty's algorithm
    n = numel(numberOfMurtyAssignments);
    rKlM = zeros(numberOfTrials, n);
    WKlM = zeros(numberOfTrials, n);
    rHM = zeros(numberOfTrials, n);
    WHM = zeros(numberOfTrials, n);
    % Save all current variables
    save(filename);
end
%% Run the experiments
startingPoint = trialIndex;
for t = (startingPoint+1):numberOfTrials
    %% Generate association matrices for new scenario
    model = generateSimplifiedModel(numberOfObjects, detectionProbability, expectedNumberOfClutterReturns);
    associationMatrices = generateAssociationMatrices(model);
    %% Run LBP
    [rLbp, WLbp] = loopyBeliefPropagation(associationMatrices, model.lbpCovergenceTolerance, model.maximumNumberOfLbpIterations);
    %% Run the Gibbs sampler
    for i = 1:m
        fprintf('Trial %d, Gibbs experiment %d, Number of Gibbs samples %d\n', t, i, numberOfGibbsSamples(i));
        [rGibbs, WGibbs, VGibbs] = lmbGibbsSampling(associationMatrices, numberOfGibbsSamples(i));
        % Kullback-Leibler
        rKlG(t, i) = averageKullbackLeiblerDivergence([1 - rGibbs rGibbs], [1 - rLbp rLbp]);
        WKlG(t, i) = averageKullbackLeiblerDivergence(WGibbs, WLbp);
        % Hellinger
        rHG(t, i) = averageHellingerDistance([1 - rGibbs rGibbs], [1 - rLbp rLbp]);
        WHG(t, i) = averageHellingerDistance(WGibbs, WLbp);
        % Distinct samples
        numberOfDistinctSamples(t, i) = size(VGibbs, 1);
    end
    %% Run Murty's algorithm
    % It is (probably) more efficient to determine all assignments first, then work out the marginals for each subset of assignments
    for i = 1:n
        fprintf('Trial %d, Murtys experiment %d, Number of assignments %d\n', t, i, numberOfMurtyAssignments(i));
        [rMurty, WMurty] = lmbMurtysAlgorithm(associationMatrices, numberOfMurtyAssignments(i));
        % Kullback-Keibler
        rKlM(t, i) = averageKullbackLeiblerDivergence([1 - rMurty rMurty], [1 - rLbp rLbp]);
        WKlM(t, i) = averageKullbackLeiblerDivergence(WMurty, WLbp);
        % Hellinger
        rHM(t, i) = averageHellingerDistance([1 - rMurty rMurty], [1 - rLbp rLbp]);
        WHM(t, i) = averageHellingerDistance(WMurty, WLbp);
    end
    %% Save the current results (loadshedding)
    trialIndex = t;
    save(filename, 'rKlG', 'WKlG', 'rHG', 'WHG', 'rKlM', 'WKlM', 'rHM', 'WHM', 'numberOfDistinctSamples', 'trialIndex', '-append')
end
%% Average Errors
% Gibbs average error
rKlGError = mean(rKlG, 1);
WKlGError = mean(WKlG, 1);
rHGError = mean(rHG, 1);
WHGError = mean(WHG, 1);
% Murty's average errors
rKlMError = mean(rKlM, 1);
WKlMError = mean(WKlM, 1);
rHMError = mean(rHM, 1);
WHMError = mean(WHM, 1);
% Number of distinct samples
distinctSamples = round(mean(numberOfDistinctSamples, 1));
%% Plot KL divergence 
figure(); hold on; grid on; box on;
set(gca, 'YScale', 'log');
loglog(numberOfGibbsSamples, rKlGError, 'Color', 'r', 'LineWidth', 2);
loglog(numberOfGibbsSamples, WKlGError, 'Color', 'b', 'LineWidth', 2);
xlabel('Number of samples');
ylabel('log(D_{KL} (p || q)');
% legend('Existence probability', 'Marginal association probabilities');
% title('Gibbs: Average KLD errors');
matlab2tikz([filename(1:end-4) 'GibbsKLD.tex']);
%% Plot KL divergence 
figure(); hold on; grid on; box on;
set(gca, 'YScale', 'log');
loglog(numberOfMurtyAssignments, rKlMError, 'Color', 'r', 'LineWidth', 2);
loglog(numberOfMurtyAssignments, WKlMError, 'Color', 'b', 'LineWidth', 2);
xlabel('Number of assignments');
ylabel('log(D_{KL} (p || q)');
% legend('Existence probability', 'Marginal association probabilities');
% title('Murty: Average KLD errors');
matlab2tikz([filename(1:end-4) 'MurtyKLD.tex']);
%% Plot Hellinger distance
figure(); hold on; grid on; box on;
set(gca, 'YScale', 'log');
loglog(numberOfGibbsSamples, rHGError, 'Color', 'r', 'LineWidth', 2);
loglog(numberOfGibbsSamples, WHGError, 'Color', 'b', 'LineWidth', 2);
xlabel('Number of samples');
ylabel('log(H(p, q))');
% legend('Existence probability', 'Marginal association probabilities');
% title('Gibbs: Average Hellinger distance errors');
matlab2tikz([filename(1:end-4) 'GibbsHellinger.tex']);
%% Plot Hellinger distance
figure(); hold on; grid on; box on;
set(gca, 'YScale', 'log');
loglog(numberOfMurtyAssignments, rHMError, 'Color', 'r', 'LineWidth', 2);
loglog(numberOfMurtyAssignments, WHMError, 'Color', 'b', 'LineWidth', 2);
xlabel('Number of assignments');
ylabel('log(H(p, q))');
% legend('Existence probability', 'Marginal association probabilities');
% title('Murty: Average Hellinger distance errors');
matlab2tikz([filename(1:end-4) 'MurtyHellinger.tex']);
%% Distinct sample plot
figure(); hold on; grid on; box on;
set(gca, 'YScale', 'log');
loglog(numberOfGibbsSamples, distinctSamples, 'Color', 'k', 'LineWidth', 2);
xlabel('Number of samples');
ylabel('Number of distinct samples');
matlab2tikz([filename(1:end-4) 'DistinctSamples.tex']);
%% Average KLD
function kl = averageKullbackLeiblerDivergence(p, q)
    logPQ = log(p ./ q);
    logPQ(isinf(logPQ)) = 0;
    kl = mean(sum(p .* logPQ, 2));
end
%% Average Hellinger distance
function h = averageHellingerDistance(p, q)
    hDist = sqrt(1 - sum( sqrt(p .* q), 2));
    h = mean(hDist);
end