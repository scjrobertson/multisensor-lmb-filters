% EVALUATESMALLEXAMPLES -- This script compares the LBP algorithm's output
% to the exact marginals porduced by Murty's algorithm. Murty's algortihm
% is used to exhaustively generate all data association events for data
% association problems.
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
    numberOfTrials = 5000;
    numberOfObjects = 7;
    trialIndex = 0;
    % Errors
    rKl = zeros(numberOfTrials, numberOfObjects);
    WKl = zeros(numberOfTrials, numberOfObjects);
    rH = zeros(numberOfTrials, numberOfObjects);
    WH = zeros(numberOfTrials, numberOfObjects);
    % Save all current variables
    save(filename);
end
%% Run the experiments
startingPoint = trialIndex;
for n = (startingPoint+1):numberOfObjects
    for t = 1:numberOfTrials
        %% Generate association matrices for new scenario
        model = generateSimplifiedModel(n, 0.95, 0);
        associationMatrices = generateAssociationMatrices(model);
        %% Run LBP
        [rLbp, WLbp] = loopyBeliefPropagation(associationMatrices, model.lbpCovergenceTolerance, model.maximumNumberOfLbpIterations);
        %% Cache association events, Murty's algorithm should exhaust all association events
        if (t == 1)
            numberOfEvents = calculateNumberOfAssociationEvents(n, n);
            [~, ~, V] = lmbMurtysAlgorithm(associationMatrices, numberOfEvents);
            ell = n * V + (1:n);
            W = repmat(V, 1, 1, n+1) == reshape(0:n, 1, 1, n+1);
        end
        %% Determine marginals
        J = reshape(associationMatrices.L(ell), numberOfEvents, n);
        L = permute(sum(prod(J, 2) .* W, 1), [2 1 3]);
        Sigma = reshape(L, n, n+1);
        % Normalise
        Tau = (Sigma .* associationMatrices.R) ./ sum(Sigma, 2);
        %% Determine existence probabilities
        rMurty = sum(Tau, 2);
        %% Determine marginal association probabilities
        WMurty =  Tau ./ rMurty;
        %% Println
        fprintf('Trial %d, Murtys experiment %d, Number of assignments %d\n', n, t, numberOfEvents);
        %% Determine errors
        % Kullback-Keibler
        rKl(t, n) = averageKullbackLeiblerDivergence([1 - rMurty rMurty], [1 - rLbp rLbp]);
        WKl(t, n) = averageKullbackLeiblerDivergence(WMurty, WLbp);
        % Hellinger
        rH(t, n) = averageHellingerDistance([1 - rMurty rMurty], [1 - rLbp rLbp]);
        WH(t, n) = averageHellingerDistance(WMurty, WLbp);
    end
    %% Save the current results (loadshedding)
    trialIndex = n;
    save(filename, 'rKl', 'WKl', 'rH', 'WH', 'trialIndex', '-append')
end
%% Determine mean errors
rKlError = mean(rKl, 1);
WKlError = mean(WKl, 1) + eps();
rHError = mean(rH, 1);
WHError = mean(WH, 1) + eps();
%% Plot KL divergence
figure(); hold on; grid on; box on;
semilogy(1:numberOfObjects, rKlError, 'Color', 'r', 'LineWidth', 2);
semilogy(1:numberOfObjects, WKlError, 'Color', 'b', 'LineWidth', 2);
xlabel('Number of objects');
ylabel('log(D_{KL} (p || q)');
legend('Existence probability', 'Marginal association probabilities');
title('Murty: Average KLD errors');
% Set plot up
ax = gca;
ax.XTick = unique( round(ax.XTick) );
set(ax, 'YScale', 'log');
% matlab2tikz([filename(1:end-4) 'KLD.tex']);
%% Plot Hellinger distance
figure(); hold on; grid on; box on;
semilogy(1:numberOfObjects, rHError, 'Color', 'r', 'LineWidth', 2);
semilogy(1:numberOfObjects, WHError, 'Color', 'b', 'LineWidth', 2);
xlabel('Number of objects');
ylabel('log(H(p, q))');
legend('Existence probability', 'Marginal association probabilities');
title('Gibbs: Average Hellinger distance errors');
% Set plot up
ax = gca;
ax.XTick = unique( round(ax.XTick) );
set(ax, 'YScale', 'log');
% matlab2tikz([filename(1:end-4) 'Hellinger.tex']);
%% Calculate number of association events
function numberOfEvents = calculateNumberOfAssociationEvents(n, m)
numberOfEvents = 0;
for k = 0:min(n, m)
    numberOfEvents = numberOfEvents + factorial(k) * nchoosek(n, k) * nchoosek(m, k);
end
end
%% Average KLD
function kl = averageKullbackLeiblerDivergence(p, q)
logPQ = log(p ./ q);
logPQ(isinf(logPQ) | isnan(logPQ)) = 0;
kl = mean(sum(p .* logPQ, 2));
end
%% Average Hellinger distance
function h = averageHellingerDistance(p, q)
hDist = sqrt(1 - sum( sqrt(p .* q), 2));
h = mean(real(hDist));
end