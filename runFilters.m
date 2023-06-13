% RUNFILTERS - Run the single-sensor LMB or LMBM filters

%% Admin
close all; clc;
setPath;
%% Generate model
useLmbFilter = true; % Use LMB filter, or use LMBM filter
model = generateModel(10, 0.95, 'LBP', 'Fixed');
%% Generate observations
[groundTruth, measurements, groundTruthRfs] = generateGroundTruth(model);
%% Run a filter
if (useLmbFilter)
    stateEstimates = runLmbFilter(model, measurements);
else
    stateEstimates = runLmbmFilter(model, measurements);
end
%% Plotting
plotResults(model, measurements, groundTruth, stateEstimates, groundTruthRfs);
