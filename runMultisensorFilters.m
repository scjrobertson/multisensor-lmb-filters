% RUNMULTISENSORFILTERS - Run the multi-sensor LMB or LMBM filters

%% Admin
close all; clc;
setPath;
%% Select type of filter
filterType = 'PU'; % 'IC', 'PU', 'LMBM'
%% Generate the model
numberOfSensors = 3;
clutterRates = [5 5 5];
detectionProbabilities = [0.67 0.70 0.73];
q = [4 3 2];
model = generateMultisensorModel(numberOfSensors, clutterRates, detectionProbabilities, q, 'PU', 'LBP', 'Fixed');
%% Generate observations
[groundTruth, measurements, groundTruthRfs] = generateMultisensorGroundTruth(model);
%% Run filters
if(strcmp(filterType, 'IC'))
    % Iterated-corrector LMB (IC-LMB) filter
    stateEstimates = runIcLmbFilter(model, measurements);
elseif(strcmp(filterType, 'PU'))
    % Parallel measurement update: PU-, GA-, or AA-LMB filters
    stateEstimates = runParallelUpdateLmbFilter(model, measurements);
else
    % Multisensor LMBM filter
    stateEstimates = runMultisensorLmbmFilter(model, measurements);
end
%% Plotting
plotMultisensorResults(model, measurements, groundTruth, stateEstimates, groundTruthRfs);
 