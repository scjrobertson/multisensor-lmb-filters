function [eOspa, hOspa, cardinality] = computeSimulationOspa(model, groundTruthRfs, stateEstimateRfs)
% COMPUTESIMULATIONOSPA - Compute the OSPA between a simulation's ground
%   truth and a filter's state estimates
%    [eOspa, hOspa] = computeSimulationOspa(model, groundTruthRFS, stateEstimateRFS)
%
%   Compute the distance between the groundtruth RFS and a filter's RFS
%   estimates for an entire simulation using the OSPA metric. This returns
%   the Euclidean-OSPA (E-OSPA) and Hellinger-OSPA (H-OSPA) for each
%   time-step of the simulation.
%
%   See also generateModel, generateGroundtruth, Hungarian, ospa
%
%   Inputs
%       model - struct. A struct with the fields declared in generateModel.
%       groundTruthRfs - struct. Groundtruth RFS and optimal Kalman filter
%           output in RFS form.
%       stateEstimateRfs - struct. State estimates in RFS form.
%
%   Outputs
%       eOspa - (1, n) array. The E-OSPA components for each time-step.
%       hOspa - (1, n) array. The H-OSPA components for each time-step.
%       cardinality - (1, n) array. The cardinality of the state estimates.
%           Returned for convinience.

simulationLength = numel(groundTruthRfs.x);
eOspa = zeros(1, simulationLength);
hOspa = zeros(1, simulationLength);
cardinality = zeros(1, simulationLength);
for i = 1:simulationLength
    [euclidean, hellinger] = ospa(groundTruthRfs.x{i}, groundTruthRfs.mu{i}, groundTruthRfs.Sigma{i}, stateEstimateRfs.mu{i}, stateEstimateRfs.Sigma{i}, model.ospaParameters);
    eOspa(i) = euclidean(1);
    hOspa(i) = hellinger(1);
    cardinality(i) = numel(stateEstimateRfs.mu{i});
end

end