function hypothesis = lmbmPredictionStep(hypothesis, model, t)
% LMBMPREDICTIONSTEP -- Complete the LMBM filter's prediction step.
%   objects = lmbmPredictionStep(objects, model, t)
%
%   Computes predicted prior for the current time-step using the 
%   Chapman-Kolmogorov equation, assuming an LMBM prior and the standard 
%   multi-object motion model. 
%
%   See also runLmbmFilter, generateModel
%
%   Inputs
%       hypothesis - struct. A struct containing a posterior LMBM hypothesis' Bernoulli components.
%       model - struct. A struct with the fields declared in generateModel.
%       t - integer. An integer representing the simulation's current
%           time-step
%
%   Output
%       hypothesis - struct. A struct containing the prior LMBM hyopthesis' Bernoulli components.

%% Put existing Bernoulli componenents through the motion model
numberOfObjects = numel(hypothesis.r);
for i = 1:numberOfObjects
    hypothesis.r(i) = model.survivalProbability * hypothesis.r(i);
    hypothesis.mu{i} = model.A * hypothesis.mu{i} + model.u;
    hypothesis.Sigma{i} = model.A * hypothesis.Sigma{i} * model.A' + model.R;
end
%% Add in Bernoulli components for newly appearing objects
stride = (numberOfObjects + 1):(numberOfObjects + model.numberOfBirthLocations);
hypothesis.birthLocation(stride) = model.birthLocationLabels;
hypothesis.birthTime(stride) = t;
hypothesis.r(stride, :) = model.rBLmbm;
hypothesis.mu(stride) = model.muB;
hypothesis.Sigma(stride) = model.SigmaB;
end