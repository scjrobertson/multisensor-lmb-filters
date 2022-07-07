function objects = lmbPredictionStep(objects, model, t)
% LMBPREDICTIONSTEP -- Complete the LMB filter's prediction step.
%   objects = lmbPredictionStep(objects, model, t)
%
%   Computes predicted prior for the current time-step using the 
%   Chapman-Kolmogorov equation, assuming an LMB prior and the standard 
%   multi-object motion model. 
%
%   See also runLmbFilter, generateModel
%
%   Inputs
%       objects - struct. A struct containing the posterior LMB's Bernoulli components.
%       model - struct. A struct with the fields declared in generateModel.
%       t - integer. An integer representing the simulation's current
%           time-step
%
%   Output
%       objects - struct. A struct containing the prior LMB's Bernoulli components.

%% Put existing Bernoulli componenents through the motion model
numberOfObjects = numel(objects);
for i = 1:numberOfObjects
    objects(i).r = model.survivalProbability * objects(i).r;
    for j = 1:objects(i).numberOfGmComponents
        objects(i).mu{j} = model.A * objects(i).mu{j} + model.u;
        objects(i).Sigma{j} = model.A * objects(i).Sigma{j} * model.A' + model.R;
    end
end
%% Add in Bernoulli components for newly appearing objects
newNumberOfObjects = numberOfObjects + model.numberOfBirthLocations;
objects(numberOfObjects+1:newNumberOfObjects) = model.birthParameters;
[objects(numberOfObjects+1:newNumberOfObjects).birthTime] = deal(t);
end