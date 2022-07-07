function objects = computePosteriorLmbSpatialDistributions(objects, r, W, posteriorParameters, model)
% COMPUTEPOSTERIORLMBSPATIALDISTRIUBUTIONS -- Complete the LMB filter's measurement update
%    objects = computePosteriorLmbSpatialDistributions(objects, r, W, posteriorParameters, model)
%
%   This function computes each object's posterior spatial distrubtion. 
%
%   See also generateModel, runLmbFilter, lmbPredictionStep, 
%            loopyBeliefPropagation, generateLmbAssociationMatrices
%
%   Inputs
%       objects - struct. A struct containing the prior LMB's Bernoulli
%           components. This struct is produced by lmbPredictionStep.
%       r - array. Each object's posterior existence probability.
%       W - array. An array of marginal association probabilities, where
%           each row is an object's marginal association probabilities.
%       posteriorParameters - struct. A struct whose fields are an object's
%           posterior spatial distribution parameters.
%       model - struct. A struct with the fields declared in generateModel.
%
%   Output
%       objects - struct. A struct containing the posterior LMB's Bernoulli
%           components.

for i = 1:numel(objects)
    %% Update posterior existence probability
    objects(i).r = r(i);
    %% Reweight each the measurement-updated Gaussian mixtures using the marginal association probabilities
    numberOfPosteriorComponents = numel(posteriorParameters(i).w);
    posteriorWeights = reshape(W(i, :)' .* posteriorParameters(i).w, 1, numberOfPosteriorComponents);
    posteriorWeights = posteriorWeights ./ sum(posteriorWeights);
    %% Crude mixture reduction algorithm
    % Sort the weights
    [posteriorWeights, sortedIndices] = sort(posteriorWeights, 'descend');
    % Discard insignificant components
    significantComponents = posteriorWeights > model.gmWeightThreshold;
    significantWeights = posteriorWeights(significantComponents);
    objects(i).w = significantWeights ./ sum(significantWeights);
    sortedIndices = sortedIndices(significantComponents);
    objects(i).numberOfGmComponents = numel(objects(i).w);
    % Impose hard limit if there are too many components
    if (objects(i).numberOfGmComponents > model.maximumNumberOfGmComponents)
        objects(i).w = objects(i).w(1:model.maximumNumberOfGmComponents);
        objects(i).w = objects(i).w ./ sum(objects(i).w);
        sortedIndices = sortedIndices(1:model.maximumNumberOfGmComponents);
        objects(i).numberOfGmComponents = model.maximumNumberOfGmComponents;
    end
    %% Select the mixture components with the largest weights
    objects(i).mu = reshape(posteriorParameters(i).mu(sortedIndices), 1, objects(i).numberOfGmComponents);
    objects(i).Sigma = reshape(posteriorParameters(i).Sigma(sortedIndices), 1, objects(i).numberOfGmComponents);
end

end