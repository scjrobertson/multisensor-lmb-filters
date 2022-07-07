function objects = aaLmbTrackMerging(measurementUpdatedDistributions, model)
%  AALMBTRACKMERGING -- Merge the objects' measurement-updated distributions together using the AA-fusion rule.
%   objects = aaLmbTrackMerging(measurementUpdatedDistributions, model)
%
%   Merge the objects' measurement-updated distributions together using the AA-fusion rule.
%
%   See also generateMultisensorModel, loopyBeliefPropagation, lmbGibbsSampling, lmbMurtysAlgorithm
%
%   Inputs
%       measurementUpdatedDistributions - (1, numberOfSensors) cell array. 
%           Each is an object struct produced by computePosteriorLmbSpatialDistributions.   
%       model - struct. A struct with the fields declared in generateMultisensorModel.
%
%   Output
%       objects - struct. A struct containing the posterior LMB's Bernoulli
%           components.

objects = measurementUpdatedDistributions{1};
for i = 1:numel(objects)
    %% Initialise the merge distribution
    objects(i).r = model.aaSensorWeights(1) * objects(i).r;
    objects(i).w = model.aaSensorWeights(1) * objects(i).w;
    %% Merge the remaining sensors
    for s = 2:model.numberOfSensors
        objects(i).r = objects(i).r + model.aaSensorWeights(s) * measurementUpdatedDistributions{s}(i).r;
        objects(i).w = horzcat(objects(i).w, model.aaSensorWeights(s) * measurementUpdatedDistributions{s}(i).w);
        objects(i).mu = horzcat(objects(i).mu, measurementUpdatedDistributions{s}(i).mu);
        objects(i).Sigma = horzcat(objects(i).Sigma, measurementUpdatedDistributions{s}(i).Sigma);
    end
    %% Sort the resulting Gaussian mixture, and reduce it if necessary
    [~, sortedIndices] = sort(objects(i).w, 'descend');
    numberOfGmComponents = numel(objects(i).w);
    sortedIndices = sortedIndices(1:min(model.maximumNumberOfGmComponents, numberOfGmComponents));
    objects(i).numberOfGmComponents = numel(sortedIndices);
    objects(i).w = objects(i).w(sortedIndices) ./ sum(objects(i).w(sortedIndices));
    objects(i).mu = objects(i).mu(sortedIndices);
    objects(i).Sigma = objects(i).Sigma(sortedIndices);
end


end