function [hypotheses, objectsLikelyToExist] = lmbmNormalisationAndGating(posteriorHypotheses, model)
% LMBMNORMALISATIONANDGATING -- Gate the number of posterior hypotheses, and their parameters.
%    [hypotheses, objectsLikelyToExist] = lmbmNormalisationAndGating(posteriorHypotheses, model)
%
%   This function discards unlikely posterior parameters, and discards
%   Bernoulli with low existence probabilities from each hypothesis.
%
%   See also runLmbmFilter
%
%   Inputs
%       posteriorHypotheses - struct. A struct containing posterior LMBM
%           hypotheses, but with unnormalised hypothesis weights.
%       model - struct. A struct with the fields declared in generateModel.
%
%   Output
%       hypotheses - struct. A struct containing likely posterior LMBM
%           hypotheses, but with normalised hypothesis weights.
%       objectsLikelyToExist - array. An array of boolean indicating which objects
%           have been kept, and which have been discarded.

%% Normalise posterior hypothesis weights
logW = [posteriorHypotheses.w];
maxW = max(logW, [], 2);
w = exp(logW - maxW) ./ sum(exp(logW - maxW));
%% Gate posterior hypotheses
likelyHypotheses = w > model.posteriorHypothesisWeightThreshold;
hypotheses = posteriorHypotheses(likelyHypotheses);
w = w(likelyHypotheses) ./ sum(w(likelyHypotheses));
%% Sort the hypotheses according to weight
[w, sortedIndices] = sort(w, 'descend');
hypotheses = hypotheses(sortedIndices);
%% If number of posterior parameters exceed the maximum, discard the hypotheses with lowest weights
numberOfHypotheses = numel(w);
if (numberOfHypotheses > model.maximumNumberOfPosteriorHypotheses)
    w = w(1:model.maximumNumberOfPosteriorHypotheses);
    w = w ./ sum(w);
    hypotheses = hypotheses(1:model.maximumNumberOfPosteriorHypotheses);
    numberOfHypotheses = model.maximumNumberOfPosteriorHypotheses;
end
%% Determine total existence probability for each object, and discard unlikely components from each posterior hypothesis
r = sum(w .* [hypotheses.r], 2);
objectsLikelyToExist = r > model.existenceThreshold;
for i = 1:numberOfHypotheses
    % Truncate components
    hypotheses(i).birthLocation = hypotheses(i).birthLocation(objectsLikelyToExist);
    hypotheses(i).birthTime = hypotheses(i).birthTime(objectsLikelyToExist);
    hypotheses(i).w = w(i);
    hypotheses(i).r = hypotheses(i).r(objectsLikelyToExist, :);
    hypotheses(i).mu = hypotheses(i).mu(objectsLikelyToExist);
    hypotheses(i).Sigma = hypotheses(i).Sigma(objectsLikelyToExist);
end
%% Finally, this might be necessary sometimes
if (numberOfHypotheses == 0)
    hypotheses = model.hypotheses;
end
end