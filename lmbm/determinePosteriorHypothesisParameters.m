function posteriorHypotheses = determinePosteriorHypothesisParameters(V, L, posteriorParameters, priorHypothesis)
% DETERMINEPOSTERIORHYPOTHESISPARAMETERS -- Determine a parameters for a new set of posterior LMBM hypotheses.
%   posteriorHypotheses = determinePosteriorHypothesisParameters(V, L, posteriorParameters, priorHypothesis)
%
%   Determine a parameters for a new set of posterior LMBM hypotheses.
%   These hypotheses will have unnormalised hypothesis weights.
%
%   See also runLmbmFilter, generateLmbmAssociationMatrices, lmbmGibbsSampling
%
%   Inputs
%       V - array. An array of distinct association events, where each row of the
%           array is an association event. See lmbmGibbsSampling.
%       L - array. An array of marginal log likelihood ratios. See
%           generateLmbmGibbsMatrices.
%       posteriorParameters - struct. A struct whose fields are a
%           hypothesis posterior distribution parameters. See
%           generateLmbmGibbsMatrices.
%       priorHypothesis - struct. A struct containing the prior LMBM hypotheses' Bernoulli components.
%       model - struct. A struct with the fields declared in generateModel.
%
%   Output
%       posteriorHypotheses - struct. A struct containing posterior LMBM
%           hypotheses, but with unnormalised hypothesis weights.

%% Declare the out variables
numberOfObjects = numel(priorHypothesis.r);
eta = 1:numberOfObjects;
numberOfPosteriorHypotheses = size(V, 1);
priorHypothesis.r = posteriorParameters.r;
posteriorHypotheses = repmat(priorHypothesis, 1, numberOfPosteriorHypotheses);
%% Generate a set of posterior hypotheses
for i = 1:numberOfPosteriorHypotheses
    % Association event
    v = V(i, :);
    % Linear indices
    ell = numberOfObjects * v + eta;
    % Missed detection events
    generatedMeasurement = v > 0;
    % Hypothesis weight
    posteriorHypotheses(i).w = log(priorHypothesis.w) + sum(L(ell));
    % Existence probabilities
    posteriorHypotheses(i).r(generatedMeasurement, :) = 1;
    % Means
    posteriorHypotheses(i).mu = posteriorParameters.mu(ell);
    % Covariance
    posteriorHypotheses(i).Sigma(generatedMeasurement) = posteriorParameters.Sigma(generatedMeasurement);
end
end