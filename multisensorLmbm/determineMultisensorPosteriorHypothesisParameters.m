function posteriorHypotheses = determineMultisensorPosteriorHypothesisParameters(A, L, posteriorParameters, priorHypothesis)
% DETERMINEMULTISENSORPOSTERIORHYPOTHESISPARAMETERS -- Determine a parameters for a new set of posterior LMBM hypotheses.
%  posteriorHypotheses = determineMultisensorPosteriorHypothesisParameters(V, L, posteriorParameters, priorHypothesis)
%
%   Determine a parameters for a new set of posterior LMBM hypotheses.
%   These hypotheses will have unnormalised hypothesis weights.
%
%   See also runLmbmFilter, generateMultisensorLmbmAssociationMatrices, lmbmGibbsSampling
%
%   Inputs
%       A - array. An array of distinct association events, where each row of the
%           array is an association event. See lmbmGibbsSampling.
%       L - array. An array of marginal log likelihood ratios. See
%           generateLmbmGibbsMatrices.
%       posteriorParameters - struct. A struct whose fields are a
%           hypothesis posterior distribution parameters. See
%           generateLmbmGibbsMatrices.
%       priorHypothesis - struct. A struct containing the prior LMBM hypothesis' Bernoulli components.
%       model - struct. A struct with the fields declared in generateModel.
%
%   Output
%       posteriorHypotheses - struct. A struct containing posterior LMBM
%           hypotheses, but with unnormalised hypothesis weights.

%% Get sizes
d = size(L);
m = d(1:end-1) - 1;
numberOfObjects = d(end);
numberOfSensors = length(m);
%% Association matrices
eta = reshape(1:numberOfObjects, numberOfObjects, 1);
U = [zeros(numberOfObjects, numberOfSensors) eta];
A = A + 1;
%% Declare the output variables
numberOfPosteriorHypotheses = size(A, 1);
posteriorHypotheses = repmat(priorHypothesis, 1, numberOfPosteriorHypotheses);
%% Generate a set of posterior hypotheses
for i = 1:numberOfPosteriorHypotheses
    % Association event
    U(:, 1:numberOfSensors) = reshape(A(i, :), numberOfObjects, numberOfSensors);
    % Linear indices
    ell = determineLinearIndex(U, d);
    % Hypothesis weight
    posteriorHypotheses(i).w = log(priorHypothesis.w) + sum(L(ell));
    % Existence probabilities
    posteriorHypotheses(i).r = posteriorParameters.r(ell);
    % Means
    posteriorHypotheses(i).mu = posteriorParameters.mu(ell);
    % Covariance
    posteriorHypotheses(i).Sigma = posteriorParameters.Sigma(ell);
end
end
%% Determine a linear index using an association vector
function ell = determineLinearIndex(U, d)
ell = U(:, 1);
Pi = 1;
for i = 2:size(U, 2)
    Pi = Pi * d(i-1);
    ell  = ell  + Pi * (U(:, i) - 1);
end
end