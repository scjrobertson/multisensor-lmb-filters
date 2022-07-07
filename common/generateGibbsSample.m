function [v, w] = generateGibbsSample(P, v, w)
% GENERATEGIBBSSAMPLE -- Generate a new association event using Gibbs sampling
%   [v, w] = generateGibbsSample(P, v, w)
%
%   This function generates a new association event using Gibbs sampling, 
%   and the previous state of the association vectors v and w.
%
%   See also runLmbFilter, lmbGibbsSampling, lmbGibbsFrequencySampling
%            lmbmGibbsSampling
%
%   Inputs
%       P - array. An array of sampling probabilites. See
%           generateGibsAssociationMatrices.
%       v - array. The object-to-measurement association vector.
%       w - array. The measurement-to-object association vector.
%
%   Output
%       v - array. The updated object-to-measurement association vector.
%       w - array. The updated measurement-to-object association vector.

% Declare variables
[n, m] = size(P);
% Loop through each row of P
for i = 1:n
    % Start at first a_i^j that equals one in row i
    k = v(i) + (v(i) == 0);
    for j = k:m
        % Sample from a_i^j if column j is otherwise unnoccupied
        if ((w(j) == 0) || (w(j) == i))
            if (rand() < P(i, j))
                % Object i generated measurement z_j, we onto object i+1
                v(i) = j;
                w(j) = i;
                break;
            else
                % Object i does not exist, or missed its detection
                v(i) = 0;
                w(j) = 0;
            end
        end
    end
end
end