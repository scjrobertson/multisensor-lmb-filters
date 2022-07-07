function [eOspa, hOspa] = ospa(x, mu, Sigma, nu, T, p)
% OSPA - Optimal Subpattern Assigment
%   [eOspa, hOspa] = ospa(x, mu, Sigma, nu, T, p)
%
%   Compute the distance between the groundtruth RFS and a filter's RFS
%   estimates using the OSPA metric. This returns the Euclidean-OSPA (E-OSPA) and
%   Hellinger-OSPA (H-OSPA).
%
%   See also generateGroundtruth, Hungarian, computeSimulationOspa.
%
%   Inputs
%       x - cell. A finite set of groundtruth state estimates.
%       mu - cell array. A finite set of optimal Gaussian means.
%       Sigma - cell array. A finite set of optimal Gaussian covariances, in
%           mu's respective order.
%       nu - cell array. A finite set of Gaussian means produced by a filter.
%       T - cell array. A finite set of Gaussian covariances produced by a filter, in
%           nu's respective order.
%       p - struct. A struct whose fileds are the E- and H-OSPA parameters.
%
%   Outputs
%       eOspa - (3, 1) array. The E-OSPA components.
%       hOspa - (3, 1) array. The H-OSPA components.

%% Case 1 - Both sets are empty
if isempty(mu) && isempty(nu)
    eOspa = zeros(3, 1);
    hOspa = zeros(3, 1);
    return;
end
%% Case 2 - A single set is empty
if isempty(mu) || isempty(nu)
    eOspa = [p.eC; 0; p.eC];
    hOspa = [p.hC; 0; p.hC];
    return;
end
%% Case 3 - Non-empty sets
% Calculate sizes of the input sets
n =  numel(mu);
m = numel(nu);
ell = max(m, n);
q = abs(m-n);
% Populate the weight matrices
euclideanDistances = zeros(n, m);
hellingerDistances = zeros(n, m);
for i = 1:m
    for j = 1:n
        euclideanDistances(j, i) = norm(x{j} - nu{i});
        hellingerDistances(j, i) = computeHellingerDistance(mu{j}, Sigma{j}, nu{i}, T{i});
    end
end
% Cut off the matrices
euclideanCutOff = min(p.eC, euclideanDistances).^(p.eP);
hellingerCutOff = min(p.hC, hellingerDistances).^(p.hP);
%Compute optimal assignment and cost using the Hungarian algorithm
[~, euclideanCost] = Hungarian(euclideanCutOff);
[~, hellingerCost] = Hungarian(hellingerCutOff);
%Calculate final distance
eOspa(1, 1) = ((1/ell)*((p.eC^p.eP)*q + euclideanCost) ) ^(1/p.eP);
eOspa(2, 1) = ((1/ell)*euclideanCost)^(1/p.eP);
eOspa(3, 1) = ((1/ell)*(p.eC^p.eP)*q)^(1/p.eP);
%Calculate final distance
hOspa(1, 1) = ((1/ell)*((p.hC^p.hP)*q + hellingerCost) ) ^(1/p.hP);
hOspa(2, 1) = ((1/ell)*hellingerCost)^(1/p.hP);
hOspa(3, 1) = ((1/ell)*(p.hC^p.hP)*q)^(1/p.hP);
%% Hellinger distance
    function h = computeHellingerDistance(mu, Sigma, nu, T)
        Z = 0.5 * (Sigma + T);
        zeta = nu - mu;
        exponentialArgument = -0.125 * (zeta)' * (Z \ zeta) + 0.25 * log(det(Sigma)) + 0.25 * log(det(T)) - 0.5 * log(det(Z));
        h = sqrt(max(1 - exp(exponentialArgument), 0));
    end
end