function [assignments,costs]= murtysAlgorithmWrapper(P0,m)
% MURTYSALGORITHMWRAPPER -- Determine the m most likely data association
% event using Murty's algorithm.
%   [assignments,costs] = murtysAlgorithmWrapper(P0,m)
%
%   This is Vo and Vo's implementation of Murty's algorithm. It determines
%   the m most likely data association events for a given cost matrix.
%   It is a wrapper for a mex function that appears to have a memory leak.
%   If this code is run for a long time, you will run out of memory. If you
%   have a large number of objects, you will run out of memory.
%
%   See also generateLmbAssociationMatrices,
%       generateLmbmAssociationMatrices, murtysAlgorithm
%
%   Inputs
%       P0 - (n, ell) array. A negative log likelihood matrix C populated in 
%           generateLmbAssociationMatrices and
%           generateLmbmAssociationMatrices.
%       m - integer. The number of assoc=iation events you wish to
%           generate.
%   Output
%       assignments - (m, n) array. An array of the m most likely
%           association events. The array's rows are the association events.
%       costs - (m, 1) array. The cost of each assignment determined using
%           the negative loglikelihood. This parameter is never used in our
%               code.

if m==0
    assignments= [];
    costs= [];
    return;
end
    
  n1 = size(P0,1);
  n2 = size(P0,2);

  % Padding blocks for dummy variables
  blk1 = -log(ones(n1,n1));

  P0 = [P0 blk1];

  % Make costs non-negative (required by 'assignmentoptimal')
  x = min(min(P0));
  P0 = P0 - x;
  
  % Murty
  [assignments, costs] = murtysAlgorithm(P0,m);

  % Restore correct costs to assignments
  costs = costs + (x.*sum(assignments>0,2))';
      
  % Strip dummy variables
  assignments = assignments(:,1:n1);

  % Dummy assignments are clutter
  assignments(assignments>n2) = 0;

      
end
    