function model = generateModel(clutterRate, detectionProbability, dataAssociationMethod, varargin)
% GENERATEMODEL -- Generates a structure containing all simulation info.
%   model = generateModel(clutterRate, detectionProbabilty);
%
%   Declares all simulation information, except the ground truth.
%
%   See also generateGroundTruth.
%
%   Inputs
%       clutterRate - integer. The number of expected clutter per
%           time-step. Clutter is assumed to generated by a Poisson point
%           process.
%       detectionProbability - double. The probability of the sensor
%           detecting a target-generated measurement.
%       dataAssociationMethod - string. Type of data association method for the filters:
%           'LBP', 'Gibbs', 'Murty', 'LBPFixed'. LBP does not apply to LMBM filter.
%           'LBPFixed' uses a fixed number of iterations and is used to
%            experimentally verify the filter's asymptotic computational
%            complexity.
%       scenarioType - char array. Optional input: Type of simulated scenario.
%           'Fixed' - Four fixed birth locations.
%           'Random' - Randomly generated birth locations.
%       numberOfBirthLocations - integer. Optional input: Number of birth 
%           locations for the 'Random' scenario. 
%
%   Output
%       model - struct. A struct with the fields declared in this function.

%% Very dodgy input checking
if (nargin > 3)
    % Get input
    if (ischar(varargin{1}))
        model.scenarioType = varargin{1};
        if (strcmp(model.scenarioType, 'Random'))
            if (nargin > 4)
                numberOfBirthLocations = varargin{2};
            else
                error('You must specify the number of birth locations');
            end
        end
    else
        error('You must specify a scenario using a char array.');
    end
    % Check input
    if ~(strcmp(model.scenarioType, 'Fixed') || (strcmp(model.scenarioType, 'Random')) || (strcmp(model.scenarioType, 'Coalescence')))
        error('Scenario type must be Fixed, Random, or Coalescence');
    end
else
    model.scenarioType = 'Fixed';
end
%% State and measurement space dimensions
model.xDimension = 4;
model.zDimension = 2;
%% Sampling period
model.T = 1;
%% Linear motion model
% State transition matrix
model.survivalProbability = 0.95; % Existing target survival probability - could be state dependent.
model.existenceThreshold = 1e-2;
model.A = [ eye(model.xDimension/2) model.T*eye(model.xDimension/2); ...
    zeros(model.xDimension/2) eye(model.xDimension/2)];
model.u = zeros(model.xDimension, 1);
% Process noise
r0 = 1; %0.1
model.R = r0*[ (1/3)*(model.T^3)*eye(model.xDimension/2) 0.5*(model.T^2)*eye(model.xDimension/2);
    0.5*(model.T^2)*eye(model.xDimension/2) model.T*eye(model.xDimension/2)];
%% Linear observation model
% Observation matrix
model.C = [ eye(model.zDimension) zeros(model.zDimension) ];
% Measurement noise
q0 = 3.^2;
model.Q = q0*eye(model.zDimension);
%% Detection probability
model.detectionProbability = detectionProbability; % Could be state dependent.
%% Observation space
model.observationSpaceLimits = 100*[-1 1; -1 1];
model.observationSpaceVolume = prod(model.observationSpaceLimits(:, 2) - model.observationSpaceLimits(:, 1));
model.clutterRate = clutterRate;
model.clutterPerUnitVolume = clutterRate/model.observationSpaceVolume;
%% Birth parameters
% Determine spawning locations
if (strcmp(model.scenarioType, 'Fixed'))
    % Four fixed birth locations
    model.numberOfBirthLocations = 4;
    birthLocations = [-80.0 -20.0 0.0 40.0;
                -20.0 80.0 0.0 -60.0;
                0.0 0.0 0.0 0.0;
                0.0 0.0 0.0 0.0];
else
    % Random birth locations
    model.numberOfBirthLocations = numberOfBirthLocations;
    birthLocations = zeros(model.xDimension, numberOfBirthLocations);
    birthLocations(1:2, :) = 0.5 * model.observationSpaceLimits(:, 1) + model.observationSpaceLimits(:, 2) .* rand(model.zDimension, numberOfBirthLocations);
end
model.birthLocationLabels = 1:model.numberOfBirthLocations;
model.rB = 0.03 * ones(model.numberOfBirthLocations, 1);
model.rBLmbm = 0.045 * ones(model.numberOfBirthLocations, 1); %0.3 for a higher number of clutter returns
model.muB = repmat({zeros(model.xDimension, 1)}, model.numberOfBirthLocations, 1);
model.SigmaB = repmat( {diag( 0.5 * ones(model.xDimension, 1)  ).^2}, model.numberOfBirthLocations, 1); % 10
% Copy matrix into cell
for i = 1:model.numberOfBirthLocations
    model.muB{i} = birthLocations(:, i);
end
%% Object struct
object.birthLocation = 0;
object.birthTime = 0;
object.r = 0;
object.numberOfGmComponents = 0;
object.w = zeros(0, 1);
object.mu = repmat({}, 0, 1);
object.Sigma = repmat({}, 0, 1);
object.trajectoryLength = 0;
object.trajectory = [];
object.timestamps = zeros(1, 0);
object = repmat(object, 0, 1);
%% LMBM trajectory structs
trajectory.birthLocation = 0;
trajectory.birthTime = 0;
% trajectory.r = [];
trajectory.trajectory = [];
trajectory.trajectoryLength = 0;
trajectory.timestamps = zeros(1, 0);
birthTrajectory = repmat(trajectory, 1, model.numberOfBirthLocations);
%% Birth object struct
birthParameters = repmat(object, model.numberOfBirthLocations, 1);
for i = 1:model.numberOfBirthLocations
    birthParameters(i).birthLocation = model.birthLocationLabels(i);
    birthParameters(i).birthTime = 0;
    birthParameters(i).r = model.rB(i);
    birthParameters(i).numberOfGmComponents = 1;
    birthParameters(i).w = ones(1, 1);
    birthParameters(i).mu = model.muB(i);
    birthParameters(i).Sigma = model.SigmaB(i);
    % Trajectory
    birthParameters(i).trajectoryLength = 0;
    birthParameters(i).trajectory = repmat(80 * ones(model.xDimension, 1), 1, 100);
    birthParameters(i).timestamps = zeros(1, 100);
    % LMBM Birth trajectory
    birthTrajectory(i).birthLocation = model.birthLocationLabels(i);
    birthTrajectory(i).trajectoryLength = 0;
    birthTrajectory(i).trajectory = repmat(80 * ones(model.xDimension, 1), 1, 100);
    birthTrajectory(i).timestamps = zeros(1, 100);
end
%% Hypothesis struct
hypotheses.birthLocation = zeros(0, 1);
hypotheses.birthTime = zeros(0, 1);
hypotheses.w = 1; % Hypothesis weight, lowercase sigma in the theory
hypotheses.r = zeros(0, 1);
hypotheses.mu = repmat({}, 0, 1);
hypotheses.Sigma = repmat({}, 0, 1);
%% Object structs
model.object = object;
model.birthParameters = birthParameters;
model.hypotheses = hypotheses;
%% LMBM trajectory structs
model.trajectory = repmat(trajectory, 1, 0);
model.birthTrajectory = birthTrajectory;
%% GM parameters
model.gmWeightThreshold = 1e-6;
model.maximumNumberOfGmComponents = 5;
%% Track parameters
model.minimumTrajectoryLength = 20;
%% Data association method
model.dataAssociationMethod = dataAssociationMethod;
%% Loopy belief propagation parameters
model.maximumNumberOfLbpIterations = 1e3;
model.lbpConvergenceTolerance = 1e-6;
%% Gibbs sampling parameters
model.numberOfSamples = 2.5e3;
%% Murty's algorithm assignments
model.numberOfAssignments = 25;
%% Number of posterior hypothesis
model.maximumNumberOfPosteriorHypotheses = 25;
%% Posterior hypothesis weight threshold
model.posteriorHypothesisWeightThreshold = 1e-3;
%% State extraction method
model.useEapOnLmbm = false;
%% OSPA parameters
model.ospaParameters.eC = 5; % Euclidean cut-off
model.ospaParameters.eP = 2;
model.ospaParameters.hC = 0.5; % Hellinger cut-off
model.ospaParameters.hP = 2;
end