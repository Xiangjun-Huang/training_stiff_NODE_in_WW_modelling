% NODE Create Mini-Batches Function
% The |createMiniBatch| function creates a batch of observations of the target 
% dynamics. It takes as input the total number of time steps of the ground truth 
% data |numTimesteps|, the number of consecutive time steps to be returned for 
% each observation |numTimesPerObs|, the number of observations |miniBatchSize|, 
% and the ground truth data |X|.

function [x0, targets] = createMiniBatch_NODE(numTimesteps,numTimesPerObs,miniBatchSize,X)

% Create batches of trajectories.
s = randperm(numTimesteps - numTimesPerObs, miniBatchSize);

x0 = dlarray(X(s, :));
targets = zeros([miniBatchSize size(X,2) numTimesPerObs]);

for i = 1:numTimesPerObs
    targets(1:miniBatchSize, :, i) = X(s+i, :);
end

end