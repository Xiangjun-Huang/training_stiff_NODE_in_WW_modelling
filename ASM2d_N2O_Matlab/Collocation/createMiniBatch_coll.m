% Collocation Create Mini-Batches Function
% The |createMiniBatch| function creates a batch of collocate y and dy to be 
% trained. shuffle input as a text "shuffle" will generate a miniBatchSize sequence 
% with each randomly selected pair, otherwise it will generate a randomly selected 
% first pair, followed by miniBatchSize seqence in its original order.

function [yInput, dyTargets] = createMiniBatch_coll(numTimesteps, miniBatchSize, y, dy, shuffle)

if shuffle == "shuffle"
    s = randperm(numTimesteps, miniBatchSize);
else
    s1 = randperm(numTimesteps - miniBatchSize, 1);
    s = s1:1:(s1+miniBatchSize);
end

yInput = dlarray(y(s,:));
dyTargets = dlarray(dy(s,:));

end