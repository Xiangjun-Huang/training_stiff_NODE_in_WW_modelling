% Collocation Model Loss Function
% This function takes as inputs a batch of y and dy and net parameters.  it 
% computes the loss and the gradient of the loss with respect to the learnable 
% parameters of the network.

function [loss,gradients] = modelLoss_coll(y, dy, neuralOdeParameters)

% Compute predictions.
dyPred = Model_coll(y, neuralOdeParameters);

% Compute L1 loss.
loss = l1loss(dyPred, dy, NormalizationFactor="all-elements",DataFormat="TC");

% Compute L2 loss.
% loss = l2loss(dyPred, dy, NormalizationFactor="batch-size",DataFormat="TC");

% Compute Huber loss.
% loss = huber(dyPred, dy, NormalizationFactor="batch-size",DataFormat="TC);

% Compute gradients.
gradients = dlgradient(loss,neuralOdeParameters);

end