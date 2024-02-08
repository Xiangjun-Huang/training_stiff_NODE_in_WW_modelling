% NODE Model Loss Function
% This function takes as inputs a vector |tspan|, a set of initial conditions 
% |X0|, the learnable parameters |neuralOdeParameters|, and target sequences |targets|. 
% It computes the predictions with the |model| function, and compares them with 
% the given targets sequences. Finally, it computes the loss and the gradient 
% of the loss with respect to the learnable parameters of the neural ODE.

function [loss,gradients] = modelLoss_NODE(tspan,X0,neuralOdeParameters,targets)

% Compute predictions.
X = Model_NODE(tspan,X0,neuralOdeParameters);

% Compute L1 loss.
loss = l1loss(X,targets,NormalizationFactor="all-elements",DataFormat="BCT");

% Compute L2 loss.
% loss = l2loss(X,targets,NormalizationFactor="batch-size",DataFormat="BCT");

% Compute Huber loss.
% loss = huber(X,targets,NormalizationFactor="batch-size",DataFormat="BCT");

% Compute gradients.
gradients = dlgradient(loss,neuralOdeParameters);

end