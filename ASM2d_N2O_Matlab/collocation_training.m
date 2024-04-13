%% Train the network using a custom training loop with collocated dataset.
% For each iteration:
% * Construct a mini-batch data from the training data with the createMiniBatch_coll function.
% * Evaluate the loss/gradients using matlab dlfeval and the customised modelLoss_coll function.
% * Update the model parameters using the adamupdate function.
% * Update the training progress plot.

function neuralOdeParameters = collocation_training(t, yColl, dyColl, neuralOdeParameters, opt)    
    arguments 
        t (:,1) double {mustBeNumeric};
        yColl (:,:) double {mustBeNumeric};
        dyColl (:,:) double {mustBeNumeric};
        neuralOdeParameters;

        opt.numIter (1,1) {mustBePositive, mustBeInteger } = 3000; 
        opt.miniBatchSize (1,1) {mustBePositive, mustBeInteger } = 200;
        opt.plotFrequency (1,1) {mustBePositive, mustBeInteger } = 100; 
        
        % options for Adam optimization.
        opt.gradDecay (1,1)  = 0.9;
        opt.sqGradDecay (1,1)  = 0.999;
        opt.learnRate (1,1)  = 0.01;
    end

    %% Initialize loss figure
    [fLoss, lineLossTrain] = plot_loss();

    %% Initialize the |averageGrad| and |averageSqGrad| parameters for the Adam solver.
    averageGrad = [];
    averageSqGrad = [];

    %% training loop    
    start = tic;
    fData=figure;
    
    for iter = 1:opt.numIter
        
        % Create batch 
        [yBatch, dyBatch] = createMiniBatch_coll(opt.miniBatchSize, yColl, dyColl, "shuffle");
    
        % Evaluate network and compute loss and gradients
        [loss,gradients] = dlfeval(@modelLoss_coll, yBatch, dyBatch, neuralOdeParameters);
        
        % Update network 
        [neuralOdeParameters,averageGrad,averageSqGrad] = adamupdate(neuralOdeParameters,gradients,averageGrad,averageSqGrad,iter, ...
            opt.learnRate,opt.gradDecay,opt.sqGradDecay);
        
        % Plot loss
        titleTxt = "Collocation - Elapsed: " + string(duration(0,0,toc(start),Format="hh:mm:ss"));
        plot_loss_update(fLoss,lineLossTrain,loss,iter,titleTxt);
        
        % Plot data
        if mod(iter,opt.plotFrequency) == 0  || iter == 1
    
            % evaluate the NN 
            dyPred = model_coll(yColl, neuralOdeParameters);
            plot_comp(t, dyColl, dyPred, f=figure(fData), yLabels={'dyColl','dyPred'}, isDy=true, figureTitle="Collocation training");
            drawnow
        end
    end
end

%% Collocation Create Mini-Batches Function

% The createMiniBatch function creates a batch of collocate y and dy to be 
% trained. shuffle input as a text "shuffle" will generate a miniBatchSize sequence 
% with each randomly selected pair, otherwise it will generate a randomly selected 
% first pair, followed by miniBatchSize seqence in its original order.

function [yInput, dyTargets] = createMiniBatch_coll(miniBatchSize, y, dy, shuffle)
    
    numTimesteps = size(y,1);
    if shuffle == "shuffle"
        s = randperm(numTimesteps, miniBatchSize);
    else
        s1 = randperm(numTimesteps - miniBatchSize, 1);
        s = s1:1:(s1+miniBatchSize);
    end
    
    yInput = dlarray(y(s,:));
    dyTargets = dlarray(dy(s,:));

end

%% Collocation Model Loss Function

% This function takes as inputs a batch of y and dy and net parameters.  it 
% computes the loss and the gradient of the loss with respect to the learnable 
% parameters of the network.

function [loss,gradients] = modelLoss_coll(y, dy, neuralOdeParameters)

    % Compute predictions.
    dyPred = model_coll(y, neuralOdeParameters);
    
    % Compute L1 loss.
    loss = l1loss(dyPred, dy, NormalizationFactor="all-elements",DataFormat="TC");
    
    % Compute L2 loss.
    % loss = l2loss(dyPred, dy, NormalizationFactor="batch-size",DataFormat="TC");
    
    % Compute Huber loss.
    % loss = huber(dyPred, dy, NormalizationFactor="batch-size",DataFormat="TC);
    
    % Compute gradients.
    gradients = dlgradient(loss,neuralOdeParameters);

end

%% Collocation model function

% y is a Matrix! The model_coll function, which defines the neural network used 
% to make predictions. simply output matrix from the input matrix of the network 
% for the purpose to reuse the network structure.
function dy = model_coll(y,neuralOdeParameters)

    dummyT = 0;
    dy =odeModel(dummyT, y, neuralOdeParameters);

end
