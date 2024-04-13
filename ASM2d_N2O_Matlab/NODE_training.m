%% Train the network using a custom training loop with training dataset by NODE method.
% For each iteration:
% * Construct a mini-batch data from the training data with the createMiniBatch_NODE function.
% * Evaluate the loss/gradients using matlab dlfeval and the customised modelLoss_NODE function.
% * Update the model parameters using the adamupdate function.
% * Update the training progress plot.

function neuralOdeParameters = NODE_training(t, y, neuralOdeParameters, opt)    
    arguments 
        t (:,1) double {mustBeNumeric};
        y (:,:) double {mustBeNumeric};
        neuralOdeParameters;

        opt.nodeTimeSteps (1,1) {mustBePositive, mustBeInteger } = 800;
        opt.numIter (1,1) {mustBePositive, mustBeInteger } = 1000; 
        opt.miniBatchSize (1,1) {mustBePositive, mustBeInteger } = 200;
        opt.plotFrequency (1,1) {mustBePositive, mustBeInteger } = 100;   
        
        % Specify options for Adam optimization.
        opt.gradDecay (1,1)  = 0.9;
        opt.sqGradDecay (1,1)  = 0.999;
        opt.learnRate (1,1)  = 0.01;
    end

    %% Initialize loss figure    
    [fLoss, lineLossTrain] = plot_loss();

    %% Initialize the |averageGrad| and |averageSqGrad| parameters for the Adam solver    
    averageGrad = [];
    averageSqGrad = [];

    %% custom training by direct NODE approach

    dt = t(2);
    timeSteps = (0:opt.nodeTimeSteps)*dt;
    
    start = tic;
    fData=figure;
    
    for iter = 1:opt.numIter
        
        % Create batch 
        [y0Batch, yTarget] = createMiniBatch_NODE(opt.nodeTimeSteps, opt.miniBatchSize, y);
    
        % Evaluate network and compute loss and gradients
        [loss,gradients] = dlfeval(@modelLoss_NODE,timeSteps,y0Batch,neuralOdeParameters,yTarget);
        
        % Update network 
        [neuralOdeParameters,averageGrad,averageSqGrad] = adamupdate(neuralOdeParameters,gradients,averageGrad,averageSqGrad,iter,...
            opt.learnRate,opt.gradDecay,opt.sqGradDecay);
        
        % Plot loss   
        titleTxt = "NODE - Elapsed: " + string(duration(0,0,toc(start),Format="hh:mm:ss"));
        plot_loss_update(fLoss,lineLossTrain,loss,iter,titleTxt);
        
        % Plot predicted vs. real 
        if mod(iter, opt.plotFrequency) == 0  || iter == 1
    
            yPred = prediction(t,neuralOdeParameters);
            plot_comp(t, y, yPred, f=figure(fData), yLabels={'yTrain','yPred'}, figureTitle="NODE training");
            drawnow
        end
    end
end

%% NODE Create Mini-Batches Function

function [y0Batch, yBatch] = createMiniBatch_NODE(numTimesPerObs,miniBatchSize,y)

    numTimesteps = size(y,1);
    s = randperm(numTimesteps - numTimesPerObs, miniBatchSize);
    
    y0Batch = dlarray(y(s, :));
    yBatch = zeros([miniBatchSize size(y,2) numTimesPerObs]);
    
    for i = 1:numTimesPerObs
        yBatch(1:miniBatchSize, :, i) = y(s+i, :);
    end
end

%% NODE model function

function y = model_NODE(tspan,y0,neuralOdeParameters)

    y = dlode45(@odeModel,tspan,y0,neuralOdeParameters,DataFormat="BC");
    y(y<0) = 0; % nonNegtive restriction

end

%% NODE Model Loss Function

function [loss,gradients] = modelLoss_NODE(tspan,y0,neuralOdeParameters,yTarget)

    % Compute predictions.
    yPred = model_NODE(tspan,y0,neuralOdeParameters);

    % Compute L1 loss.
    loss = l1loss(yPred,yTarget,NormalizationFactor="all-elements",DataFormat="BCT");

    % Scaled L1 loss (need yMean)
    % yWeight = [1 1 1 1 1 1 1 1 1 1 1 1 1 yMean(14)/10 yMean(15)/10 yMean(16)/10 1 1 1 1 1 1 1 1];
    % loss = l1loss(yPredict./yWeight,dlarray(yTarget)./yWeight,NormalizationFactor="all-elements",DataFormat="BCT");
        
    % Compute L2 loss.
    % loss = l2loss(X,targets,NormalizationFactor="batch-size",DataFormat="BCT");
    
    % Compute Huber loss.
    % loss = huber(X,targets,NormalizationFactor="batch-size",DataFormat="BCT");
    
    % Compute gradients.
    gradients = dlgradient(loss,neuralOdeParameters);

end

