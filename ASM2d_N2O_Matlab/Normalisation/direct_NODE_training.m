%% custom training by direct NODE approach
% Define and Initialize Model Parameters
% The model function consists of a single call to |dlode45| to solve the ODE 
% defined by the approximated dynamics $f(t,y,\theta)$ for specified time steps.

neuralOdeTimesteps = 800;
dt = t(2);
timesteps = (0:neuralOdeTimesteps)*dt;
%% 
% Same learnable parameters |neuralOdeParameters| to use in the call to |dlode45|. 
% 
% |odeModel| takes as input the time input (unused), the corresponding solution, 
% and the ODE function parameters. 
% The function applies a fully connected operation, a gelu operation, and another 
% fully connected operation to the input data using the weights and biases given 
% by the parameters.
% Define NODE Model Function
% Create the function M|odel_NODE,| which computes the outputs of the deep learning 
% model. The function M|odel_NODE| takes as input the model parameters and the 
% input data. The function outputs the solution of the neural ODE.
% Define NODE Model Loss Function
% Create the function |modelLoss_NODE,| which takes as input the model parameters, 
% a mini-batch of input data with corresponding targets, and returns the loss 
% and the gradients of the loss with respect to the learnable parameters.


% Specify Training Options
% Specify options for Adam optimization.

gradDecay = 0.9;
sqGradDecay = 0.999;
learnRate = 0.01;
%% 
% Train for the specified iterations with a specified mini-batch-size.

numIter = 3000;
miniBatchSize = 200;
%% 
% Every specified iterations, solve the learned dynamics and display them against 
% the ground truth in a phase diagram to show the training path.

plotFrequency = 50;
% Train Model Using Custom Training Loop
% Initialize the training progress plot.

fLoss=figure;
% fLoss.Position(3) = 2*fLoss.Position(3);
C = colororder;
lineLossTrain = animatedline(Color=C(2,:));
set(gca, 'YScale', 'log')
ylim([0 inf])
xlabel("Iteration")
ylabel("Loss")
grid on
%% 
% Initialize the |averageGrad| and |averageSqGrad| parameters for the Adam solver.

averageGrad = [];
averageSqGrad = [];
%% 
% Train the network using a custom training loop.
% 
% For each iteration:
%% 
% * Construct a mini-batch of data from the synthesized data with the |createMiniBatch_NODE| 
% function.
% * Evaluate the model loss and gradients and loss using the |dlfeval| function 
% and the |modelLoss_NODE| function.
% * Update the model parameters using the |adamupdate| function.
% * Update the training progress plot.

numTrainingTimesteps = numTimeSteps;
trainingTimesteps = 1:numTrainingTimesteps;
plottingTimesteps = 2:numTimeSteps;

start = tic;
fData=figure;

for iter = 1:numIter
    
    % Create batch 
    [X, targets] = createMiniBatch_NODE(numTrainingTimesteps, neuralOdeTimesteps, miniBatchSize, ySmooth);

    % Evaluate network and compute loss and gradients
    [loss,gradients] = dlfeval(@modelLoss_NODE,timesteps,X,neuralOdeParameters,targets);
    
    % Update network 
    [neuralOdeParameters,averageGrad,averageSqGrad] = adamupdate(neuralOdeParameters,gradients,averageGrad,averageSqGrad,iter,...
        learnRate,gradDecay,sqGradDecay);
    
    % Plot loss
    figure(fLoss);
    currentLoss = double(loss);
    addpoints(lineLossTrain,iter,currentLoss);
    D = duration(0,0,toc(start),Format="hh:mm:ss");
    title("Elapsed: " + string(D))
    drawnow
    
    % Plot predicted vs. real dynamics
    if mod(iter,plotFrequency) == 0  || iter == 1

        % Use ode45 to compute the solution 
        % y = Model_NODE(t,dlarray(y0),neuralOdeParameters);
        y = dlode45(@odeModel,t,dlarray(y0),neuralOdeParameters,DataFormat="BC");
        y = squeeze(y)';

        figure(fData);
        for k = 1:numFeature
            subplot(4,6,k)
            plot(tHour, ySol(:,k),"r--")
            hold on
            plot(tHour(2:end), y(:,k), "b-")
            hold off
            title(nameFeature(k))
            xlim([0 6])
            xticks(0:2:6)
            xticklabels({'0','2','4','6h'})
            grid on
        end
        drawnow
    end
end

%% 

% Evaluate Model
% Use the model to compute approximated solutions with different initial conditions.
% 
% Define a new initial condition different from the one used for training the 
% model.

tPred = t;

y0Pred1 = [2,48.5,32.3,20,0,0,0,0,2.6,9,48.5,5,0,40.4,202,2500,250,70,100,200,100,189,50,220];
% y0Pred1 = [1.5,40,37,22,0,0,0,0,2.8,11,61,5.1,0,45,243,2760,2225,74,102,210,110,195,53,235];
% y0Pred1 = [1.2,32,28,15,0,0,0,0,2,7,39,4,0,34,189,2245,210,59,87,187,95,162,46,191];

% Numerically solve the ODE true dynamics with |ode15s| for the new initial 
% condition.

[~, yTrue1] = ode15s(trueModel, tPred, y0Pred1, odeOptions);
%% 
% Add same level of noise as training to the true solutions.

yNoise1 = yTrue1 .* (1 + noiseSigma * randn(size(yTrue1))); 

% Numerically solve the ODE with the learned neural ODE dynamics.

yPred1 = Model_NODE(tPred, dlarray(y0Pred1), neuralOdeParameters);
yPred1 = squeeze(yPred1)';
yPred1 = extractdata(yPred1);

tPred = tHour;
%% 
% Evaluate true y and predicted y

% rmse = sqrt(mean((yTrue1(2:end,:) - yPred1).^2));
% 
% errCol = mean(abs(yTrue1(2:end,:) - yPred1));
% 
% r = corrcoef(yTrue1(2:end,:), yPred1);
% rc = r(2,1);
%% 
% Evaluate noise y and predicted y

rmse = sqrt(mean((yNoise1(2:end,:) - yPred1).^2));

errCol = mean(abs(yNoise1(2:end,:) - yPred1));

r = corrcoef(yNoise1(2:end,:), yPred1);
rc = r(2,1);
% Visualize Predictions
% Visualize the predicted solutions for different initial conditions against 
% the ground truth solutions.
% True Y and Predicted Y

% err = mean(abs(yTrue1(2:end,:) - yPred1), "all");
% 
% figure;
% for k = 1:numFeature
%     subplot(4,6,k)
%     plot(tPred, yTrue1(:, k), "r--", tPred(2:end), yPred1(:,k), "b-", LineWidth=1)
%     title(nameFeature(k))
%     xlim([0 6])
%     xticks(0:2:6)
%     xticklabels({'0','2','4','6h'})
%     ylabel(nameUnit(k))
%     xlabel('')
%     if k==1 
%         legend('Observed Y','Predicted Y','Location', 'southeast')
%     end
%     grid on
% end
% 
% titleText = sprintf('%s %.0f%% %s %.4f %s',"Observed Y (", noiseSigma*100, "\sigma noise) and Predicted Y (MAE:", err,")");
% sgtitle(titleText);

% Noise Y and Predicted Y

err = mean(abs(yNoise1(2:end,:) - yPred1), "all");

figure;
for k = 1:numFeature
    subplot(4,6,k)
    plot(tPred, yNoise1(:, k), "r--", tPred(2:end), yPred1(:,k), "b-", LineWidth=1)
    title(nameFeature(k))
    xlim([0 6])
    ylim([0 inf])
    xticks(0:2:6)
    xticklabels({'0','2','4','6h'})
    ylabel(nameUnit(k))
    xlabel('')
    if k==1 
        legend('Observed Y','Predicated Y','Location', 'southeast')
    end
    grid on
end

titleText = sprintf('%s %.0f%% %s %.4f %s %.0f%% %s',"Observed Y (", noiseSigma*100, "\sigma noise) and Predicted Y (MAE:", err," RC: ", rc, ")");
sgtitle(titleText);