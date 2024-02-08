% Specify options for Adam optimization.

gradDecay = 0.9;
sqGradDecay = 0.999;
learnRate = 0.01;
%% 
% Train for the specified iterations with a specified mini-batch-size.

numIter = 10000;
miniBatchSize = 200;
%% 
% Every specified iterations, solve the learned dynamics and display them against 
% the ground truth in a phase diagram to show the training path.

plotFrequency = 50;
% Train Model Using Custom Training Loop
% Initialize the training progress plot.

fLoss=figure;
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
% * Construct a mini-batch of data from the synthesized data with the |createMiniBatch_coll| 
% function.
% * Evaluate the model loss and gradients and loss using the |dlfeval| function 
% and the collocation |modelLoss_coll| function.
% * Update the model parameters using the |adamupdate| function.
% * Update the training progress plot.

start = tic;
fData=figure;

for iter = 1:numIter
    
    % Create batch 
    [yBatch, dyBatch] = createMiniBatch_coll(numTimeSteps, miniBatchSize, coll_y, coll_dy, "shuffle");

    % Evaluate network and compute loss and gradients
    [loss,gradients] = dlfeval(@modelLoss_coll, yBatch, dyBatch, neuralOdeParameters);
    
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

        % validate the solution 
        dyValidate = Model_coll(coll_y, neuralOdeParameters);

        figure(fData)        
        for k = 1:numFeature
            subplot(4,6,k)
            plot(tHour, dyValidate(:,k), tHour, coll_dy(:,k))
            xlim([0 6])
            xticks(0:2:6)
            xticklabels({'0','2','4','6h'})
            title(nameFeature(k)) 
            grid on
        end
        drawnow
    end
end
% Evaluate Model
% Use the model to compute approximated solutions with different initial conditions.
% 
% Define a new initial condition different from the one used for training the 
% model.

dyPred = Model_coll(coll_y,neuralOdeParameters);

dyPred = extractdata(dyPred);

rmse = sqrt(mean((coll_dy - dyPred).^2));

errCol = mean(abs(coll_dy - dyPred));

r = corrcoef(coll_dy, dyPred);
rc = r(2,1);

% Visualize Results
% Visualize the trained coll_dy and predicated dy  against the ground truth 
% solutions. 

% err = mean(abs(coll_dy - dyPred), "all");
% 
% f=figure;
% f.Position(3) = 2*f.Position(3);
% f.Position(4) = 1.8*f.Position(4);
% for k = 1:numFeature
%     subplot(4,6,k)
%     plot(tHour, dyPred(:,k), 'r-')
%     hold on
%     plot(tHour, coll_dy(:,k), 'b-.')
%     hold off
%     title(nameFeature(k)) 
%     xlim([0 6])
%     xticks(0:2:6)
%     xticklabels({'0','2','4','6h'})
%     ylabel(nameUnit(k)+'/d')
%     xlabel('')
%     if k==1 
%         legend('Train Y''','Collocate Y''','Location', 'southeast')
%     end
%     grid on
% end
% titleText = sprintf('%s %.4f %s %.0f%% %s',"Train Y' and Collocate Y' (MAE:", err,"), Y with", noiseSigma*100, "\sigma noise");
% sgtitle(titleText);