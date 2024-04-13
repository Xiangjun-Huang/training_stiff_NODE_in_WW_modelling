%% an IVP (initial value problem) of ASM2d_N2O model by collocation training and further NODE training

% Xiangjun Huang, Brunel University London, Apr 2024

% The scenario is set for a CSTR running for 6 hours with fixed DO at 2 mg/l. refer "Development of an ASM2d-N2O model
% to describe N2O emissions in municipal WWTPs under dynamic conditions" (http://dx.doi.org/10.1016/j.cej.2017.10.119).

%% 0. Start from scratch
clear;
clc;
close all;

%% 1. Generate trajectory data by ASM2d_N2O model
% simulate yTrue, dyTrue trajectory data 
[t, yTrue, dyTrue] = simulate_data();

% Add noise to the yTrue with noiseStd (5%)
yNoise = noise_data(yTrue,noiseStd=0.01); % options: noiseStd = 0.05

% Smooth noise data
ySmooth = smooth_data(yNoise); % options: method = 'gaussian', window =50

% Visualise true and noise and smooth data
plot_comp(t, yNoise, ySmooth, yTrue, yLabels={'Noise','Smooth','True'}, showRMSE=false, showUnit=false); 

%% 2. Initialize NODE network parameters
NodeParameters = initialize_NODE_parameters(); % options: stateSize=24, hiddenSize = 50

%% 3. Normalisation parameters
% choose the dataset as training data
yTrain = yNoise;

global yMean yStd dyMean dyStd
% [yMean, yStd, dyMean, dyStd] = mean_std_from_train_data(t, yTrain); % options: fromSaved = false

%% 4. Collocation training

% Calculate collocated pair (y, y').
[yColl, dyColl] = collocate_data(yTrain,t*24,"EpanechnikovKernel"); 
[yMean, yStd, dyMean, dyStd] = mean_std_from_collocated(yColl, dyColl); % options: fromSaved = false

% plot_comp(t,yTrue, yColl, yLabels={'yTrue','yColl'},showUnit=false,showRMSE=true,showMAE=true,showR2=true);
% plot_comp(t,dyTrue, dyColl, yLabels={'dyTrue','dyColl'},isDy=true, showUnit=false);

% Collocation training 
NodeParameters = collocation_training(t,yColl,dyColl, NodeParameters);

yPred = prediction(t,NodeParameters);
plot_comp(t,yNoise,yPred,yLabels={'Noise','Pred'},showRMSE=true,showMAE=true,showR2=true,showUnit=true,figureTitle="Collocation | ");

yTrain = yColl;

%% 5. Further direct NODE training
NodeParameters = NODE_training(t,yTrain,NodeParameters);

yPred = prediction(t,NodeParameters);
plot_comp(t,yNoise,yPred,yLabels={'Noise','Pred'},showRMSE=true,showMAE=true,showR2=true,showUnit=true,figureTitle="Noise and Pred | ");