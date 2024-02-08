% Initialize the parameters structure. 

neuralOdeParameters = struct;
%% 
% Initialize the parameters for the fully connected operations in the ODE model. 
% The first fully connected operation takes as input a vector of size |stateSize| 
% and increases its length to |hiddenSize|. Conversely, the second fully connected 
% operation takes as input a vector of length |hiddenSize| and decreases its length 
% to |stateSize|.

stateSize = numFeature;
hiddenSize = 50;

neuralOdeParameters.fc1 = struct;
neuralOdeParameters.fc1.Weights = initialize_Glorot_V2(stateSize, hiddenSize);
neuralOdeParameters.fc1.Bias = initialize_Zeros_V2(hiddenSize);

neuralOdeParameters.fc2 = struct;
neuralOdeParameters.fc2.Weights = initialize_Glorot_V2(hiddenSize, hiddenSize);
neuralOdeParameters.fc2.Bias = initialize_Zeros_V2(hiddenSize);

neuralOdeParameters.fc3 = struct;
neuralOdeParameters.fc3.Weights = initialize_Glorot_V2(hiddenSize, hiddenSize);
neuralOdeParameters.fc3.Bias = initialize_Zeros_V2(hiddenSize);

neuralOdeParameters.fc4 = struct;
neuralOdeParameters.fc4.Weights = initialize_Glorot_V2(hiddenSize, stateSize);
neuralOdeParameters.fc4.Bias = initialize_Zeros_V2(stateSize);