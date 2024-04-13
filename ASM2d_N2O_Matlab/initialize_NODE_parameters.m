%% Define and Initialize NODE network Parameters
% Define and initialize the parameters for the multilayer perceptron in the ODE model. 
% input/output: stateSize, hidden: hiddenSize, number of hidden layers: numHiddenLayer
% The initialization functions are attached to this file as supporting files. 

function neuralOdeParameters = initialize_NODE_parameters(opt)
    arguments
        opt.stateSize (1,1) {mustBePositive, mustBeInteger }  = 24;
        opt.hiddenSize (1,1)  {mustBePositive, mustBeInteger} = 50;
        opt.numHiddenLayer (1,1) {mustBePositive, mustBeInteger} = 2;
    end
    
    neuralOdeParameters = struct;

    % input layer
    neuralOdeParameters.fc1 = struct;
    neuralOdeParameters.fc1.Weights = initialize_Glorot_V2(opt.stateSize, opt.hiddenSize);
    neuralOdeParameters.fc1.Bias = initialize_Zeros_V2(opt.hiddenSize);
    
    % hidden layer
    for i = 1:opt.numHiddenLayer
        fld = "fc" + num2str(i+1);
        vlu = struct ('Weights',initialize_Glorot_V2(opt.hiddenSize, opt.hiddenSize),'Bias',initialize_Zeros_V2(opt.hiddenSize));
        neuralOdeParameters.(fld) = vlu;
    end

    % output layer
    fld = "fc" + num2str(opt.numHiddenLayer + 2);
    vlu = struct ('Weights',initialize_Glorot_V2(opt.hiddenSize, opt.stateSize),'Bias',initialize_Zeros_V2(opt.stateSize));
    neuralOdeParameters.(fld) = vlu;

end

%% Initialise weights

function weights = initialize_Glorot_V2(numIn, numOut,className)

    arguments
        numIn
        numOut
        className = 'single'
    end
    
    Z = 2*rand([numIn,numOut],className) - 1;
    bound = sqrt(6 / (numIn + numOut));
    
    weights = bound * Z;
    weights = dlarray(weights);

end

%% initialise bias

function parameter = initialize_Zeros_V2(numState,className)

    arguments
        numState
        className = 'single'
    end
    
    parameter = zeros([1, numState],className);
    parameter = dlarray(parameter);

end