%% Evaluate trained model, produce a prediction with a given initial condition by current neuralOdeParameters.

function yPred = prediction(t, neuralOdeParameters, opt)
    arguments
    t (:,:) double {mustBeNumeric};
    neuralOdeParameters;
    opt.y0 (1,:) double {mustBeNumeric}= [2,48.5,32.3,20,0,0,0,0,2.6,9,48.5,5,0,40.4,202,2500,250,70,100,200,100,189,50,220];
    %  [1.5,40,37,22,0,0,0,0,2.8,11,61,5.1,0,45,243,2760,2225,74,102,210,110,195,53,235];
    %  [1.2,32,28,15,0,0,0,0,2,7,39,4,0,34,189,2245,210,59,87,187,95,162,46,191];
    end
    yPred = dlode45(@odeModel,t,dlarray(opt.y0),neuralOdeParameters,DataFormat="BC");
    yPred(yPred<0) = 0; % nonNegtive restriction
    yPred = squeeze(yPred)';
    yPred = extractdata(yPred);
    yPred = [opt.y0; yPred];
end
