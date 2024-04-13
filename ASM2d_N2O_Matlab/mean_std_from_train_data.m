%% Initialize yMean, yStd, dyMean, dyStd parameters for normalisation method from train data by difference quotient 

% default keyword arguments:
% t: (data_size, 1) 
% y: (data_size, :)
% fromSaved: false

function [yMean, yStd, dyMean, dyStd] = mean_std_from_train_data(t, y, opt)
    arguments
        t (:, 1) double {mustBeNumeric};
        y (:, :) double {mustBeNumeric};
        opt.fromSaved (1,1) logical {mustBeNumericOrLogical} = false;
    end

    if opt.fromSaved && isfile("globalMeanStd.mat")
        load data/globalMeanStd.mat yMean yStd dyMean dyStd;
    else
        yMean = mean(y);
        yStd = std(y,1);
        
        [~, dy] = gradient(y);
        dy = dy ./ gradient(t);
        
        dyMean = mean(dy);
        dyStd = std(dy,1);
    end
end