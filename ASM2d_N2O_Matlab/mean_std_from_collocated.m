%% Initialize yMean, yStd, dyMean, dyStd parameters for normalisation method from collocated dataset. 

% default keyword arguments:
% yColl: (data_size, :) 
% dyColl: (data_size, :)
% fromSaved: false

function [yMean, yStd, dyMean, dyStd] = mean_std_from_collocated(yColl, dyColl, opt)
    arguments
        yColl (:, :) double {mustBeNumeric};
        dyColl (:, :) double {mustBeNumeric};
        opt.fromSaved (1,1) logical {mustBeNumericOrLogical} = false;
    end

    if opt.fromSaved && isfile("globalMeanStd.mat")
        load data/globalMeanStd.mat yMean yStd dyMean dyStd;
    else
        yMean = mean(yColl);
        yStd = std(yColl,1);
        dyMean = mean(dyColl);
        dyStd = std(dyColl,1);
    end
end