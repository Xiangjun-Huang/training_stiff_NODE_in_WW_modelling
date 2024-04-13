% add white noise with mean =0 and std = noise_std

function yNoise = noise_data(y, opt)
    arguments
        y (:, :) double {mustBeNumeric};
        opt.noiseStd (1,1) double {mustBeNumeric} = 0.05;
    end
      
    % yMean = mean(y);
    % yNoise = y +  randn(size(y,1),1) * (opt.noiseStd * yMean);
    yNoise = y .* (1 + opt.noiseStd * randn(size(y)) );

end
