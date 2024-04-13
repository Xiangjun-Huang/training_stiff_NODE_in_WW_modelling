%% Smooth data function 
% method = "gaussian" (better than smoothingspline), "movmean", "movmedian", 
%          "gaussian", "lowess", "loess", "rlowess", "rloess", "sgolay" 
% window =50
% 
% !Important: try as possible as it can be to smooth the data firstly, as more 
% noise you bring into neural network, the more instable the training will be.

function ySmooth = smooth_data(yNoise, opt)
    arguments
        yNoise (:, :) double {mustBeNumeric};
        opt.method (1,:) char {mustBeMember(opt.method,{'movmean','movmedian',...
                               'gaussian','lowess','loess','rlowess','rloess','sgolay'})} = 'gaussian';

        opt.window (1,1) {mustBeNumeric} = 50;
    end

    ySmooth = smoothdata(yNoise, opt.method, opt.window);

end