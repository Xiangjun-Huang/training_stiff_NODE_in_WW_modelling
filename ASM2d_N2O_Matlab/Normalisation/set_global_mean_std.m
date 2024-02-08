% calculate Mean and Std for input data sequence
% estimate Mean and Std for its derivatives by difference quotient

global yMean yStd dyMean dyStd

yMean = mean(ySmooth);
yStd = std(ySmooth);

dySmooth = gradient(ySmooth')./gradient(t');
dySmooth = dySmooth';

dyMean = mean(dySmooth);
dyStd = std(dySmooth);

save global_mean_std.mat yMean yStd dyMean dyStd;