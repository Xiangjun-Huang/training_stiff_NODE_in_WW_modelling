global yMean yStd dyMean dyStd
yMean = mean(coll_y);
yStd = std(coll_y,1);
dyMean = mean(coll_dy);
dyStd = std(coll_dy,1);
save global_mean_std.mat yMean yStd dyMean dyStd;