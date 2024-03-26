import torch
import gc

def estimate_mean_std(t, true_y):
    noFeature = true_y.size(dim=true_y.ndim-1)

    diff_true_y = torch.diff(true_y, dim=0)

    yMu = torch.mean(true_y, dim=0)
    yStd = torch.std(true_y, dim=0)
    yStd[yStd == 0.0] = 1e-12
    dyMu = torch.mean(diff_true_y, dim=0)
    dyStd = torch.std(diff_true_y, dim=0)

    yMax = torch.max(true_y, dim=0).values
    yMin = torch.min(true_y, dim=0).values
    dyMax = torch.max(diff_true_y, dim=0).values
    dyMin = torch.min(diff_true_y, dim=0).values

    yScale = yMax - yMin
    yScale[yScale == 0.0] = yMax[yScale == 0.0]
    tScale = (t[-1]-t[0]).numpy()
    ytScale = yScale/tScale

    # release the memory of unused variable
    del diff_true_y
    gc.collect()

    # print('mu===\n',mu, '\nstd===\n',std, '\ndymu===\n',dymu, '\ndystd===\n',dystd)

    return noFeature, yMu, yStd, dyMu, dyStd, yMax, yMin, dyMax, dyMin, ytScale