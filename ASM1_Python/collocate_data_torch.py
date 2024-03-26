import numpy as np
import torch

def collocate_data_torch(t, true_y, kernel_str="EpanechnikovKernel"):
    """
    Computes a non-parametrically smoothed estimate of `y` and `y'` given the
    `data`.

    Args:
        t: A torch list of time points.
        true_y: A torch matrix with dimension (length(tpoints), 1, number_of_features).   
        kernel_str: The kernel function to use (default: "EpanechnikovKernel").

    Returns:
        y: The smoothed data.
        dy: The smoothed derivative of the data.
    """
    data = torch.squeeze(true_y).detach().numpy()
    tpoints = t.detach().numpy()

    nt = len(tpoints)
    nd = data.shape[0]

    if nt > nd:
        n = nd
        tpoints = tpoints[:n]
    elif nt < nd:
        n = nt
        data = data[:n,:]
    else:
        n = nt       

    y = np.zeros_like(data)
    dy = np.zeros_like(data)

    h = n**(-1 / 5) * n**(-3 / 35) * (np.log(n))**(-1 / 16)

    e1 = np.array([1, 0])
    e2 = np.array([0, 1, 0])

    for i in range(n):
        t1 = np.hstack((np.ones((n, 1)), np.reshape(tpoints - tpoints[i], (n, 1))))
        t2 = np.hstack((t1, np.reshape((tpoints - tpoints[i])**2, (n, 1))))
        
        w1 = np.zeros(n)
        for j in range(n):
            w1[j] = calculate_kernel(kernel_str,(tpoints[j]-tpoints[i])/h)/h
        
        w = np.diag(w1)

        # if not 'singular matrix'
        try:
            y[i, :] = e1 @ np.linalg.solve(t1.T @ w @ t1, t1.T) @ w @ data
            dy[i, :] = e2 @ np.linalg.solve(t2.T @ w @ t2, t2.T) @ w @ data
        # if is 'singular matrix', use pseudoinverse
        except np.linalg.LinAlgError:
            y[i, :] = e1 @ (np.linalg.pinv(t1.T @ w @ t1) @ t1.T) @ w @ data
            dy[i, :] = e2 @ (np.linalg.pinv(t2.T @ w @ t2) @ t2.T) @ w @ data      

        # or use least square optimisation
        # y[i, :] = e1 @ (np.linalg.lstsq(t1.T @ w @ t1, t1.T, rcond=None)[0]) @ w @ data
        # dy[i, :] = e2 @ (np.linalg.lstsq(t2.T @ w @ t2, t2.T, rcond=None)[0]) @ w @ data

    return torch.unsqueeze(torch.from_numpy(y),dim=1), torch.unsqueeze(torch.from_numpy(dy),dim=1)


def calculate_kernel(kernel_str, t):
    """
    Calculates the value of the specified kernel function.

    Args:
        kernel_str: The name of the kernel function.
        t: The input value of timepoints[i].

    Returns:
        The value of the kernel function.
    """

    kernel_functions = {
        "EpanechnikovKernel": lambda t: 0.75 * (1 - t**2) if abs(t) <= 1 else 0,
        "UniformKernel": lambda t: 0.5 if abs(t) <= 1 else 0,
        "TriangularKernel": lambda t: 1 - abs(t) if abs(t) <= 1 else 0,
        "QuarticKernel": lambda t: (15 * (1 - t**2)**2) / 16 if abs(t) <= 0 else 0,
        "TriweightKernel": lambda t: (35 * (1 - t**2)**3) / 32 if abs(t) <= 0 else 0,
        "TricubeKernel": lambda t: (70 * (1 - abs(t)**3)**3) / 80 if abs(t) <= 0 else 0,
        "CosineKernel": lambda t: (np.pi * np.cos(np.pi * t / 2)) / 4 if abs(t) <= 0 else 0,
        "GaussianKernel": lambda t: np.exp(-0.5 * t**2) / np.sqrt(2 * np.pi),
        "LogisticKernel": lambda t: 1 / (np.exp(t) + 2 + np.exp(-t)),
        "SigmoidKernel": lambda t: 2 / (np.pi * (np.exp(t) + np.exp(-t))),
        "SilvermanKernel": lambda t: (
            0.5 * np.sin(abs(t) / 2 + np.pi / 4) * np.exp(-abs(t) / np.sqrt(2))
        ),
    }

    if kernel_str not in kernel_functions:
        raise ValueError("Wrong kernel function name!")

    return kernel_functions[kernel_str](t)