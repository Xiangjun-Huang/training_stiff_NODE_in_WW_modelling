import numpy as np
import torch

## load data from external file in the directory
def load_data(file_name='ASM_CSTR/ASM1_NODE_package/data/true_y.csv'):
    read_data = np.loadtxt(fname=file_name, delimiter=",", skiprows=1)
    t = torch.tensor(read_data[:,0], dtype=torch.float32)
    true_y = torch.unsqueeze(torch.tensor(read_data[:,1:],dtype=torch.float32),dim=1)

    return t, true_y