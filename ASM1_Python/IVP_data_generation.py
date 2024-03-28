"""
This function generate the trajectory data of an IVP for simulating a CSTR for 6 hours 
based on Extended ASM1 model (from Sumo22) with an initial condition as shown below:
'S_I', 'S_S', 'X_I', 'X_S', 'X_BH', 'X_BA', 'X_P', 'S_O', 'S_NO', 'S_NH', 'S_ND', 'X_ND', 'S_ALK', 'S_N2', 'X_INORG"
20.2,  59.8,  58.8, 260.1,   21.,   0.1,    0.1,    2,     0,      23.0,   1.8,    7.8,    0.007,    16,      35

It take the default values of parameters from sumo22.

# t_span: unit: day, 0.25 day = 6/24 day = 6 hours
# Initial condition

# true_y0 = torch.tensor([[20.2, 59.8, 58.8, 260.1, 2552.0, 148, 449, 2, 0, 23.0, 1.8, 7.8, 0.007, 0, 35]])
# true_y0 = torch.tensor([[30.0, 63.6, 58.5, 224.4, 1800., 110.0, 401, 2, 0, 20.0, 1.2, 7.1, 0.005, 0, 30]])
# true_y0 = torch.tensor([[15.6, 44.9, 32.8, 183.4, 1500., 86.1, 359, 2, 0, 32.0, 4.5, 11.9, 0.014, 0, 45]])
# true_y0 = torch.tensor([[49.1, 88.9, 72.1, 320.7, 2352., 56.3, 357, 2, 0, 62.0, 7.4, 9.88, 0.012, 0, 65]])
# true_y0 = torch.tensor([[25.6, 35.8, 18.4, 126.1, 1346., 68.6, 321, 2, 0, 21.0, 1.2, 5.75, 0.004, 0, 29]])

"""

import torch
import numpy as np
from mechanistic_model import mechanisticODE
from default_stoichiometric_kinetic_value import default_stoichimetric_kinetic_value

from torchdiffeq import odeint
# from torchdiffeq import odeint_adjoint as odeint # for adjoint method!

def IVP_data_generation(data_size = 1000, t_span = 0.25, 
    t_y_name = ['time','S_I', 'S_S', 'X_I', 'X_S', 'X_BH', 'X_BA', 'X_P', 'S_O', 'S_NO', 'S_NH', 'S_ND', 'X_ND', 'S_ALK', 'S_N2', 'X_INORG'],
    true_y0 = torch.tensor([[20.2, 59.8, 58.8, 260.1, 2552.0, 148, 449, 2, 0, 23.0, 1.8, 7.8, 0.007, 0, 35]]),
    save_csv=False):

    noFeature = true_y0.size(dim=true_y0.ndim-1)

    t = torch.linspace(0., t_span, data_size)

    # Generating training data from mechanistic ODEs without autograd
    # print("Generating data.")
    stoiP, kineP, contunity_check = default_stoichimetric_kinetic_value()
    with torch.no_grad():
        model_mechanistic = mechanisticODE(noFeature, kineP, stoiP)
        true_y = odeint(model_mechanistic, true_y0, t, method='dopri5')
        
        true_dy =torch.zeros_like(true_y)
        for i in range(data_size):
            true_dy[i,:,:] = model_mechanistic(t,true_y[i,:,:])
    print("Data generated.")

    if save_csv:
        save_y = torch.hstack((t.unsqueeze(-1),torch.squeeze(true_y))).numpy()
        save_dy = torch.hstack((t.unsqueeze(-1),torch.squeeze(true_dy))).numpy()
        np.savetxt(fname='ASM_CSTR/ASM1_NODE_package/data/true_y.csv', X=save_y, fmt='%.6f', delimiter=',', header=','.join(t_y_name), comments='')
        np.savetxt(fname='ASM_CSTR/ASM1_NODE_package/data/true_dy.csv', X=save_dy, fmt='%.6f', delimiter=',', header=','.join(t_y_name), comments='')

    return t, true_y, true_dy