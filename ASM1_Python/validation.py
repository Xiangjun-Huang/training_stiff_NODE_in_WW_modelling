### Validation ###
import torch
from torchdiffeq import odeint
# from torchdiffeq import odeint_adjoint as odeint # for adjoint method!

def validation(model, t, 
    true_y0 = torch.tensor([[20.2, 59.8, 58.8, 260.1, 2552.0, 148, 449, 2, 0, 23.0, 1.8, 7.8, 0.007, 0, 35]]),
    solver='dopri5'):
    print('Starting testing...')
    noFeature = true_y0.size(dim=true_y0.ndim-1)
    model.eval()
    with torch.no_grad():
        pred_y = odeint(model, true_y0.view(1, 1, noFeature), t, method=solver).view(-1, 1, noFeature)

    print('testing end')
    return pred_y