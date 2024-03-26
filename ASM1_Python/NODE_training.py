import torch
import numpy as np
from torchdiffeq import odeint
# from torchdiffeq import odeint_adjoint as odeint # for adjoint method!

def NODE_training(t, true_y, model, n_iters=2000, batch_time=16, batch_size=512):

    # Training initialization
    ilr = 1e-1  # learning rate
    # lr_adjust = { 2: 5e-5, 4: 1e-5, 6: 5e-6, 8: 1e-6,
    #     10: 5e-7, 15: 1e-7, 20: 5e-8 }
    gamma = 0.5
    lr_adjust = {int(0.2 * n_iters): gamma * ilr, int(0.4 * n_iters): gamma ** 2 * ilr, int(0.5 * n_iters): gamma ** 3 * ilr, int(0.6 * n_iters): gamma ** 4 * ilr,
                int(0.7 * n_iters): gamma ** 5 * ilr, int(0.8 * n_iters): gamma ** 6 * ilr, int(0.9 * n_iters): gamma ** 7 * ilr}
    # lr_adjust = {int(0.2 * n_iters): gamma * ilr, int(0.8 * n_iters): gamma ** 2 * ilr}
    # print(lr_adjust)
    optimizer = torch.optim.Adam(model.parameters(), lr=ilr)
    # scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=15, gamma=0.1)  # optional learning rate scheduler

    loss_list = []  # record loss data
    lr_list = []  # record learning rate data
    show_all_loss = False

    grad_norm = []

    # Define the loss function: Mean Squared Error
    # The sum of the squares of the differences between prediction and ground truth
    loss_mse = torch.nn.MSELoss(reduction='mean')

    # Define a Huber Loss Parameter: Delta
    delta_loss = 0.03

    # Train iteration
    print("Starting NODE training.")

    for it in range(1, n_iters + 1):
        model.train(mode=True)
        optimizer.zero_grad()
        batch_y0, batch_t, batch_y = get_batch(batch_time, batch_size, true_y, t)
        # print('batch_y0, batch_t batch_y size:===', batch_y0.size(),batch_t.size(),batch_y.size())

        # torchDiffEq provided solvers
        pred_y = odeint(model, batch_y0, batch_t, method='dopri5')  # you can change solver here

        # SciPy solver
        # pred_y = odeint(model, batch_y0, batch_t, method='scipy_solver', options={"solver": "BDF"})

        # Stable solver (need to extend the torchdiffeq package!)
        # pred_y, error_dopri5, error_dopri8, error_ssolver = odeint(model, batch_y0, batch_t, method='ssolver')

        # define a weighted loss function，also called MAPE (mean absolute percentage error)
        # mask = batch_y != 0
        # loss = torch.mean(torch.abs((pred_y[mask] - batch_y[mask])/batch_y[mask]))

        # define RMSE as the loss function (L2)
        # loss = torch.sqrt(loss_mse(pred_y, batch_y))

        # define MAE as the loss function (L1)
        loss = torch.mean(torch.abs(pred_y - batch_y))
        # loss_L1 = nn.L1Loss(reduction='mean')
        # loss = loss_L1(pred_y, batch_y)
            
        # define MAE with scaling, refer to Kim: "stiff NODEs"
        # loss = torch.mean(torch.abs(pred_y / yScale - batch_y / yScale))
        # loss_L1 = nn.L1Loss(reduction='mean')
        # loss = loss_L1(pred_y / yScale, batch_y / yScale)

        # define MSE as the loss function
        # loss = loss_mse(pred_y, batch_y)

        # define a Huber loss function
        # loss = torch.mean(torch.where(torch.abs(pred_y - batch_y) < delta_loss, 0.5 * ((pred_y - batch_y) ** 2), delta_loss * torch.abs(pred_y - batch_y) - 0.5 * (delta_loss ** 2)))

        # define a log cosh loss function
        # loss = torch.mean(torch.log(torch.cosh(pred_y - batch_y)))

        loss.backward()

        loss_list.append(loss.item())
        lr_list.append(optimizer.param_groups[-1]['lr'])

        # ========================== gradient norm ==================
        total_norm = 0.0
        for p in model.parameters():
            param_norm = p.grad.detach().data.norm(2)
            total_norm += param_norm.item() ** 2
        total_norm = total_norm ** 0.5
        # print(total_norm)
        grad_norm.append(total_norm)

        # ==========================gradient norm end ==============

        optimizer.step()
        # scheduler.step()

        # ====================adjust lr========================

        if it in lr_adjust.keys():
            lr = lr_adjust[it]
            for param_group in optimizer.param_groups:
                param_group['lr'] = lr
            # print('Updating learning rate to {}'.format(lr))

        if show_all_loss:
            print('Iter {:04d} | Loss {:.7f}'.format(it, loss.item()))
            # plot_loss(loss_list)
        else:
            if it % 250 == 0:
                print('Iteration: ', it, '/', n_iters)
    return model, loss_list, grad_norm, lr_list

# Random Batch function
def get_batch(batch_time, batch_size, true_y, t):
    data_size = t.size(dim=t.ndim-1)
    s = torch.from_numpy(np.random.choice(np.arange(data_size - batch_time, dtype=np.int64), batch_size, replace=True))
    batch_y0 = true_y[s]
    batch_t = t[:batch_time]
    batch_y = torch.stack([true_y[s + i] for i in range(batch_time)], dim=0)
    return batch_y0, batch_t, batch_y