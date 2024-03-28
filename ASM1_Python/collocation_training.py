import torch
import numpy as np

def collocation_training(t, y, dy, model, n_iters=2000, batch_time=16, batch_size=512):

    # Training initialization
    ilr = 1e-1  # learning rate

    gamma = 0.5
    lr_adjust = {int(0.2 * n_iters): gamma * ilr, int(0.4 * n_iters): gamma ** 2 * ilr, int(0.5 * n_iters): gamma ** 3 * ilr, int(0.6 * n_iters): gamma ** 4 * ilr,
                int(0.7 * n_iters): gamma ** 5 * ilr, int(0.8 * n_iters): gamma ** 6 * ilr, int(0.9 * n_iters): gamma ** 7 * ilr}
    # lr_adjust = {int(0.2 * n_iters): gamma * ilr, int(0.8 * n_iters): gamma ** 2 * ilr}
    # print(lr_adjust)
    optimizer = torch.optim.Adam(model.parameters(), lr=ilr)

    loss_list = []  # record loss data
    grad_norm = [] # record gardient norm data
    lr_list = []  # record learning rate data
    show_all_loss = False

    # Define the loss function: Mean Squared Error
    # The sum of the squares of the differences between prediction and ground truth
    loss_mse = torch.nn.MSELoss(reduction='mean')

    # Define a Huber Loss Parameter: Delta
    delta_loss = 0.03

    # Train iteration
    print("Starting collocation training.")

    for it in range(1, n_iters + 1):
        model.train(mode=True)
        optimizer.zero_grad()
        batch_y, batch_dy = collocation_batch(batch_size, y, dy)
        # print('batch_y0, batch_t batch_y size:===', batch_y0.size(),batch_t.size(),batch_y.size())

        pred_dy = model(t, batch_y)

        # define a weighted loss function，also called MAPE (mean absolute percentage error)
        # mask = batch_dy != 0
        # loss = torch.mean(torch.abs((pred_dy[mask] - batch_dy[mask])/batch_dy[mask]))

        # define RMSE as the loss function (L2)
        # loss = torch.sqrt(loss_mse(pred_dy, batch_dy))

        # define MAE as the loss function (L1)
        loss = torch.mean(torch.abs(pred_dy - batch_dy))
        # loss_L1 = nn.L1Loss(reduction='mean')
        # loss = loss_L1(pred_dy, batch_dy)
            
        # define MSE as the loss function
        # loss = loss_mse(pred_dy, batch_dy)

        # define a Huber loss function
        # loss = torch.mean(torch.where(torch.abs(pred_dy - batch_dy) < delta_loss, 0.5 * ((pred_dy - batch_dy) ** 2), delta_loss * torch.abs(pred_dy - batch_dy) - 0.5 * (delta_loss ** 2)))

        # define a log cosh loss function
        # loss = torch.mean(torch.log(torch.cosh(pred_dy - batch_dy)))

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
        #     # print('Updating learning rate to {}'.format(lr))

        if show_all_loss:
            print('Iter {:04d} | Loss {:.7f}'.format(it, loss.item()))
            # plot_loss(loss_list)
        else:
            if it % 250 == 0:
                print('Iteration: ', it, '/', n_iters)
    return model, loss_list, grad_norm, lr_list

# Random Batch function
def collocation_batch(batch_size, y, dy, shuffle=True):
    data_size = y.size(dim=0)
    if shuffle:
        s = torch.from_numpy(np.random.choice(np.arange(data_size, dtype=np.int64), batch_size, replace=True))
        batch_y = y[s]
        batch_dy = dy[s]
    else:
        s = torch.from_numpy(np.random.choice(np.arange(data_size - batch_size, dtype=np.int64), 1, replace=True))
        batch_y = torch.stack([y[s + i] for i in range(batch_size)], dim=0)
        batch_dy = torch.stack([dy[s + i] for i in range(batch_size)], dim=0)
    return batch_y, batch_dy