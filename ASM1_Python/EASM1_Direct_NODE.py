"""
This program simulate a CSTR for 6 hours based on Extended ASM1 model (from Sumo22) with an initial condition as shown below:
'S_I', 'S_S', 'X_I', 'X_S', 'X_BH', 'X_BA', 'X_P', 'S_O', 'S_NO', 'S_NH', 'S_ND', 'X_ND', 'S_ALK', 'S_N2', 'X_INORG"
 20.2,  59.8,  58.8, 260.1,   21.,   0.1,    0.1,    2,     0,      23.0,   1.8,    7.8,    0.007,    16,      35

The program uses mechanistic mode, neural ODE network and hybrid model.

The default values of parameters for mechanistic mode are taken from Sumo22.

The neural ODE network uses torchdiffeq package from Ricky TQ Chen, as same as claimed by Ward Quaghebeur
in his paper: "Hybrid differential equations_integrating mechanistic and data-driven techniques for modelling
of water systems".
"""

# Import packages
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
# import matplotlib.pyplot as plt
import gc
from torchmetrics.functional import r2_score

# Install torchdiffeq
# Import odeint with automatic differentiation or adjoint method
adjoint = False
if adjoint:
    from torchdiffeq import odeint_adjoint as odeint
else:
    from torchdiffeq import odeint

# If GPU acceleration is available
gpu = 0
global device
device = torch.device('cuda:' + str(gpu) if torch.cuda.is_available() else 'cpu')
if torch.cuda.is_available():
    torch.set_default_tensor_type('torch.cuda.FloatTensor')

import matplotlib
matplotlib.use('TkAgg',force=True)
from matplotlib import pyplot as plt
print("Switched to:",matplotlib.get_backend())

### Define model class ###

# EASM1 Mechanistic Model #
class mechanisticODE(nn.Module):

    def __init__(self, nf, k0, s0):
        super(mechanisticODE, self).__init__()
        self.nf = nf

        self.kinetic = nn.Parameter(k0)
        self.u_A = self.kinetic[0]
        self.b_A = self.kinetic[1]
        self.k_a = self.kinetic[2]
        self.K_O_A = self.kinetic[3]
        self.K_NH = self.kinetic[4]

        self.u_H = self.kinetic[5]
        self.n_g = self.kinetic[6]
        self.K_S = self.kinetic[7]
        self.b_H = self.kinetic[8]
        self.K_O_H = self.kinetic[9]
        self.K_NO = self.kinetic[10]
        self.K_NH_H = self.kinetic[11]

        self.k_h = self.kinetic[12]
        self.K_X = self.kinetic[13]
        self.n_h = self.kinetic[14]

        self.K_ALK = self.kinetic[15]

        self.stoiP = nn.Parameter(s0)

    def forward(self, t, y):
        S_I = y.view(-1, self.nf)[:, 0]
        S_S = y.view(-1, self.nf)[:, 1]
        X_I = y.view(-1, self.nf)[:, 2]
        X_S = y.view(-1, self.nf)[:, 3]
        X_BH = y.view(-1, self.nf)[:, 4]
        X_BA = y.view(-1, self.nf)[:, 5]
        X_P = y.view(-1, self.nf)[:, 6]
        S_O = y.view(-1, self.nf)[:, 7]
        S_NO = y.view(-1, self.nf)[:, 8]
        S_NH = y.view(-1, self.nf)[:, 9]
        S_ND = y.view(-1, self.nf)[:, 10]
        X_ND = y.view(-1, self.nf)[:, 11]
        S_ALK = y.view(-1, self.nf)[:, 12]
        S_N2 = y.view(-1, self.nf)[:, 13]
        X_INORG = y.view(-1, self.nf)[:, 14]

        rho = torch.zeros(8)
        rho[0] = self.u_H * (S_S / (self.K_S + S_S)) * (S_O / (self.K_O_H + S_O)) * (S_NH / (self.K_NH_H + S_NH)) * (S_ALK / (self.K_ALK + S_ALK)) * X_BH
        rho[1] = self.u_H * (S_S / (self.K_S + S_S)) * (self.K_O_H / (self.K_O_H + S_O)) * (S_NO / (self.K_NO + S_NO)) * (S_NH / (self.K_NH_H + S_NH)) * self.n_g * X_BH
        rho[2] = self.u_A * (S_NH / (self.K_NH + S_NH)) * (S_O / (self.K_O_A + S_O)) * (S_ALK / (self.K_ALK + S_ALK)) * X_BA
        rho[3] = self.b_H * X_BH
        rho[4] = self.b_A * X_BA
        rho[5] = self.k_a * S_ND * X_BH
        rho[6] = self.k_h * ((X_S / X_BH) / (self.K_X + (X_S / X_BH))) * (S_O / (self.K_O_H + S_O) + self.n_h * (self.K_O_H / (self.K_O_H + S_O)) * (S_NO / (self.K_NO + S_NO))) * X_BH
        rho[7] = rho[6] * (X_ND / X_S)

        r = torch.matmul(rho, self.stoiP)
        r[7] = 0  # keep oxygen at 2 mg/l or maneuver oxygen level as expected

        # dx.append(torch.stack((y,r),dim=1))
        return r.view(-1, 1, self.nf).to(device)


# Neural ODE network model #
class neuralODE(nn.Module):

    def __init__(self, nf, m, s, dm, ds, mx, mn, dmx, dmn, yts):
        super(neuralODE, self).__init__()
        self.nf = nf
        self.ymu = m
        self.ystd = s
        self.dymu = dm
        self.dystd = ds
        self.ymax = mx
        self.ymin = mn
        self.dymax = dmx
        self.dymin = dmn
        self.ytscale = yts

        self.net = nn.Sequential(
            nn.Linear(self.nf, 50),
            nn.GELU(),
            nn.Linear(50, 50),
            nn.GELU(),
            nn.Linear(50, 50),
            nn.GELU(),
            nn.Linear(50, self.nf))

        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, mean=0, std=0.5)

    def forward(self, t, y):

        # choose the methods: direct, standardisation, max-min normalisation, or equation scaling
        # direct output, without normalisation
        # y = self.net(y)

        # z-score standardisation
        y = (y - self.ymu) / self.ystd   
        y = self.net(y)
        y = y * self.dystd + self.dymu

        # max-min normalisation
        # ymax_ymin = self.ymax - self.ymin
        # y[:, ymax_ymin == 0.0] = 0.0
        # y[:, ymax_ymin != 0.0] = (y[:, ymax_ymin != 0.0] - self.ymin[ymax_ymin != 0]) / ymax_ymin[ymax_ymin != 0]
        # y = self.net(y)
        # y = y * (self.dymax - self.dymin) + self.dymin

        # equation scaling, refer to Kim:"stiff neural ordinary equations"
        # loss function must be matched with dividend : /yscale
        # y = self.net(y) 
        # y = y * self.ytscale

        return y


### Default value settings ###
torch.set_default_tensor_type(torch.DoubleTensor)

# Define number of components in the model (no. of features in the dataset, excluding t)
noFeature = 15

# Define stoichiometric and kinetic parameters at temperature 20 degree and neutral pH, following ASM1
Y_A = 0.24  # (g cell COD formed)/(g N oxidised)
Y_H = 0.67  # (g COD)/(g COD oxidised)
f_P = 0.08  # dimensionless
i_XB = 0.086  # (g N)/(g COD) in biomass
i_XP = 0.06  # (g N)/(g COD) in endogenous mass
i_COD_N2 = -12 / 7  # (g COD)/(g N) conversion factor for N2 into COD
i_COD_NO3 = -32 / 7  # (g COD)/(g N) conversion factor for NO3 into COD
i_NO3_N2 = 20 / 7  # (g COD)/(g N) conversion factor for NO3 reduction to N2
i_Charge_SNHx = 1 / 14000  # (kCharge)/(g N) conversion factor for NHx into charge
i_Charge_SNOx = -1 / 14000  # (kCharge)/(g N) conversion factor for NO3 into charge

stoiP = torch.tensor([[0, -1 / Y_H, 0, 0, 1, 0, 0, -(1 - Y_H) / Y_H, 0, -i_XB, 0, 0, -i_XB * i_Charge_SNHx, 0, 0],
                      [0, -1 / Y_H, 0, 0, 1, 0, 0, 0, -(1 - Y_H) / (i_NO3_N2 * Y_H), -i_XB, 0, 0, -(1 - Y_H) / (i_NO3_N2 * Y_H) * i_Charge_SNOx - i_XB * i_Charge_SNHx, (1 - Y_H) / (i_NO3_N2 * Y_H),
                       0],
                      [0, 0, 0, 0, 0, 1, 0, -(-i_COD_NO3 - Y_A) / Y_A, 1 / Y_A, -i_XB - 1 / Y_A, 0, 0, -(i_XB + 1 / Y_A) * i_Charge_SNHx + (1 / Y_A) * i_Charge_SNOx, 0, 0],
                      [0, 0, 0, 1 - f_P, -1, 0, f_P, 0, 0, 0, 0, i_XB - f_P * i_XP, 0, 0, 0],
                      [0, 0, 0, 1 - f_P, 0, -1, f_P, 0, 0, 0, 0, i_XB - f_P * i_XP, 0, 0, 0],
                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, i_Charge_SNHx, 0, 0],
                      [0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0]]).to(device)
kineP = torch.tensor([0.8,  # u_A, /day
                      0.15,  # b_A, /day
                      0.08,  # k_a, m3 /(g COD.day)
                      0.4,  # K_O_A, g O2/m3
                      1,  # K_NH, g N /m3

                      6,  # u_H, /day
                      0.8,  # n_g, dimensionless
                      20,  # K_S, g COD/m3
                      0.62,  # b_H, /day
                      0.2,  # K_O_H, g O2/m3
                      0.5,  # K_NO, g NO-N/m3
                      0.05,  # K_NH_H, g NH-N/m3

                      3,  # k_h, g X_S/(g COD X_BH.day)
                      0.03,  # K_X, g X_S/(g X_BH)
                      0.4,  # dimensionless

                      0.001  # eq/L
                      ]).to(device)

contiP = torch.tensor([[1., 0., 0.],
                       [1., 0., 0.],
                       [1., 0., 0.],
                       [1., 0., 0.],
                       [1., i_XB, 0.],
                       [1., i_XB, 0.],
                       [1., i_XP, 0.],
                       [-1., 0., 0.],
                       [i_COD_NO3, 1., i_Charge_SNOx],
                       [0., 1., i_Charge_SNHx],
                       [0., 1., 0.],
                       [0., 1., 0.],
                       [0., 0., -1.],
                       [i_COD_N2, 1, 0],
                       [0, 0, 0]
                       ]).to(device)
# Continuity check #
continuity = torch.matmul(stoiP, contiP)
# print('Stoip \n', stoiP.dtype, '\n', stoiP)
# print('ContiP \n', contiP)
# print('continuity=\n',continuity)

if torch.any(torch.abs(continuity) > 10 ** -15):
    print('Continuity check NOT pass!')
else:
    print('Continuity check pass!!!')

### Training part ###
# Training parameters and data generation
n_iters = 1000  # training iterations
data_size = 1000  # samples in dataset
batch_time = 16  # steps in batch
batch_size = 512  # samples per batch

# Initial condition, time span & parameters
y_name = ['S_I', 'S_S', 'X_I', 'X_S', 'X_BH', 'X_BA', 'X_P', 'S_O', 'S_NO', 'S_NH', 'S_ND', 'X_ND', 'S_ALK', 'S_N2', 'X_INORG']  # Name of components
true_y0 = torch.tensor([[20.2, 59.8, 58.8, 260.1, 2552.0, 148, 449, 2, 0, 23.0, 1.8, 7.8, 0.007, 0, 35]]).to(device)
# true_y0 = torch.tensor([[30.0, 63.6, 58.5, 224.4, 1800., 110.0, 401, 2, 0, 20.0, 1.2, 7.1, 0.005, 0, 30]]).to(device) # other initial conditions
# true_y0 = torch.tensor([[15.6, 44.9, 32.8, 183.4, 1500., 86.1, 359, 2, 0, 32.0, 4.5, 11.9, 0.014, 0, 45]]).to(device)
# true_y0 = torch.tensor([[49.1, 88.9, 72.1, 320.7, 2352., 56.3, 357, 2, 0, 62.0, 7.4, 9.88, 0.012, 0, 65]]).to(device)
# true_y0 = torch.tensor([[25.6, 35.8, 18.4, 126.1, 1346., 68.6, 321, 2, 0, 21.0, 1.2, 5.75, 0.004, 0, 29]]).to(device)

HRT = 0.25  # unit: day, 0.25 day = 6/24 day = 6 hours
t = torch.linspace(0., HRT, data_size).to(device)
kp0 = torch.tensor([0.8, 0.15, 0.08, 0.4, 1, 6, 0.8, 20, 0.62, 0.2, 0.5, 0.05, 3, 0.03, 0.4, 0.001]).to(device)  # manoeuvre to have an incomplete knowledge model

# Generating training data from mechanistic ODEs without autograd
# print("Generating data.")
with torch.no_grad():
    true_y = odeint(mechanisticODE(noFeature, kineP, stoiP), true_y0, t, method='dopri5')
print("Data generated.")

## Read training data from external file in the model/true_y_data.csv
# true_y=torch.tensor(np.loadtxt(fname='model/true_y_data.csv', delimiter=",", skiprows=0)).float().view(-1,1,noFeature)
# print('true_y =========== \n', true_y)
# np.savetxt(fname='model/true_y.csv', X=true_y.detach().numpy().reshape(-1, noFeature), fmt='%.3f', delimiter=',', header=','.join(y_name), comments='')


# plot true_y before adding noise
# fig = plt.figure(figsize=(12, 12))
# for i in range(noFeature):
#     ax = fig.add_subplot(4, 4, i + 1)
#     ax.plot((t.detach().cpu().numpy()) * 24, true_y[:, 0][:, i].detach().cpu().numpy(), label='True')
#     ax.set_title(y_name[i], pad=10)
#     ax.grid()
# lines, labels = fig.axes[-1].get_legend_handles_labels()
# fig.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.9, 0.1))
# fig.suptitle(str(batch_time) + ' steps of ' + str(batch_size) + ' batches for ' + str(n_iters) + ' iterations', fontsize=16, fontweight='bold')
# fig.tight_layout()
# plt.show()

# Add noise (mean = 0, std = 0.05)
# true_y *= (1 + torch.randn(data_size, 1, noFeature) / 20.)

# plot true_y after adding noise
# fig = plt.figure(figsize=(12, 12))
# for i in range(noFeature):
#     ax = fig.add_subplot(4, 4, i + 1)
#     ax.plot((t.detach().cpu().numpy()) * 24, true_y[:, 0][:, i].detach().cpu().numpy(), label='True')
#     ax.set_title(y_name[i], pad=10)
#     ax.grid()
# lines, labels = fig.axes[-1].get_legend_handles_labels()
# fig.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.9, 0.1))
# fig.suptitle(str(batch_time) + ' steps of ' + str(batch_size) + ' batches for ' + str(n_iters) + ' iterations', fontsize=16, fontweight='bold')
# fig.tight_layout()
# plt.show()
#
# print('true_y =========== \n', true_y)
# print('true_y size ===========', true_y.size())
# np.savetxt(fname='model/true_y_noise.csv', X=true_y.detach().numpy().reshape(-1, noFeature), fmt='%.3f', delimiter=',', header=','.join(y_name), comments='')

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
tScale = HRT
ytScale = yScale/tScale

# release the memory of unused variable
del diff_true_y
gc.collect()


# print('mu===\n',mu, '\nstd===\n',std, '\ndymu===\n',dymu, '\ndystd===\n',dystd)

# Random Batch function
def get_batch(batch_time, batch_size, data_size, true_y, t):
    s = torch.from_numpy(np.random.choice(np.arange(data_size - batch_time, dtype=np.int64), batch_size, replace=True)).to(device)
    batch_y0 = true_y[s]
    batch_t = t[:batch_time]
    batch_y = torch.stack([true_y[s + i] for i in range(batch_time)], dim=0).to(device)
    return batch_y0, batch_t, batch_y


# Training initialization
ilr = 1e-1  # learning rate
# lr_adjust = { 2: 5e-5, 4: 1e-5, 6: 5e-6, 8: 1e-6,
#     10: 5e-7, 15: 1e-7, 20: 5e-8 }
gamma = 0.5
lr_adjust = {int(0.2 * n_iters): gamma * ilr, int(0.4 * n_iters): gamma ** 2 * ilr, int(0.5 * n_iters): gamma ** 3 * ilr, int(0.6 * n_iters): gamma ** 4 * ilr,
             int(0.7 * n_iters): gamma ** 5 * ilr, int(0.8 * n_iters): gamma ** 6 * ilr, int(0.9 * n_iters): gamma ** 7 * ilr}
# lr_adjust = {int(0.2 * n_iters): gamma * ilr, int(0.8 * n_iters): gamma ** 2 * ilr}
# print(lr_adjust)
model = neuralODE(noFeature, yMu, yStd, dyMu, dyStd, yMax, yMin, dyMax, dyMin, ytScale).to(device)  # choose type of model to train (with signature neuralODE(), hybridODE(noFeature, kp0, stoiP))
optimizer = optim.Adam(model.parameters(), lr=ilr)
# scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=15, gamma=0.1)  # optional learning rate scheduler

loss_list = []  # record loss data
lr_list = []  # record learning rate data
show_all_loss = True

# Define the loss function: Mean Squared Error
# The sum of the squares of the differences between prediction and ground truth
loss_mse = nn.MSELoss(reduction='mean')

# Define a Huber Loss Parameter: Delta
delta_loss = 0.03

def plot_loss(loss_list):
	line1, = plt.plot(range(0,len(loss_list)),loss_list,'r.-')
	plt.legend(handles=[line1], labels=["train_loss"], loc="upper right", fontsize=7)
	plt.ylabel('log LOSS')
	# plt.show()

# Train iteration
print("Starting training.")

for it in range(1, n_iters + 1):
    model.train(mode=True)
    optimizer.zero_grad()
    batch_y0, batch_t, batch_y = get_batch(batch_time, batch_size, data_size, true_y, t)
    # print('batch_y0, batch_t batch_y size:===', batch_y0.size(),batch_t.size(),batch_y.size())

    # torchDiffEq provided solvers
    pred_y = odeint(model, batch_y0, batch_t, method='dopri5')  # you can change solver here

    # SciPy solver
    # pred_y = odeint(model, batch_y0, batch_t, method='scipy_solver', options={"solver": "BDF"})

    # Stable solver (need to extend the torchdiffeq package!)
    # pred_y, error_dopri5, error_dopri8, error_ssolver = odeint(model, batch_y0, batch_t, method='ssolver')

    # print('pred_y size:===', pred_y.size())
    # np.savetxt(fname='model/pred_y.csv', X=pred_y.detach().numpy().reshape(-1, noFeature), fmt='%.3f', delimiter=',', header=','.join(y_name), comments='')
    # np.savetxt(fname='model/batch_y.csv', X=batch_y.detach().numpy().reshape(-1, noFeature), fmt='%.3f', delimiter=',', header=','.join(y_name), comments='')

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

    # print('pred_y shape =====', pred_y.shape)
    # print('batch_y shape =====', batch_y.shape)

    # print ('pred_y ============== \n', pred_y)
    # print('batch_y ============== \n', batch_y)
    loss.backward()

    loss_list.append(loss.item())
    lr_list.append(optimizer.param_groups[-1]['lr'])

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

### Prediction ###

print('Starting testing...')
model.eval()
with torch.no_grad():
    pred_y = odeint(model, true_y0.view(1, 1, noFeature), t, method='dopri5').view(-1, 1, noFeature)

### Plot result ###
# plot data
fig = plt.figure(figsize=(12, 12))
for i in range(noFeature):
    ax = fig.add_subplot(4, 4, i + 1)
    ax.plot((t.detach().cpu().numpy()) * 24, pred_y[:, 0][:, i].detach().cpu().numpy(), '*', label='Pred')
    ax.plot((t.detach().cpu().numpy()) * 24, true_y[:, 0][:, i].detach().cpu().numpy(), label='True')
    ax.set_title(y_name[i], pad=10)
    ax.grid()
lines, labels = fig.axes[-1].get_legend_handles_labels()
fig.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.9, 0.1))
fig.suptitle(str(batch_time) + ' steps of ' + str(batch_size) + ' batches for ' + str(n_iters) + ' iterations', fontsize=16, fontweight='bold')
fig.tight_layout()
plt.show()

# Plot loss
fig, ax1 = plt.subplots(figsize=(8, 5))
ax2 = ax1.twinx()
ax1.plot(loss_list, 'b-')
ax2.plot(lr_list, 'g--')

ax1.set_ylabel('Log Loss', color='b')
ax1.set_yscale('log')
ax2.set_ylabel('Learning Rate', color='g')
ax2.set_yscale('log')
plt.show()
