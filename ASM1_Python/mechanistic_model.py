import torch.nn as nn
import torch 

### Define extended ASM1 mechanistic model class ###

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
        
        return r.view(-1, 1, self.nf)