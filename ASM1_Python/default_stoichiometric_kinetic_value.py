import torch

def default_stoichimetric_kinetic_value():

    # Define stoichiometric and kinetic parameters at temperature 20 degree and neutral pH, following ASM1

    # stoichimetric parameters
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
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0]],dtype = torch.float64)

    # kinetic parameter
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
                        ])

    # continuity matrix
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
                        ], dtype=torch.float64)

    # Continuity check #
    continuity = torch.matmul(stoiP, contiP)
    # print('Stoip \n', stoiP.dtype, '\n', stoiP)
    # print('ContiP \n', contiP)
    # print('continuity=\n',continuity)

    if torch.any(torch.abs(continuity) > 10 ** -15):
        print('Continuity check NOT pass!')
        continuity_check = False
    else:
        print('Continuity check pass!!!')
        continuity_check = True
    
    stoiP = stoiP.float()
    return stoiP, kineP, continuity_check