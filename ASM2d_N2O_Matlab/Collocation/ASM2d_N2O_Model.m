% This function defines ASM2d-N2O model. 
% t the unit depends on the unit of kinetic value, must be consistent.
% y denotes 24 components
% Sto denotes stoichimetric matrix
% Kine denotes kinetic value in a vector, must be in the defined order.

function dydt=ASM2d_N2O_Model(t,y,Sto,Kine)

Kine_cell = num2cell(Kine);

[K_H, K_O2_H, K_x_H, n_NO3_H, n_NO2_H, K_NO3_H, K_NO2_H, n_fe_H,...
 u_H, K_O2, K_F, K_NH4, K_P, K_ALK, K_A, K_NO3, K_NO2, n_NO3_D, q_fe, K_fe_H, b_H, u_H_Den,... 
 n_G3, n_G4, n_G5, K_S3, K_S4, K_S5, K_NO2_Den, K_OH4, K_N2O_Den, K_OH3, K_NO_Den, K_OH5, K_I3NO, K_I4NO, K_I5NO,...
 q_PHA, K_A_P, K_ALK_P, q_PP, K_O2_P, K_P_P, K_PHA_P, K_MAX_P, K_PP_P, K_IPP_P, K_PO4_P, n_NO3_P, n_NO2_P,... 
 K_NO3_P, K_NO2_P, u_PAO, b_PAO, b_PP, b_PHA,...
 u_AOB_HAO, q_AOB_AMO, K_O2_AOB1, K_NH4_AOB, K_O2_AOB2, K_NH2OH_AOB, q_AOB_HAO, K_NO_AOB_HAO, q_AOB_N2O_NN,...
 K_NO_AOB_NN, K_O2_AOB_ND, K_I_O2_AOB, K_HNO2_AOB, q_AOB_N2O_ND, K_ALK_AOB, K_P_AOB, u_NOB, K_O2_NOB, K_ALK_NOB,...
 K_NO2_NOB, K_P_NOB, b_AOB, b_NOB,...
 k_PRE, k_RED, K_ALK_PR] = Kine_cell{:};

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Part A: Hydrolysis processes
% K_H, K_O2_H, K_x_H, n_NO3_H, n_NO2_H, K_NO3_H, K_NO2_H, n_fe_H
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Part B: Heterotrophic biomass XH
% u_H, K_O2, K_F, K_NH4, K_P, K_ALK, K_A, K_NO3, K_NO2, n_NO3_D, q_fe, K_fe_H, b_H, u_H_Den,
% n_G3, n_G4, n_G5, K_S3, K_S4, K_S5, K_NO2_Den, K_OH4, K_N2O_Den, K_OH3, K_NO_Den, K_OH5, K_I3NO, K_I4NO, K_I5NO,
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Part C: PAO
% q_PHA, K_A_P, K_ALK_P, q_PP, K_O2_P, K_P_P, K_PHA_P, K_MAX_P, K_PP_P, K_IPP_P, K_PO4_P, n_NO3_P, n_NO2_P,
% K_NO3_P, K_NO2_P, u_PAO, b_PAO, b_PP, b_PHA
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Part D: Nitrifying Organisms
% u_AOB_HAO, q_AOB_AMO, K_O2_AOB1, K_NH4_AOB, K_O2_AOB2, K_NH2OH_AOB, q_AOB_HAO, K_NO_AOB_HAO, q_AOB_N2O_NN,
% K_NO_AOB_NN, K_O2_AOB_ND, K_I_O2_AOB, K_HNO2_AOB, q_AOB_N2O_ND, K_ALK_AOB, K_P_AOB, u_NOB, K_O2_NOB, K_ALK_NOB,
% K_NO2_NOB, K_P_NOB, b_AOB, b_NOB
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Part E: Precipitation of P with Fe(OH)3
% k_PRE, k_RED, K_ALK_PR
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% Process rate Rho calculation
% Part A: Hydrolysis processes
    Rho(1,1) = K_H*y(1)/(K_O2_H+y(1))*(y(15)/y(16))/(K_x_H+(y(15)/y(16)))*y(16);
    Rho(2,1) = K_H*n_NO3_H*K_O2_H/(K_O2_H+y(1))*y(9)/(K_NO3_H+y(9))*(y(15)/y(16))/(K_x_H+(y(15)/y(16)))*y(16);
    Rho(3,1) = K_H*n_NO2_H*K_O2_H/(K_O2_H+y(1))*y(8)/(K_NO2_H+y(8))*(y(15)/y(16))/(K_x_H+(y(15)/y(16)))*y(16);
    Rho(4,1) = K_H*n_fe_H*K_O2_H/(K_O2_H+y(1))*(K_NO2_H/(K_NO2_H+(y(9)+y(8))))*(y(15)/y(16))/(K_x_H+(y(15)/y(16)))*y(16);
% Part B: Heterotrophic biomass X_H
    Rho(5,1) = u_H*y(2)/(K_F+y(2))*y(2)/(y(2)+y(3))*y(1)/(K_O2+y(1))*y(4)/(K_NH4+y(4))*y(10)/(K_P+y(10))*y(12)/(K_ALK+y(12))*y(16);
    Rho(6,1) = u_H*y(3)/(K_A+y(3))*y(3)/(y(2)+y(3))*y(1)/(K_O2+y(1))*y(4)/(K_NH4+y(4))*y(10)/(K_P+y(10))*y(12)/(K_ALK+y(12))*y(16);
    Rho(7,1) = u_H*n_NO3_D*y(2)/(K_F+y(2))*y(2)/(y(2)+y(3))*K_O2/(K_O2+y(1))*y(9)/(K_NO3+y(9))*y(4)/(K_NH4+y(4))*y(10)/(K_P+y(10))*y(12)/(K_ALK+y(12))*y(16);
    Rho(8,1) = u_H*n_G3*y(2)/(K_S3+y(2))*y(2)/(y(2)+y(3))*y(8)/(K_NO2_Den+y(8))*K_OH3/(K_OH3+y(1))*y(4)/(K_NH4+y(4))*y(10)/(K_P+y(10))*y(12)/(K_ALK+y(12))*y(16);
    Rho(9,1) = u_H*n_G4*y(2)/(K_S4+y(2))*y(2)/(y(2)+y(3))*y(7)/(K_NO_Den+y(7)+y(7)*y(7)/K_I4NO)*K_OH4/(K_OH4+y(1))*y(4)/(K_NH4+y(4))*y(10)/(K_P+y(10))*y(12)/(K_ALK+y(12))*y(16);
    Rho(10,1) = u_H*n_G5*y(2)/(K_S5+y(2))*y(2)/(y(2)+y(3))*y(6)/(K_N2O_Den+y(6))*K_OH5/(K_OH5+y(1))*y(4)/(K_NH4+y(4))*y(10)/(K_P+y(10))*y(12)/(K_ALK+y(12))*y(16);
    Rho(11,1) = u_H*n_NO3_D*y(3)/(K_A+y(3))*y(3)/(y(2)+y(3))*K_O2/(K_O2+y(1))*y(9)/(K_NO3+y(9))*y(4)/(K_NH4+y(4))*y(10)/(K_P+y(10))*y(12)/(K_ALK+y(12))*y(16);
    Rho(12,1) = u_H*n_G3*y(3)/(K_S3+y(3))*y(3)/(y(2)+y(3))*y(8)/(K_NO2_Den+y(8))*K_OH3/(K_OH3+y(1))*y(4)/(K_NH4+y(4))*y(10)/(K_P+y(10))*y(12)/(K_ALK+y(12))*y(16);
    Rho(13,1) = u_H*n_G4*y(3)/(K_S4+y(3))*y(3)/(y(2)+y(3))*y(7)/(K_NO_Den+y(7)+y(7)*y(7)/K_I4NO)*K_OH4/(K_OH4+y(1))*y(4)/(K_NH4+y(4))*y(10)/(K_P+y(10))*y(12)/(K_ALK+y(12))*y(16);
    Rho(14,1) = u_H*n_G5*y(3)/(K_S5+y(3))*y(3)/(y(2)+y(3))*y(6)/(K_N2O_Den+y(6))*K_OH5/(K_OH5+y(1))*y(4)/(K_NH4+y(4))*y(10)/(K_P+y(10))*y(12)/(K_ALK+y(12))*y(16);
    Rho(15,1) = q_fe*K_O2/(K_O2+y(1))*K_NO2/(K_NO2+(y(9)+y(8)))*y(2)/(K_fe_H+y(2))*y(12)/(K_ALK+y(12))*y(16);
    Rho(16,1) = b_H*y(16);
% Part C: PAO
    Rho(17,1) = q_PHA*y(3)/(K_A_P+y(3))*y(12)/(K_ALK_P+y(12))*(y(18)/y(17))/(K_PP_P+(y(18)/y(17)))*y(17);
    Rho(18,1) = q_PP*y(1)/(K_O2_P+y(1))*y(10)/(K_P_P+y(10))*y(12)/(K_ALK_P+y(12))*(y(19)/y(17))/(K_PHA_P+(y(19)/y(17)))*(K_MAX_P-(y(18)/y(17)))/(K_IPP_P+K_MAX_P-(y(18)/y(17)))*y(17);
    Rho(19,1) = q_PP*n_NO3_P*y(9)/(K_NO3_P+y(9))*K_O2_P/(K_O2_P+y(1))*y(10)/(K_P_P+y(10))*y(12)/(K_ALK_P+y(12))*(y(19)/y(17))/(K_PHA_P+(y(19)/y(17)))*(K_MAX_P-(y(18)/y(17)))/(K_IPP_P+K_MAX_P-(y(18)/y(17)))*y(17);
    Rho(20,1) = q_PP*n_G3*y(8)/(K_NO2_Den+y(8))*K_OH3/(K_OH3+y(1))*y(10)/(K_P_P+y(10))*y(12)/(K_ALK_P+y(12))*(y(19)/y(17))/(K_PHA_P+(y(19)/y(17)))*(K_MAX_P-(y(18)/y(17)))/(K_IPP_P+K_MAX_P-(y(18)/y(17)))*y(17);
    Rho(21,1) = q_PP*n_G4*y(7)/(K_NO_Den+y(7)+y(7)*y(7)/K_I4NO)*K_OH4/(K_OH4+y(1))*y(10)/(K_P_P+y(10))*y(12)/(K_ALK_P+y(12))*(y(19)/y(17))/(K_PHA_P+(y(19)/y(17)))*(K_MAX_P-(y(18)/y(17)))/(K_IPP_P+K_MAX_P-(y(18)/y(17)))*y(17);
    Rho(22,1) = q_PP*n_G5*y(6)/(K_N2O_Den+y(6))*K_OH5/(K_OH5+y(1))*y(10)/(K_P_P+y(10))*y(12)/(K_ALK_P+y(12))*(y(19)/y(17))/(K_PHA_P+(y(19)/y(17)))*(K_MAX_P-(y(18)/y(17)))/(K_IPP_P+K_MAX_P-(y(18)/y(17)))*y(17);
    Rho(23,1) = u_PAO*y(1)/(K_O2_P+y(1))*y(10)/(K_P_P+y(10))*y(4)/(K_NH4+y(4))*y(12)/(K_ALK_P+y(12))*(y(19)/y(17))/(K_PHA_P+(y(19)/y(17)))*y(17);
    Rho(24,1) = u_PAO*n_NO3_P*y(9)/(K_NO3_P+y(9))*K_O2_P/(K_O2_P+y(1))*y(10)/(K_P_P+y(10))*y(4)/(K_NH4+y(4))*y(12)/(K_ALK_P+y(12))*(y(19)/y(17))/(K_PHA_P+(y(19)/y(17)))*y(17);
    Rho(25,1) = u_PAO*n_G3*y(8)/(K_NO2_Den+y(8))*K_OH3/(K_OH3+y(1))*y(10)/(K_P_P+y(10))*y(4)/(K_NH4+y(4))*y(12)/(K_ALK_P+y(12))*(y(19)/y(17))/(K_PHA_P+(y(19)/y(17)))*y(17);
    Rho(26,1) = u_PAO*n_G4*y(7)/(K_NO_Den+y(7)+y(7)*y(7)/K_I4NO)*K_OH4/(K_OH4+y(1))*y(10)/(K_P_P+y(10))*y(4)/(K_NH4+y(4))*y(12)/(K_ALK_P+y(12))*(y(19)/y(17))/(K_PHA_P+(y(19)/y(17)))*y(17);
    Rho(27,1) = u_PAO*n_G5*y(6)/(K_N2O_Den+y(6))*K_OH5/(K_OH5+y(1))*y(10)/(K_P_P+y(10))*y(4)/(K_NH4+y(4))*y(12)/(K_ALK_P+y(12))*(y(19)/y(17))/(K_PHA_P+(y(19)/y(17)))*y(17);
    Rho(28,1) = b_PAO*y(12)/(K_ALK_P+y(12))*y(17);
    Rho(29,1) = b_PP*y(12)/(K_ALK_P+y(12))*y(18);
    Rho(30,1) = b_PHA*y(12)/(K_ALK_P+y(12))*y(19);
% Part D: Denitrifying organisms (process 35 is at 20 degree, PH 7)
    Rho(31,1) = q_AOB_AMO*y(1)/(K_O2_AOB1+y(1))*y(4)/(K_NH4_AOB+y(4))*y(20);
    Rho(32,1) = u_AOB_HAO*y(1)/(K_O2_AOB2+y(1))*y(5)/(K_NH2OH_AOB+y(5))*y(4)/(y(4)+10^(-12))*y(10)/(K_P_AOB+y(10))*y(12)/(K_ALK_AOB+y(12))*y(20);
    Rho(33,1) = q_AOB_HAO*y(1)/(K_O2_AOB2+y(1))*y(7)/(K_NO_AOB_HAO+y(7))*y(20);
    Rho(34,1) = q_AOB_N2O_NN*y(5)/(K_NH2OH_AOB+y(5))*y(7)/(K_NO_AOB_NN+y(7))*y(20);
    Rho(35,1) = q_AOB_N2O_ND*y(5)/(K_NH2OH_AOB+y(5))*y(8)*0.000861/(K_HNO2_AOB+y(8)*0.000861)*y(1)/(K_O2_AOB_ND+(1-2*sqrt(K_O2_AOB_ND/K_I_O2_AOB))*y(1)+y(1)*y(1)/K_I_O2_AOB)*y(20);
    Rho(36,1) = u_NOB*y(1)/(K_O2_NOB+y(1))*y(8)/(K_NO2_NOB+y(8))*y(10)/(K_P_NOB+y(10))*y(12)/(K_ALK_NOB+y(12))*y(21);
    Rho(37,1) = b_AOB*y(20);
    Rho(38,1) = b_NOB*y(21);
% Part E: precipitation of P with Fe(OH)3
    Rho(39,1) = k_PRE*y(10)*y(23);
    Rho(40,1) = k_RED*y(12)/(K_ALK_PR+y(12))*y(24);
% ODE group
    dydt = Sto'*Rho;
    dydt(1) = 0; % DO set to constant 2 mg/L
end
