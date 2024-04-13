%% This function simulate the trajectory data for an IVP of ASM2d-N2O model
% The initial condition y0, size (1, num of components)
% time length of the simulation: tLen 
% number of discretised points for the trajectory: numTimeSteps

function [t, yTrue, dyTrue] = simulate_data(opt)
    arguments
        opt.y0 (1,:) double = [2,48.5,32.3,20,0,0,0,0,2.6,9,48.5,5,0,40.4,202,2500,250,70,100,200,100,189,50,220];
        opt.tLen (1,1) double = 1/4;
        opt.numTimeSteps (1,1) double = 1000;
    end
    
    %% Initialise the default ASM2d_N2O parameters
    if isfile("data/stoiMatrix.mat") && isfile("data/kinetVector.mat") 
        load data/stoiMatrix.mat stoiMatrix;
        load data/kinetVector.mat kinetVector;
    else
        [stoiMatrix, kinetVector, ~] = initialize_ASM2d_N2O_parameters();
    end

    %% ODE15s solver
    
    tSpan = [0 opt.tLen]; % unit: day, 1/4 day = 6 hours
    odeOptions = odeset(RelTol=1e-12, AbsTol=1e-12, NonNegative=1);
    t = linspace(0, tSpan(2), opt.numTimeSteps);
    
    trueModel = @(t,y) ASM2d_N2O_model(t, y, stoiMatrix, kinetVector);
    [t,yTrue] = ode15s(trueModel, t, opt.y0, odeOptions);
    
    % Calculate true y' (dyTrue) by ASM2d_N2O model
    
    dyTrue = zeros(size(yTrue));
    for k = 1:opt.numTimeSteps
        dyTrue(k,:) = ASM2d_N2O_model(t(k),yTrue(k,:), stoiMatrix, kinetVector);
    end
    
end

%% Mathematical ASM2d_N2O ODE model
% t, time, the unit depends on the unit of kinetic value, it is day here.  
% y denotes 24 components. StoiMatrix denotes stoichimetric matrix,
% KinetVector denotes kinetic value in a vector, must be in the defined order.

function dydt=ASM2d_N2O_model(t,y,stoiMatrix,kinetVector)
    % The order of all the values in the kinetic vector:
    % Part A: Hydrolysis processes
    % K_H, K_O2_H, K_x_H, n_NO3_H, n_NO2_H, K_NO3_H, K_NO2_H, n_fe_H
    % Part B: Heterotrophic biomass XH
    % u_H, K_O2, K_F, K_NH4, K_P, K_ALK, K_A, K_NO3, K_NO2, n_NO3_D, q_fe, K_fe_H, 
    % b_H, u_H_Den, 
    % 
    % n_G3, n_G4, n_G5, K_S3, K_S4, K_S5, K_NO2_Den, K_OH4, K_N2O_Den, K_OH3, K_NO_Den, 
    % K_OH5, K_I3NO, K_I4NO, K_I5NO,
    % Part C: PAO
    % q_PHA, K_A_P, K_ALK_P, q_PP, K_O2_P, K_P_P, K_PHA_P, K_MAX_P, K_PP_P, K_IPP_P, 
    % K_PO4_P, n_NO3_P, n_NO2_P,
    % 
    % K_NO3_P, K_NO2_P, u_PAO, b_PAO, b_PP, b_PHA
    % Part D: Nitrifying Organisms
    % u_AOB_HAO, q_AOB_AMO, K_O2_AOB1, K_NH4_AOB, K_O2_AOB2, K_NH2OH_AOB, q_AOB_HAO, 
    % K_NO_AOB_HAO, q_AOB_N2O_NN,
    % 
    % K_NO_AOB_NN, K_O2_AOB_ND, K_I_O2_AOB, K_HNO2_AOB, q_AOB_N2O_ND, K_ALK_AOB, 
    % K_P_AOB, u_NOB, K_O2_NOB, K_ALK_NOB,
    % 
    % K_NO2_NOB, K_P_NOB, b_AOB, b_NOB
    % Part E: Precipitation of P with Fe(OH)3
    % k_PRE, k_RED, K_ALK_PR
    
    Kine_cell = num2cell(kinetVector);
    
    [K_H, K_O2_H, K_x_H, n_NO3_H, n_NO2_H, K_NO3_H, K_NO2_H, n_fe_H,...
     u_H, K_O2, K_F, K_NH4, K_P, K_ALK, K_A, K_NO3, K_NO2, n_NO3_D, q_fe, K_fe_H, b_H, u_H_Den,... 
     n_G3, n_G4, n_G5, K_S3, K_S4, K_S5, K_NO2_Den, K_OH4, K_N2O_Den, K_OH3, K_NO_Den, K_OH5, K_I3NO, K_I4NO, K_I5NO,...
     q_PHA, K_A_P, K_ALK_P, q_PP, K_O2_P, K_P_P, K_PHA_P, K_MAX_P, K_PP_P, K_IPP_P, K_PO4_P, n_NO3_P, n_NO2_P,... 
     K_NO3_P, K_NO2_P, u_PAO, b_PAO, b_PP, b_PHA,...
     u_AOB_HAO, q_AOB_AMO, K_O2_AOB1, K_NH4_AOB, K_O2_AOB2, K_NH2OH_AOB, q_AOB_HAO, K_NO_AOB_HAO, q_AOB_N2O_NN,...
     K_NO_AOB_NN, K_O2_AOB_ND, K_I_O2_AOB, K_HNO2_AOB, q_AOB_N2O_ND, K_ALK_AOB, K_P_AOB, u_NOB, K_O2_NOB, K_ALK_NOB,...
     K_NO2_NOB, K_P_NOB, b_AOB, b_NOB,...
     k_PRE, k_RED, K_ALK_PR] = Kine_cell{:};

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
    dydt = stoiMatrix'*Rho;
    dydt(1) = 0; % DO set to constant 2 mg/L
end

%% initialize ASM2d-N2O model parameters with default values

function [Stoi_Mat, Kinet_Vec, Continuity_Check] = initialize_ASM2d_N2O_parameters()
% Default values of stoichiometric coefficients
% Refer to the paper for details.

Y_H = 0.625; % gCOD/gCOD Yield coefficient of Heterotrophic bacteria
Y_PHA = 0.2; % gCOD/gP PHA requirement for PP storage
Y_PAO = 0.625; % gCOD/gCOD Yield coefficient (biomass/PHA)
Y_PO4 = 0.4; % gP/gCOD PP requirement (PO43- release) per PHA stored
Y_AOB = 0.18; % gCOD/gN 
Y_NOB = 0.08; % gCOD/gN
f_SI = 0; % gCOD/gCOD
f_XI = 0.1; % gCOD/gCOD
n_G = 1; % dimensionless, Anoxic growth factor

i_NSI = 0.01; % gN/gCOD
i_NSF = 0.03; % gN/gCOD
i_NXI = 0.02; % gN/gCOD
i_NXS = 0.04; % gN/gCOD
i_NBM = 0.07; % gN/gCOD

i_PSI = 0; % gP/gCOD
i_PSF = 0.01; % gP/gCOD
i_PXI = 0.01; % gP/gCOD
i_PXS = 0.01; % gP/gCOD
i_PBM = 0.02; % gP/gCOD

i_TSSXI = 0.75; % gTSS/gCOD
i_TSSXS = 0.75; % gTSS/gCOD
i_TSSBM = 0.9; % gTSS/gCOD
% Initialise composition conversion factor matrix
% Dimension of composition converstion factor matrix is 24 components x 5 composition 
% factors (COD N P Charge TSS), refer to the paper for details.

Comp_Mat = zeros(24,5);
% C1: S_O2

Comp_Mat(1,1) = -1;
% C2: S_F

Comp_Mat(2,1) = 1;
Comp_Mat(2,2) = i_NSF;
Comp_Mat(2,3) = i_PSF;
% C3: S_A

Comp_Mat(3,1) = 1;
Comp_Mat(3,4) = -1/64;
% C4: S_NH4

Comp_Mat(4,2) = 1;
Comp_Mat(4,4) = 1/14;
% C5: S_NH2OH

Comp_Mat(5,1) = -8/7;
Comp_Mat(5,2) = 1;
% C6: S_N2O

Comp_Mat(6,1) = -16/7;
Comp_Mat(6,2) = 1;
% C7: S_N2O

Comp_Mat(7,1) = -20/7;
Comp_Mat(7,2) = 1;
% C8: S_NO2

Comp_Mat(8,1) = -24/7;
Comp_Mat(8,2) = 1;
Comp_Mat(8,4) = -1/14;
% C9: S_NO3

Comp_Mat(9,1) = -32/7;
Comp_Mat(9,2) = 1;
Comp_Mat(9,4) = -1/14;
% C10: S_PO4

Comp_Mat(10,3) = 1;
Comp_Mat(10,4) = -1.5/31;
% C11: S_I

Comp_Mat(11,1) = 1;
Comp_Mat(11,2) = i_NSI;
Comp_Mat(11,3) = i_PSI;
% C12: S_ALK

Comp_Mat(12,4) = -1;
% C13: S_N2

Comp_Mat(13,1) = -24/14;
Comp_Mat(13,2) = 1;
% C14: X_I

Comp_Mat(14,1) = 1;
Comp_Mat(14,2) = i_NXI;
Comp_Mat(14,3) = i_PXI;
Comp_Mat(14,5) = i_TSSXI;
% C15: X_S

Comp_Mat(15,1) = 1;
Comp_Mat(15,2) = i_NXS;
Comp_Mat(15,3) = i_PXS;
Comp_Mat(15,5) = i_TSSXS;
% C16: X_H

Comp_Mat(16,1) = 1;
Comp_Mat(16,2) = i_NBM;
Comp_Mat(16,3) = i_PBM;
Comp_Mat(16,5) = i_TSSBM;
% C17: X_PAO

Comp_Mat(17,1) = 1;
Comp_Mat(17,2) = i_NBM;
Comp_Mat(17,3) = i_PBM;
Comp_Mat(17,5) = i_TSSBM;
% C18: X_PP

Comp_Mat(18,3) = 1;
Comp_Mat(18,4) = -1/31;
Comp_Mat(18,5) = 3.23;
% C19: X_PHA

Comp_Mat(19,1) = 1;
Comp_Mat(19,5) = 0.6;
% C20: X_AOB

Comp_Mat(20,1) = 1;
Comp_Mat(20,2) = i_NBM;
Comp_Mat(20,3) = i_PBM;
Comp_Mat(20,5) = i_TSSBM;
% C21: X_NOB

Comp_Mat(21,1) = 1;
Comp_Mat(21,2) = i_NBM;
Comp_Mat(21,3) = i_PBM;
Comp_Mat(21,5) = i_TSSBM;
% C22: X_TSS

Comp_Mat(22,5) = -1;
% C23: X_MeOH

Comp_Mat(23,5) = 1;
% C24: X_MeP

Comp_Mat(24,3) = 30.974/150.816;
Comp_Mat(24,5) = 1;
%% Initialise stoichiometric matrix
% Dimension of stoichiometric matrix is 40 reactions x 24 components, refer 
% to the paper for details.
% 
% Component 12 S_ALK and 22 X_TSS are dependent on other variables, so put at 
% the end.

Stoi_Mat = zeros(40,24);
% Process 1: aerobic hydrolysis 

Stoi_Mat(1,2) = 1-f_SI;
Stoi_Mat(1,4) = i_NXS-(1-f_SI)*i_NSF;
Stoi_Mat(1,10) = i_PXS-(1-f_SI)*i_PSF;
Stoi_Mat(1,11) = f_SI;
Stoi_Mat(1,15) = -1;
Stoi_Mat(1,12) = Stoi_Mat(1,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(1,22) = Stoi_Mat(1,14:21)*Comp_Mat(14:21,5);
% Process 2: Anoxic hydrolysis (1st step: NO3-)

Stoi_Mat(2,2) = 1-f_SI;
Stoi_Mat(2,4) = i_NXS-(1-f_SI)*i_NSF;
Stoi_Mat(2,10) = i_PXS-(1-f_SI)*i_PSF;
Stoi_Mat(2,11) = f_SI;
Stoi_Mat(2,15) = -1;
Stoi_Mat(2,12) = Stoi_Mat(2,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(2,22) = Stoi_Mat(2,14:21)*Comp_Mat(14:21,5);
% Process 3: Anoxic hydrolysis (2nd step: NO2-)

Stoi_Mat(3,2) = 1-f_SI;
Stoi_Mat(3,4) = i_NXS-(1-f_SI)*i_NSF;
Stoi_Mat(3,10) = i_PXS-(1-f_SI)*i_PSF;
Stoi_Mat(3,11) = f_SI;
Stoi_Mat(3,15) = -1;
Stoi_Mat(3,12) = Stoi_Mat(3,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(3,22) = Stoi_Mat(3,14:21)*Comp_Mat(14:21,5);
% Process 4: Anaerobic hydrolysis 

Stoi_Mat(4,2) = 1-f_SI;
Stoi_Mat(4,4) = i_NXS-(1-f_SI)*i_NSF;
Stoi_Mat(4,10) = i_PXS-(1-f_SI)*i_PSF;
Stoi_Mat(4,11) = f_SI;
Stoi_Mat(4,15) = -1;
Stoi_Mat(4,12) = Stoi_Mat(4,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(4,22) = Stoi_Mat(4,14:21)*Comp_Mat(14:21,5);
% Process 5: Aerobic gorwth on S_F 

Stoi_Mat(5,1) = 1-(1/Y_H);
Stoi_Mat(5,2) = -1/Y_H;
Stoi_Mat(5,4) = -i_NBM+(1/Y_H)*i_NSF;
Stoi_Mat(5,10) = -i_PBM+(1/Y_H)*i_PSF;
Stoi_Mat(5,16) = 1;
Stoi_Mat(5,12) = Stoi_Mat(5,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(5,22) = Stoi_Mat(5,14:21)*Comp_Mat(14:21,5);
% Process 6: Aerobic gorwth on S_A 

Stoi_Mat(6,1) = 1-(1/Y_H);
Stoi_Mat(6,3) = -1/Y_H;
Stoi_Mat(6,4) = -i_NBM;
Stoi_Mat(6,10) = -i_PBM;
Stoi_Mat(6,16) = 1;
Stoi_Mat(6,12) = Stoi_Mat(6,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(6,22) = Stoi_Mat(6,14:21)*Comp_Mat(14:21,5);
% Process 7: Anoxic gorwth of heterotrophs on S_F (NO3- > NO2-) 

Stoi_Mat(7,2) = -1/(Y_H*n_G);
Stoi_Mat(7,4) = -i_NBM+(1/Y_H)*i_NSF;
Stoi_Mat(7,8) = (1-Y_H*n_G)/((8/7)*Y_H*n_G);
Stoi_Mat(7,9) = -(1-Y_H*n_G)/((8/7)*Y_H*n_G);
Stoi_Mat(7,10) = -i_PBM+(1/Y_H)*i_PSF;
Stoi_Mat(7,16) = 1;
Stoi_Mat(7,12) = Stoi_Mat(7,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(7,22) = Stoi_Mat(7,14:21)*Comp_Mat(14:21,5);
% Process 8: Anoxic gorwth of heterotrophs on S_F (NO2- > NO) 

Stoi_Mat(8,2) = -1/(Y_H*n_G);
Stoi_Mat(8,4) = -i_NBM+(1/Y_H)*i_NSF;
Stoi_Mat(8,7) = (1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(8,8) = -(1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(8,10) = -i_PBM+(1/Y_H)*i_PSF;
Stoi_Mat(8,16) = 1;
Stoi_Mat(8,12) = Stoi_Mat(8,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(8,22) = Stoi_Mat(8,14:21)*Comp_Mat(14:21,5);
% Process 9: Anoxic gorwth of heterotrophs on S_F (NO- > N2O)

Stoi_Mat(9,2) = -1/(Y_H*n_G);
Stoi_Mat(9,4) = -i_NBM+(1/Y_H)*i_NSF;
Stoi_Mat(9,6) = (1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(9,7) = -(1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(9,10) = -i_PBM+(1/Y_H)*i_PSF;
Stoi_Mat(9,16) = 1;
Stoi_Mat(9,12) = Stoi_Mat(9,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(9,22) = Stoi_Mat(9,14:21)*Comp_Mat(14:21,5);
% Process 10: Anoxic gorwth of heterotrophs on S_F (N2O- > N2)

Stoi_Mat(10,2) = -1/(Y_H*n_G);
Stoi_Mat(10,4) = -i_NBM+(1/Y_H)*i_NSF;
Stoi_Mat(10,6) = -(1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(10,10) = -i_PBM+(1/Y_H)*i_PSF;
Stoi_Mat(10,13) = (1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(10,16) = 1;
Stoi_Mat(10,12) = Stoi_Mat(10,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(10,22) = Stoi_Mat(10,14:21)*Comp_Mat(14:21,5);
% Process 11: Anoxic gorwth of heterotrophs on S_A (NO3- > NO2-) 

Stoi_Mat(11,3) = -1/(Y_H*n_G);
Stoi_Mat(11,4) = -i_NBM;
Stoi_Mat(11,8) = (1-Y_H*n_G)/((8/7)*Y_H*n_G);
Stoi_Mat(11,9) = -(1-Y_H*n_G)/((8/7)*Y_H*n_G);
Stoi_Mat(11,10) = -i_PBM;
Stoi_Mat(11,16) = 1;
Stoi_Mat(11,12) = Stoi_Mat(11,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(11,22) = Stoi_Mat(11,14:21)*Comp_Mat(14:21,5);
% Process 12: Anoxic gorwth of heterotrophs on S_A (NO2- > NO) 

Stoi_Mat(12,3) = -1/(Y_H*n_G);
Stoi_Mat(12,4) = -i_NBM;
Stoi_Mat(12,7) = (1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(12,8) = -(1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(12,10) = -i_PBM;
Stoi_Mat(12,16) = 1;
Stoi_Mat(12,12) = Stoi_Mat(12,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(12,22) = Stoi_Mat(12,14:21)*Comp_Mat(14:21,5);
% Process 13: Anoxic gorwth of heterotrophs on S_A (NO- > N2O)

Stoi_Mat(13,3) = -1/(Y_H*n_G);
Stoi_Mat(13,4) = -i_NBM;
Stoi_Mat(13,6) = (1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(13,7) = -(1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(13,10) = -i_PBM;
Stoi_Mat(13,16) = 1;
Stoi_Mat(13,12) = Stoi_Mat(13,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(13,22) = Stoi_Mat(13,14:21)*Comp_Mat(14:21,5);
% Process 14: Anoxic gorwth of heterotrophs on S_A (N2O- > N2) 

Stoi_Mat(14,3) = -1/(Y_H*n_G);
Stoi_Mat(14,4) = -i_NBM;
Stoi_Mat(14,6) = -(1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(14,10) = -i_PBM;
Stoi_Mat(14,13) = (1-Y_H*n_G)/((4/7)*Y_H*n_G);
Stoi_Mat(14,16) = 1;
Stoi_Mat(14,12) = Stoi_Mat(14,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(14,22) = Stoi_Mat(14,14:21)*Comp_Mat(14:21,5);
% Process 15: Fermentation 

Stoi_Mat(15,2) = -1;
Stoi_Mat(15,3) = 1;
Stoi_Mat(15,4) = i_NSF;
Stoi_Mat(15,10) = i_PSF;
Stoi_Mat(15,12) = Stoi_Mat(15,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(15,22) = Stoi_Mat(15,14:21)*Comp_Mat(14:21,5);
% Process 16: Lysis

Stoi_Mat(16,4) = i_NBM-i_NXI*f_XI-(1-f_XI)*i_NXS;
Stoi_Mat(16,10) = i_PBM-i_PXI*f_XI-(1-f_XI)*i_PXS;
Stoi_Mat(16,14) = f_XI;
Stoi_Mat(16,15) = 1-f_XI;
Stoi_Mat(16,16) = -1;
Stoi_Mat(16,12) = Stoi_Mat(16,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(16,22) = Stoi_Mat(16,14:21)*Comp_Mat(14:21,5);
% Process 17: Storage of X_PHA 

Stoi_Mat(17,3) = -1;
Stoi_Mat(17,10) = Y_PO4;
Stoi_Mat(17,18) = -Y_PO4;
Stoi_Mat(17,19) = 1;
Stoi_Mat(17,12) = Stoi_Mat(17,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(17,22) = Stoi_Mat(17,14:21)*Comp_Mat(14:21,5);
% Process 18: Aerobic storage of X_PP 

Stoi_Mat(18,1) = -Y_PHA;
Stoi_Mat(18,10) = -1;
Stoi_Mat(18,18) = 1;
Stoi_Mat(18,19) = -Y_PHA;
Stoi_Mat(18,12) = Stoi_Mat(18,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(18,22) = Stoi_Mat(18,14:21)*Comp_Mat(14:21,5);
% Process 19: Anoxic storage of X_PP (NO3- > NO2-)

Stoi_Mat(19,8) = Y_PHA/(8/7);
Stoi_Mat(19,9) = -Y_PHA/(8/7);
Stoi_Mat(19,10) = -1;
Stoi_Mat(19,18) = 1;
Stoi_Mat(19,19) = -Y_PHA;
Stoi_Mat(19,12) = Stoi_Mat(19,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(19,22) = Stoi_Mat(19,14:21)*Comp_Mat(14:21,5);
% Process 20: Anoxic storage of X_PP (NO2- > NO)

Stoi_Mat(20,7) = Y_PHA/(4/7);
Stoi_Mat(20,8) = -Y_PHA/(4/7);
Stoi_Mat(20,10) = -1;
Stoi_Mat(20,18) = 1;
Stoi_Mat(20,19) = -Y_PHA;
Stoi_Mat(20,12) = Stoi_Mat(20,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(20,22) = Stoi_Mat(20,14:21)*Comp_Mat(14:21,5);
% Process 21: Anoxic storage of X_PP (NO > N2O)

Stoi_Mat(21,6) = Y_PHA/(4/7);
Stoi_Mat(21,7) = -Y_PHA/(4/7);
Stoi_Mat(21,10) = -1;
Stoi_Mat(21,18) = 1;
Stoi_Mat(21,19) = -Y_PHA;
Stoi_Mat(21,12) = Stoi_Mat(21,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(21,22) = Stoi_Mat(21,14:21)*Comp_Mat(14:21,5);
% Process 22: Anoxic storage of X_PP (N2O > N2)

Stoi_Mat(22,6) = -Y_PHA/(4/7);
Stoi_Mat(22,10) = -1;
Stoi_Mat(22,13) = Y_PHA/(4/7);
Stoi_Mat(22,18) = 1;
Stoi_Mat(22,19) = -Y_PHA;
Stoi_Mat(22,12) = Stoi_Mat(22,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(22,22) = Stoi_Mat(22,14:21)*Comp_Mat(14:21,5);
% Process 23: Aerobic growth of X_PAO

Stoi_Mat(23,1) = 1-(1/Y_PAO);
Stoi_Mat(23,4) = -i_NBM;
Stoi_Mat(23,10) = -i_PBM;
Stoi_Mat(23,17) = 1;
Stoi_Mat(23,19) = -1/Y_PAO;
Stoi_Mat(23,12) = Stoi_Mat(23,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(23,22) = Stoi_Mat(23,14:21)*Comp_Mat(14:21,5);
% Process 24: Anoxic growth of X_PAO (NO3- > NO2-)

Stoi_Mat(24,4) = -i_NBM;
Stoi_Mat(24,8) = (1-Y_PAO*n_G)/((8/7)*Y_PAO*n_G);
Stoi_Mat(24,9) = -(1-Y_PAO*n_G)/((8/7)*Y_PAO*n_G);
Stoi_Mat(24,10) = -i_PBM;
Stoi_Mat(24,17) = 1;
Stoi_Mat(24,19) = -1/Y_PAO;
Stoi_Mat(24,12) = Stoi_Mat(24,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(24,22) = Stoi_Mat(24,14:21)*Comp_Mat(14:21,5);
% Process 25: Anoxic growth of X_PAO (NO2- > NO)

Stoi_Mat(25,4) = -i_NBM;
Stoi_Mat(25,7) = (1-Y_PAO*n_G)/((4/7)*Y_PAO*n_G);
Stoi_Mat(25,8) = -(1-Y_PAO*n_G)/((4/7)*Y_PAO*n_G);
Stoi_Mat(25,10) = -i_PBM;
Stoi_Mat(25,17) = 1;
Stoi_Mat(25,19) = -1/Y_PAO;
Stoi_Mat(25,12) = Stoi_Mat(25,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(25,22) = Stoi_Mat(25,14:21)*Comp_Mat(14:21,5);
% Process 26: Anoxic growth of X_PAO (NO- > N2O)

Stoi_Mat(26,4) = -i_NBM;
Stoi_Mat(26,6) = (1-Y_PAO*n_G)/((4/7)*Y_PAO*n_G);
Stoi_Mat(26,7) = -(1-Y_PAO*n_G)/((4/7)*Y_PAO*n_G);
Stoi_Mat(26,10) = -i_PBM;
Stoi_Mat(26,17) = 1;
Stoi_Mat(26,19) = -1/Y_PAO;
Stoi_Mat(26,12) = Stoi_Mat(26,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(26,22) = Stoi_Mat(26,14:21)*Comp_Mat(14:21,5);
% Process 27: Anoxic growth of X_PAO (N2O- > N2)

Stoi_Mat(27,4) = -i_NBM;
Stoi_Mat(27,6) = -(1-Y_PAO*n_G)/((4/7)*Y_PAO*n_G);
Stoi_Mat(27,10) = -i_PBM;
Stoi_Mat(27,13) = (1-Y_PAO*n_G)/((4/7)*Y_PAO*n_G);
Stoi_Mat(27,17) = 1;
Stoi_Mat(27,19) = -1/Y_PAO;
Stoi_Mat(27,12) = Stoi_Mat(27,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(27,22) = Stoi_Mat(27,14:21)*Comp_Mat(14:21,5);
% Process 28: lysis of X_PAO

Stoi_Mat(28,4) = i_NBM-i_NXI*f_XI-(1-f_XI)*i_NXS;
Stoi_Mat(28,10) = i_PBM-i_PXI*f_XI-(1-f_XI)*i_PXS;
Stoi_Mat(28,14) = f_XI;
Stoi_Mat(28,15) = 1-f_XI;
Stoi_Mat(28,17) = -1;
Stoi_Mat(28,12) = Stoi_Mat(28,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(28,22) = Stoi_Mat(28,14:21)*Comp_Mat(14:21,5);
% Process 29: lysis of X_PP

Stoi_Mat(29,10) = 1;
Stoi_Mat(29,18) = -1;
Stoi_Mat(29,12) = Stoi_Mat(29,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(29,22) = Stoi_Mat(29,14:21)*Comp_Mat(14:21,5);
% Process 30: lysis of X_PHA

Stoi_Mat(30,3) = 1;
Stoi_Mat(30,19) = -1;
Stoi_Mat(30,12) = Stoi_Mat(30,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(30,22) = Stoi_Mat(30,14:21)*Comp_Mat(14:21,5);
% Process 31: (1) NH3 oxidation to NH2OH with oxygen consumption

Stoi_Mat(31,1) = -8/7;
Stoi_Mat(31,4) = -1;
Stoi_Mat(31,5) = 1;
Stoi_Mat(31,12) = Stoi_Mat(31,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(31,22) = Stoi_Mat(31,14:21)*Comp_Mat(14:21,5);
% Process 32: (2) NH2OH to NO coupled with oxygen reduction (X_AOB growth)

Stoi_Mat(32,1) = -(12/7-Y_AOB)/Y_AOB;
Stoi_Mat(32,4) = -i_NBM;
Stoi_Mat(32,5) = -1/Y_AOB;
Stoi_Mat(32,7) = 1/Y_AOB;
Stoi_Mat(32,10) = -i_PBM;
Stoi_Mat(32,20) = 1;
Stoi_Mat(32,12) = Stoi_Mat(32,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(32,22) = Stoi_Mat(32,14:21)*Comp_Mat(14:21,5);
% Process 33: (3) NO oxidation to NO2- coupled with oxygen reduction

Stoi_Mat(33,1) = -4/7;
Stoi_Mat(33,7) = -1;
Stoi_Mat(33,8) = 1;
Stoi_Mat(33,12) = Stoi_Mat(33,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(33,22) = Stoi_Mat(33,14:21)*Comp_Mat(14:21,5);
% Process 34: (4) NO to N2O coupled with NH2OH to NO2- (N2O from NN pathway)

Stoi_Mat(34,5) = -1;
Stoi_Mat(34,6) = 4;
Stoi_Mat(34,7) = -4;
Stoi_Mat(34,8) = 1;
Stoi_Mat(34,12) = Stoi_Mat(34,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(34,22) = Stoi_Mat(34,14:21)*Comp_Mat(14:21,5);
% Process 35: (5) HNO2 to N2O coupled with NH2OH to NO2- (N2O from ND pathway)

Stoi_Mat(35,5) = -1;
Stoi_Mat(35,6) = 2;
Stoi_Mat(35,8) = -1;
Stoi_Mat(35,12) = Stoi_Mat(35,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(35,22) = Stoi_Mat(35,14:21)*Comp_Mat(14:21,5);
% Process 36: Aerobic growth of X_NOB

Stoi_Mat(36,1) = -((8/7)-Y_NOB)/Y_NOB;
Stoi_Mat(36,4) = -i_NBM;
Stoi_Mat(36,8) = -1/Y_NOB;
Stoi_Mat(36,9) = 1/Y_NOB;
Stoi_Mat(36,10) = -i_PBM;
Stoi_Mat(36,21) = 1;
Stoi_Mat(36,12) = Stoi_Mat(36,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(36,22) = Stoi_Mat(36,14:21)*Comp_Mat(14:21,5);
% Process 37: Lysis of X_AOB

Stoi_Mat(37,4) = i_NBM-i_NXI*f_XI-(1-f_XI)*i_NXS;
Stoi_Mat(37,10) = i_PBM-i_PXI*f_XI-(1-f_XI)*i_PXS;
Stoi_Mat(37,14) = f_XI;
Stoi_Mat(37,15) = 1-f_XI;
Stoi_Mat(37,20) = -1;
Stoi_Mat(37,12) = Stoi_Mat(37,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(37,22) = Stoi_Mat(37,14:21)*Comp_Mat(14:21,5);
% Process 38: Lysis of X_NOB

Stoi_Mat(38,4) = i_NBM-i_NXI*f_XI-(1-f_XI)*i_NXS;
Stoi_Mat(38,10) = i_PBM-i_PXI*f_XI-(1-f_XI)*i_PXS;
Stoi_Mat(38,14) = f_XI;
Stoi_Mat(38,15) = 1-f_XI;
Stoi_Mat(38,21) = -1;
Stoi_Mat(38,12) = Stoi_Mat(38,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(38,22) = Stoi_Mat(38,14:21)*Comp_Mat(14:21,5);
% Process 39: Precipitation

Stoi_Mat(39,10) = -1;
Stoi_Mat(39,23) = -106.867/30.974;
Stoi_Mat(39,24) = 150.816/30.974;
Stoi_Mat(39,12) = Stoi_Mat(39,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(39,22) = Stoi_Mat(39,23)+Stoi_Mat(39,24); % Formula is changed from previous!!
% Process 40: Redissolution

Stoi_Mat(40,10) = 1;
Stoi_Mat(40,23) = 106.867/30.974;
Stoi_Mat(40,24) = -150.816/30.974;
Stoi_Mat(40,12) = Stoi_Mat(40,[3 4 8 9 10 18])*Comp_Mat([3 4 8 9 10 18],4);
Stoi_Mat(40,22) = Stoi_Mat(40,23)+Stoi_Mat(40,24); % Formula is changed from previous!!
% Default value of kinetic cofficients
% Refer to the paper for details. Where the value is temperature or pH dependant, 
% the condition is set at temperature: 20 degree, and pH: 7.
% Part A: Hydrolysis processes

K_H = 3; % /d
K_O2_H = 0.2; % gO2/m3
K_x_H = 0.1; % gXs/gXH
n_NO3_H = 0.6; 
n_NO2_H = 0.6;
K_NO3_H = 0.5; % gN/m3
K_NO2_H = 0.5; % gN/m3
n_fe_H = 0.4; 
% Part B: Heterotrophic biomass XH

u_H = 6; % gXs/(gXH.d)
K_O2 = 0.2; % gO2/m3
K_F = 4; % gCOD/m3
K_NH4 = 0.05; % gN/m3
K_P = 0.01; % gP/m3
K_ALK = 0.1; % mole HCO3-/m3
K_A = 4; % gCOD/m3
K_NO3 = 0.5; % gN/m3
K_NO2 = 0.5; % gN/m3
n_NO3_D = 0.8; 
q_fe = 3; % gSF/(gXH.d)
K_fe_H = 4; % gCOD/m3
b_H = 0.4; % /d
u_H_Den = 6.25; % /d
n_G3 = 0.16;
n_G4 = 0.35;
n_G5 = 0.35;
K_S3 = 20; % gCOD/m3
K_S4 = 20; % gCOD/m3
K_S5 = 40; % gCOD/m3
K_NO2_Den = 0.2; % gN/m3
K_OH4 = 0.1; % gO2/m3
K_N2O_Den = 0.05; % gN/m3
K_OH3 = 0.1; % gO3/m3
K_NO_Den = 0.05; % gN/m3
K_OH5 = 0.1; % gO2/m3
K_I3NO = 0.5; % gN/m3
K_I4NO = 0.3; % gN/m3
K_I5NO = 0.075; % gN/m3
% Part C: PAO

q_PHA = 3; % gXPHA/(gXPAO.d)
K_A_P = 4; % gCOD/m3
K_ALK_P = 0.1; % mole HCO3-/m3
q_PP = 1.5; % g XPP/(gXPAO.d)
K_O2_P = 0.2; % gO2/m3
K_P_P = 0.2; % gP/m3
K_PHA_P = 0.01; % gXPHA/gXPAO
K_MAX_P = 0.34; % gXPP/gXPAO
K_PP_P = 0.01; % gXPP/gXPAO
K_IPP_P = 0.02; % gXPP/gXPAO
K_PO4_P = 0.01; % gP/m3
n_NO3_P = 0.6; 
n_NO2_P = 0.6;
K_NO3_P = 0.5; % gN/m3
K_NO2_P = 0.5; % gN/m3
u_PAO = 1; % /d
b_PAO = 0.2; % /d
b_PP = 0.2; % /d
b_PHA = 0.2; % /d
% Part D: Nitrifying Organisms

u_AOB_HAO = 0.78; % /d
q_AOB_AMO = 5.2008; % gN/(gCOD.d)
K_O2_AOB1 = 1; % gO2/m3
K_NH4_AOB = 0.2; % gN/m3
K_O2_AOB2 = 0.6; % gO2/m3
K_NH2OH_AOB = 0.9; % gN/m3
q_AOB_HAO = 5.2008; % gN/(gCOD.d)
K_NO_AOB_HAO = 0.0003; % gN/m3
q_AOB_N2O_NN = 0.0078; % gN/(gCOD.d)
K_NO_AOB_NN = 0.008; % gN/m3
K_O2_AOB_ND = 0.5; % gO2/m3
K_I_O2_AOB = 0.8; % gO2/m3
K_HNO2_AOB = 0.004; % gN/m3
q_AOB_N2O_ND = 1.3008; % gN/(gCOD.d)
K_ALK_AOB = 0.1; % mole HCO3-/m3
K_P_AOB = 0.01; % gP/m3
u_NOB = 0.78; % /d first tested value
K_O2_NOB = 1.2; % gO2/m3 first tested value
K_ALK_NOB = 0.1; % mole HCO3-/m3
K_NO2_NOB = 0.5; % gN/m3
K_P_NOB = 0.01; % gP/m3
b_AOB = 0.096; % /d 
b_NOB = 0.096; % /d first tested value
% Part E: Precipitation of P with Fe(OH)3

k_PRE = 1; % m3/(g FeOH3.d) (lowcase k, not Capital K!!)
k_RED = 0.6; % /d (lowcase k, not Capital K!!!)
K_ALK_PR = 0.5; % mole HCO3-/m3
% Assemble all kinetic parameters into a vector, to pass to the solver

Kinet_Vec = [K_H, K_O2_H, K_x_H, n_NO3_H, n_NO2_H, K_NO3_H, K_NO2_H, n_fe_H,...
 u_H, K_O2, K_F, K_NH4, K_P, K_ALK, K_A, K_NO3, K_NO2, n_NO3_D, q_fe, K_fe_H, b_H, u_H_Den,... 
 n_G3, n_G4, n_G5, K_S3, K_S4, K_S5, K_NO2_Den, K_OH4, K_N2O_Den, K_OH3, K_NO_Den, K_OH5, K_I3NO, K_I4NO, K_I5NO,...
 q_PHA, K_A_P, K_ALK_P, q_PP, K_O2_P, K_P_P, K_PHA_P, K_MAX_P, K_PP_P, K_IPP_P, K_PO4_P, n_NO3_P, n_NO2_P,... 
 K_NO3_P, K_NO2_P, u_PAO, b_PAO, b_PP, b_PHA,...
 u_AOB_HAO, q_AOB_AMO, K_O2_AOB1, K_NH4_AOB, K_O2_AOB2, K_NH2OH_AOB, q_AOB_HAO, K_NO_AOB_HAO, q_AOB_N2O_NN,...
 K_NO_AOB_NN, K_O2_AOB_ND, K_I_O2_AOB, K_HNO2_AOB, q_AOB_N2O_ND, K_ALK_AOB, K_P_AOB, u_NOB, K_O2_NOB, K_ALK_NOB,...
 K_NO2_NOB, K_P_NOB, b_AOB, b_NOB,...
 k_PRE, k_RED, K_ALK_PR];
% Continuity check
% To conform the conservation law, the absolute result of Stoichiometric Matrix 
% (40 x 24) multiplying Composition Matrix (24 x 5) must be zero or less than 
% 1E-15.

Con_Mat = Stoi_Mat * Comp_Mat;
if ~isempty(Con_Mat(Con_Mat>1E-15))
    Continuity_Check = false;
    error("Continuity Check failed! Please check and correct your parameters input!");
else
    Continuity_Check = true;
end


% function end
end

