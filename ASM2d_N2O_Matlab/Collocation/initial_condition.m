%% Initial condition 
% Define all the 24 components value of initial condition. Some parameters are 
% defined in the paper, such as, total COD: 420 g/m3, NH4: 20 g/m3, NO3: 2.6 g/m3, 
% PO4: 9 g/m3, TSS: 189 g/m3. The typical value of other parameters can be find 
% from IWA AMD2d official report. Remember the unit of alkalinity is mole HCO3-/m3.
% 
% Initial condition and default values of coefficients can also be imported 
% from a external file, such as from Excel, or other file format.
% Define 24 components as a 24 x 1 column vector

numFeature = 24;
y0 = [2,48.5,32.3,20,0,0,0,0,2.6,9,48.5,5,0,40.4,202,2500,250,70,100,200,100,189,50,220];
% Define the name and unit of 24 components

nameFeature = ["S_{O_2}","S_F","S_A","S_{NH_4}","S_{NH_2OH}","S_{N_2O}",...
                "S_{NO}","S_{NO_2}","S_{NO_3}","S_{PO_4}","S_I","S_{ALK}",...
                "N_2","X_I","X_S","X_H","X_{PAO}","X_{PP}",...
                "X_{PHA}","X_{AOB}","X_{NOB}","X_{TSS}","X_{MeOH}","X_{Mep}"];
nameUnit = ["g(O_2)/m^3","g(COD)/m^3","g(COD)/m^3","g(N)/m^3","g(N)/m^3","g(N)/m^3",...
            "g(N)/m^3","g(N)/m^3","g(N)/m^3","g(P)/m^3","g(COD)/m^3","mol(HCO_3^-)/m^3",...
            "g(N)/m^3","g(COD)/m^3","g(COD)/m^3","g(COD)/m^3","g(COD)/m^3","g(P)/m^3",...
            "g(COD)/m^3","g(COD)/m^3","g(COD)/m^3","g(TSS)/m^3","g(TSS)/m^3","g(TSS)/m^3"];