Solving stiffness issue in training of neural ordinary differential equations for data-driven wastewater process modelling

The folder ASM1_python contains code in python simulating extended ASM1 model with neural ordinariy differential equations. The traing uses package of torchdiffeq and follows method domenstrated in paper "neural ordinary differential equations" by Ricky TQ Chen, with the z-score standardisation, max-min normalisation proposed by myself and equation scaling raised by Kim in his paper "stiff neural ordinary differential equations". You can "comment" or "uncomment" the lines in the code to choose the method to test.

The folder ASM2d_N2O_matlab contains code in matlab simulating ASM2d_N2O model with two methods: 1) normalisation in direct NODE training 2) collocation training followed direct NODE training with normalisation. Please click the "main" file to run.
