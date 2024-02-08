function dy = Model_coll(y,neuralOdeParameters)

% Collocation model function
% y is a Matrix! The |model| function, which defines the neural network used 
% to make predictions. simply output matrix from the input matrix of the network 
% for the purpose to reuse the network structure.

dummyT = 0;
dy =odeModel(dummyT, y, neuralOdeParameters);

end