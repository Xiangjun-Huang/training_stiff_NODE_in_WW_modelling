function X = Model_NODE(tspan,X0,neuralOdeParameters)
% NODE model function
% The |model| function, which defines the neural network used to make predictions, 
% is composed of a single neural ODE call. For each observation, this function 
% takes a vector of length |stateSize|, which is used as initial condition for 
% solving numerically the ODE with the function |odeModel|, which represents the 
% learnable right-hand side $f(t,y,\theta)$ of the ODE to be solved, as right 
% hand side and a vector of time points |tspan| defining the time at which the 
% numerical solution is output. The function uses the vector |tspan| for each 
% observation, regardless of the initial condition, since the learned system is 
% autonomous. That is, the |odeModel| function does not explicitly depend on time.

X = dlode45(@odeModel,tspan,X0,neuralOdeParameters,DataFormat="BC");

end