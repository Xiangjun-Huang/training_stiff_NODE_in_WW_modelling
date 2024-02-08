function dy = odeModel(~,y,theta)
% ODE Model
% y is a matrix! The |odeModel| function is the learnable right-hand side used 
% in the call to deep learning for collocation method and |dlode45 for NODE method|. 

global yMean yStd dyMean dyStd

idx = (yStd ~= 0);
y(:, idx) = (y(:, idx) - yMean(idx))./yStd(idx);

y = gelu(y*theta.fc1.Weights + theta.fc1.Bias);
y = gelu(y*theta.fc2.Weights + theta.fc2.Bias);
y = gelu(y*theta.fc3.Weights + theta.fc3.Bias);
y = y*theta.fc4.Weights + theta.fc4.Bias;

dy = y .* dyStd + dyMean;

end