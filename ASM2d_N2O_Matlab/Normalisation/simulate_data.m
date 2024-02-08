%% ODE15s solver

numTimeSteps = 1000;
tSpan = [0 1/4]; % unit: day, 1/4 day = 6 hours
odeOptions = odeset(RelTol=1e-12, AbsTol=1e-12, NonNegative=1);
t = linspace(0, tSpan(2), numTimeSteps);

trueModel = @(t,y) ASM2d_N2O_Model(t, y, Stoi_Mat, Kinet_Vec);
[t,ySol] = ode15s(trueModel, t, y0, odeOptions);
% Calculate true y' (dySol) by ASM2d_N2O model

numR = size(ySol,1);
dySol = zeros(numR,numFeature);
for k = 1:numR
    dySol(k,:) = ASM2d_N2O_Model(t(k),ySol(k,:), Stoi_Mat, Kinet_Vec);
end

tHour = t*24; % time unit: from day to hour

% Add noise to the true y with noise sigma (normally varying from 10% ~ 0.5%) as simulated observation data

noiseSigma = 0.05;

ySolMean = mean(ySol);
ySolNoise = ySol +  randn(size(ySol,1),1) * (noiseSigma * ySolMean);