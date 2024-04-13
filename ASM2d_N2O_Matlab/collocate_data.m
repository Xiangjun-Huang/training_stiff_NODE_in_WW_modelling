%% Define collocate_data function
% Xiangjun Huang, Brunel University London, Apr 2024

% y, y' = collocate_data(data,tpoints,kernel="EpanechnikovKernel"), computes 
% a non-parametrically smoothed estimate of `y` and `y'` given the `data`, where 
% each column is a feature of the timeseries at each tpoints, and the matrix is 
% with dimension  (length(tpoints),numberOfFeatuer). if the length of tpoints 
% and data timeseries are not same, the shortest will be used.
% 
% For kernels, the following exist: 
% - EpanechnikovKernel
% - UniformKernel
% - TriangularKernel 
% - QuarticKernel 
% - TriweightKernel
% - TricubeKernel
% - GaussianKernel
% - CosineKernel 
% - LogisticKernel
% - SigmoidKernel
% - SilvermanKernel

function [y,dy] = collocate_data(data,tpoints,kernelStr)

    arguments
        data (:,:) {mustBeNumeric,mustBeReal}
        tpoints (:, 1) {mustBeNumeric,mustBeReal}
        kernelStr (1,:) char {mustBeMember(kernelStr,{'EpanechnikovKernel',...
                                                      'UniformKernel',...
                                                      'TriangularKernel',...
                                                      'QuarticKernel',...
                                                      'TriweightKernel',...
                                                      'TricubeKernel',...
                                                      'GaussianKernel',...
                                                      'CosineKernel',...
                                                      'LogisticKernel',...
                                                      'SigmoidKernel',...
                                                      'SilvermanKernel'})} = 'EpanechnikovKernel'
    end

    nt = length(tpoints);
    nd = size(data,1);
    if nt > nd
        n = nd;
        tpoints = tpoints(1:n);
    elseif nt<nd
        n = nt;
        data = data(1:n,:);
    else
        n=nt;
    end

    y = zeros(size(data));
    dy = zeros(size(data));

    e1 = [1;0];
    e2 = [0;1;0];

    h = (n^(-1/5))*(n^(-3/35))*((log(n))^(-1/16));
    
    for i = 1:1:n
      T1 = horzcat(ones(n,1),reshape(tpoints-tpoints(i),[n,1]));
      T2 = horzcat(T1,reshape((tpoints-tpoints(i)).^2,[n,1]));

      % construct W
      w1 = zeros(1,n);
      for j = 1:1:n
          w1(j) = calculateKernel(kernelStr,(tpoints(j)-tpoints(i))/h)/h;
      end
      W = diag(w1);

      y(i,:) = e1'*((T1'*W*T1)\T1')*W*data;
      dy(i,:) = e2'*((T2'*W*T2)\T2')*W*data;

      % in Matlab, x=a\b is not same as x=inv(a)*b
      % y(i,:) = e1'*((T1'*W*T1)^(-1))*T1'*W*data;
      % dy(i,:) = e2'*((T2'*W*T2)^(-1))*T2'*W*data;
    end
end

% Define Kernel function

function kernel=calculateKernel(kernelStr, t)
    switch kernelStr
        case "EpanechnikovKernel"
            if abs(t) > 1
                kernel = 0;
            else
                kernel = 0.75*(1-t^2);
            end
        case "UniformKernel"
            if abs(t) > 1
                kernel = 0;
            else
                kernel = 0.5;
            end
        case "TriangularKernel"
            if abs(t) > 1
                kernel = 0;
            else
                kernel = (1-abs(t));
            end
        case "QuarticKernel"
            if abs(t) > 0
                kernel = 0;
            else
                kernel = (15*(1-t^2)^2)/16;
            end
        case "TriweightKernel"
            if abs(t)>0
                kernel = 0;
            else
                kernel = (35*(1-t^2)^3)/32;
            end
        case "TricubeKernel"
            if abs(t) > 0
                kernel = 0;
            else
                kernel = (70*(1-abs(t)^3)^3)/80;
            end
        case "CosineKernel"
            if abs(t)>0
                kernel = 0;
            else
                kernel = (pi*cos(pi*t/2))/4;
            end
        case "GaussianKernel"
            kernel = exp(-0.5*t^2)/(sqrt(2*pi));
        case "LogisticKernel"
            kernel = 1/(exp(t)+2+exp(-t));
        case "SigmoidKernel"
            kernel = 2/(pi*(exp(t)+exp(-t)));
        case "SilvermanKernel"
            kernel = sin(abs(t)/2+pi/4)*0.5*exp(-abs(t)/sqrt(2));
        otherwise
            disp('Wrong Kernel Function Name!!!')
    end
end