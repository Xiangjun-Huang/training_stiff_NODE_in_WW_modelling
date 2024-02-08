% u′,u = collocate_data(data,tpoints,kernel=“EpanechnikovKernel”)
% 
% Computes a non-parametrically smoothed estimate of `u'` and `u` given the
% `data`, where each column is a snapshot of the timeseries at `tpoints[i]`.
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
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2631937/

% Define collocate_data function, data is a matrix with dimension
% (length(tpoints),numberOfFeatuer)
function [dy,y] = collocate_data(data,tpoints,kernelStr)

arguments
    data
    tpoints
    kernelStr = 'EpanechnikovKernel'
end

    e1 = [1;0];
    e2 = [0;1;0];
    n = length(tpoints);
    h = (n^(-1/5))*(n^(-3/35))*((log(n))^(-1/16));
    
    dy = zeros(size(data));
    y = zeros(size(data));
    for i = 1:1:n
      T1 = construct_t1(tpoints(i),tpoints);
      T2 = construct_t2(tpoints(i),tpoints);
      W = construct_w(tpoints(i),tpoints,h,kernelStr);
      Wd = W * data;
      WT1 = W * T1;
      WT2 = W * T2;
      dy(i,:) = (e2'*((T2'*WT2)\T2'))*Wd;
      y(i,:) = (e1'*((T1'*WT1)\T1'))*Wd;
    end
end

% define Kernel function
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
% Construct t1, t2, w function
function t1 = construct_t1(t,tpoints)
    t1 = horzcat(ones(length(tpoints),1),reshape(tpoints-t,[length(tpoints),1]));
end

function t2 = construct_t2(t,tpoints)
    t2 = horzcat(ones(length(tpoints),1),reshape([tpoints-t,(tpoints-t).^2],[length(tpoints),2]));
end

function w = construct_w(t,tpoints,h,kernelStr)
    n = length(tpoints);
    w1 = zeros(1,n);
    for i = 1:1:n
        w1(i) = calculateKernel(kernelStr,(tpoints(i)-t)/h)/h;
    end
    w = diag(w1);
end

