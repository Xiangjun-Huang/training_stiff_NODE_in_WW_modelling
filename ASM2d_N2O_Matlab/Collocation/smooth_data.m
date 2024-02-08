% Smoothdata function 
% Gauss fitering (better than smoothingspline). 
% 
% _Important: try as possible as it can be to smooth the data firstly, as more 
% noise you bring into neural network, the more instable the training will be.

f=figure;
f.Position(3) = 2.5*f.Position(3);
f.Position(4) = 2*f.Position(4);
ySmooth = zeros(size(ySolNoise));
smoothWindow = 50; % set smoothing window length, adaptive length can be considered.
for k = 1:numFeature
    ySmooth(:,k) = smoothdata(ySolNoise(:,k),'gaussian',smoothWindow);
    subplot(4,6,k)
    plot(tHour, ySolNoise(:,k), '-.', tHour,ySmooth(:,k))
    title(nameFeature(k)) 
    xlim([0 6])
    ylim([0 inf])
    xticks(0:2:6)
    xticklabels({'0','2','4','6h'})
    ylabel(nameUnit(k))
    xlabel('')
    if k==1 
        legend('Noisy Y','Smooth Y','Location','southeast')
    end
    grid on
end
titleText = sprintf('%s%.0f%% %s%d%s %d %s',"Noisy Y (", noiseSigma*100, "\sigma) and Gaussian filtered Y (window=", smoothWindow,") with", numTimeSteps, "samples");
sgtitle(titleText);