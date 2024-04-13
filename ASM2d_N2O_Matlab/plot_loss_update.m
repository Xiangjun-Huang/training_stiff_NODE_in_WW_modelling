       
% Update loss plot along with the increase of iteration 
function plot_loss_update(fLoss, lineLossTrain, loss, iter, titleText)
    figure(fLoss);
    currentLoss = double(loss);
    addpoints(lineLossTrain,iter,currentLoss);
    title(titleText);
    drawnow
end