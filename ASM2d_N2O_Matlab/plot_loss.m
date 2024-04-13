    
function [fLoss, lineLossTrain] = plot_loss()
    fLoss=figure;
    C = colororder;
    lineLossTrain = animatedline(Color=C(2,:));
    set(gca, 'YScale', 'log')
    % ylim([0 inf])
    xlabel("Iteration")
    ylabel("Loss")
    grid on
end

   