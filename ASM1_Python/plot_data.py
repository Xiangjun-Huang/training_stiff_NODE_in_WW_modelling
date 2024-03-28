### Plot result, not show inactive components ###
from matplotlib import pyplot as plt
import torch

def plot_data(t, *ydata, y_labels=['True data'], active_only=True, show_RMSE = False,
    y_name = ['$S_I$', '$S_S$', '$X_I$', '$X_S$', '$X_{BH}$', '$X_{BA}$', '$X_P$', '$S_O$', '$S_{NO}$', '$S_{NH}$', '$S_{ND}$', '$X_{ND}$', '$S_{ALK}$', '$S_{N_2}$', '$X_{INORG}$'],
    unit_name = ['mg COD/l', 'mg COD/l', 'mg COD/l', 'mg COD/l', 'mg COD/l', 'mg COD/l', 'mg COD/l', '$mg O_2/l$', 'mg N/l', 'mg N/l', 'mg N/l', 'mg N/l', 'eq ALK/l', 'mg N/l', 'mg COD/l']):

    no_ydata = len(ydata)

    if show_RMSE:
        if no_ydata >= 2:
            # calcuate RMSE of first two inputs of ydata only
            loss_mse = torch.nn.MSELoss(reduction='mean')
            RMSE = torch.sqrt(loss_mse(ydata[0], ydata[1])).numpy()
        else:
            show_RMSE = False
    
    # plot 
    plt.style.use('seaborn-v0_8-paper')
    # plt.style.use('ggplot')

    # not show inactive components
    if active_only:
        selected_feature = (1,3,4,5,6,8,9,10,11,12,13)
        plot_row = 2
        plot_col = 6
        fig_size = (12,3.5)
    else: 
        selected_feature = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14)
        plot_row = 3
        plot_col = 6
        fig_size = (12, 5.5)

    fig = plt.figure(figsize=fig_size)

    j = 1
    for i in selected_feature:
        ax = fig.add_subplot(plot_row, plot_col, j)
        k = 0
        for y in ydata:
            ax.plot((t.detach().cpu().numpy()) * 24, y[:, 0][:, i].detach().cpu().numpy(),label=y_labels[k], lw=1)
            k=k+1
        if show_RMSE:
            ax.set_title('%s | RMSE = %.2f' % (y_name[i],(torch.sqrt(loss_mse(ydata[0][0:,0][:,i], ydata[1][0:,0][:,i]))).numpy() ), pad=4)
        else:
            ax.set_title(y_name[i])
        ax.margins(x=0,y=0)
        ax.set_xlabel('hour', fontsize=7)
        ax.set_ylabel(unit_name[i],fontsize=7) 
        ax.xaxis.set_label_coords(1, -0.15)
        ax.grid(which='both', linewidth=0.2)
        ax.tick_params(which="both", direction="in")

        j=j+1

    plt.subplots_adjust(left=0.07, bottom=0.1, right=0.97, top=0.9, wspace=0.5,hspace=0.5)
    lines, labels = fig.axes[-1].get_legend_handles_labels()
    fig.legend(lines, labels, loc='lower right', bbox_to_anchor=(0.97, 0.1))

    if show_RMSE:
        plt.text(1.4, 0.8, 'Overall RMSE = %.2f' % (RMSE), transform=ax.transAxes)

    plt.show()

    return