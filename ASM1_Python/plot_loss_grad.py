from matplotlib import pyplot as plt

def plot_loss_grad(loss_list,grad_norm,lr_list):
    # Plot loss
    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax2 = ax1.twinx()
    ax1.plot(loss_list, 'b-', label ='Loss')
    ax2.plot(lr_list, 'g--', label = 'Learning rate')

    ax1.set_xlabel('iterations')
    ax1.set_ylabel('Log Loss', color='b')
    ax1.set_yscale('log')
    ax1.set_xlim((0,len(loss_list)))
    ax1.grid(which='both', linewidth=0.2)
    ax1.margins(x=0)
    ax1.tick_params(which="both", direction="in")
    ax2.tick_params(which="both", axis="y", direction="in")
    ax2.set_ylabel('Learning Rate', color='g')
    ax2.set_yscale('log')

    fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
    plt.show()

    # Plot gradient norm
    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax2 = ax1.twinx()
    ax1.plot(grad_norm, 'r', label = 'Gradient norm')
    ax2.plot(lr_list, 'g--', label = 'Learning rate')

    ax1.set_xlabel('iterations')
    ax1.set_ylabel('graident norm', color='r')
    ax1.set_yscale('log')
    ax1.set_xlim((0,len(grad_norm)))
    ax1.grid(which='both',linewidth=0.2)
    ax1.margins(x=0)
    ax1.tick_params(which="both", direction="in")
    ax2.tick_params(which="both", axis="y", direction="in")
    ax2.set_ylabel('Learning Rate', color='g')
    ax2.set_yscale('log')

    fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
    plt.show()
    return