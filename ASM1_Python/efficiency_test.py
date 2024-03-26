from IVP_data_generation import IVP_data_generation
from load_data import load_data
from plot_data import plot_data
from estimate_mean_std import estimate_mean_std
from model_initiation import model_initiation
from collocate_data_torch import collocate_data_torch
from collocation_training import collocation_training
from NODE_training import NODE_training
from validation import validation
from plot_loss_grad import plot_loss_grad
import time
from torch.nn import MSELoss
from torch import sqrt
import numpy as np

# generate trajectory data from mathemetical modelling
t, true_y, true_dy = IVP_data_generation()
# plot_data(t, true_y)
# plot_data(t, true_dy, active_only=True, y_labels=['True dy'], 
#            unit_name = ['mg COD/(l.d)', 'mg COD/(l.d)', 'mg COD/(l.d)', 'mg COD/(l.d)', 'mg COD/(l.d)', 
#                         'mg COD/(l.d)', 'mg COD/(l.d)', '$mg O_2/(l.d)$', 'mg N/(l.d)', 'mg N/(l.d)', 
#                         'mg N/(l.d)', 'mg N/(l.d)', 'eq ALK/(l.d)', 'mg N/(l.d)', 'mg COD/(l.d)'])

# load data from previously saved files
# t, true_y = load_data()
# plot_data(t, true_y)
# t, true_dy = load_data('ASM_CSTR/ASM1_NODE_package/data/true_dy.csv')
# plot_data(t, true_dy, active_only=True, y_labels=['True dy'], 
#            unit_name = ['mg COD/(l.d)', 'mg COD/(l.d)', 'mg COD/(l.d)', 'mg COD/(l.d)', 'mg COD/(l.d)', 
#                         'mg COD/(l.d)', 'mg COD/(l.d)', '$mg O_2/(l.d)$', 'mg N/(l.d)', 'mg N/(l.d)', 
#                         'mg N/(l.d)', 'mg N/(l.d)', 'eq ALK/(l.d)', 'mg N/(l.d)', 'mg COD/(l.d)'])

# estimation of mean and std
noFeature, yMu, yStd, dyMu, dyStd, yMax, yMin, dyMax, dyMin, ytScale = estimate_mean_std(t, true_y)

# collocation calculation
coll_y, coll_dy = collocate_data_torch(t, true_y)
# plot_data(t, coll_y, true_y, y_labels=['Coll y','True y'], active_only=True, show_RMSE=True)
# plot_data(t, coll_dy, true_dy, y_labels=['Coll dy','True dy'], active_only=True, show_RMSE=True)

mse = MSELoss()
num_test = 1
result_test = np.zeros((num_test,8))

for i in range(num_test):

    # ============incremental strategy =======================
    # model initiation
    model = model_initiation(noFeature, yMu, yStd, dyMu, dyStd, yMax, yMin, dyMax, dyMin, ytScale)
    
    # collocation training
    coll_start_time = time.time()
    model, loss_list, grad_norm, lr_list = collocation_training(t, coll_y, coll_dy, model)
    coll_end_time = time.time()
    # print('collocation training time = %.2f'%(coll_end_time-coll_start_time))
    pred_y = validation(model, t)
    RMSE = sqrt(mse(pred_y, true_y)).numpy()
    # plot_data(t, pred_y, true_y, y_labels=['Pred','True'], show_RMSE=True)
    # plot_loss_grad(loss_list, grad_norm, lr_list)
    result_test[i,0] = i
    result_test[i,1] = coll_end_time-coll_start_time
    result_test[i,2] = RMSE

    # NODE training
    node_start_time = time.time()
    model, loss_list, grad_norm, lr_list = NODE_training(t, true_y, model)
    node_end_time = time.time()
    # print('node training time = %.2f'%(node_end_time-node_start_time))
    pred_y = validation(model, t)
    RMSE = sqrt(mse(pred_y, true_y)).numpy()
    # plot_data(t, pred_y, true_y, y_labels=['Pred','True'], show_RMSE=True)
    # plot_loss_grad(loss_list, grad_norm, lr_list)
    result_test[i,3] = node_end_time-node_start_time
    result_test[i,4] = RMSE
    result_test[i,5] = (coll_end_time-coll_start_time) + (node_end_time-node_start_time)

    # ====================node only =====================
    # model initiation
    model = model_initiation(noFeature, yMu, yStd, dyMu, dyStd, yMax, yMin, dyMax, dyMin, ytScale) 

    # NODE training
    node_start_time = time.time()
    model, loss_list, grad_norm, lr_list = NODE_training(t, true_y, model)
    node_end_time = time.time()
    # print('node training time = %.2f'%(node_end_time-node_start_time))
    pred_y = validation(model, t)
    RMSE = sqrt(mse(pred_y, true_y)).numpy()
    # plot_data(t, pred_y, true_y, y_labels=['Pred','True'], show_RMSE=True)
    # plot_loss_grad(loss_list, grad_norm, lr_list)
    result_test[i,6] = node_end_time-node_start_time
    result_test[i,7] = RMSE

np.savetxt(fname='ASM_CSTR/ASM1_NODE_package-upload/result/efficiency_test.csv', X=result_test, fmt='%.2f', delimiter=',',
           header=','.join(['No','Coll Time','Coll RMSE','NODE Time','NODE RMSE','Incre Time','NODE only T', 'NODE only RMSE' ]))