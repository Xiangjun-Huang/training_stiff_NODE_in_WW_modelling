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

# model initiation
noFeature, yMu, yStd, dyMu, dyStd, yMax, yMin, dyMax, dyMin, ytScale = estimate_mean_std(t, true_y)
model = model_initiation(noFeature, yMu, yStd, dyMu, dyStd, yMax, yMin, dyMax, dyMin, ytScale)

# collocation calculation
coll_y, coll_dy = collocate_data_torch(t, true_y)
# plot_data(t, coll_y, true_y, y_labels=['Coll Trajectory','True Trajectory'], active_only=True, show_RMSE=True)
# plot_data(t, coll_dy, true_dy, y_labels=['Coll derivative','True derivative'], active_only=True, show_RMSE=True,
#           y_name = ['$S_I\'$', '$S_S\'$', '$X_I\'$', '$X_S\'$', '$X_{BH}\'$', '$X_{BA}\'$', '$X_P\'$', '$S_O\'$', '$S_{NO}\'$', '$S_{NH}\'$', '$S_{ND}\'$', '$X_{ND}\'$', '$S_{ALK}\'$', '$S_{N_2}\'$', '$X_{INORG}\'$'],
#           unit_name = ['mg COD/(l.d)', 'mg COD/(l.d)', 'mg COD/(l.d)', 'mg COD/(l.d)', 'mg COD/(l.d)', 'mg COD/(l.d)', 'mg COD/(l.d)', '$mg O_2/(l.d)$', 'mg N/(l.d)', 'mg N/(l.d)', 'mg N/(l.d)', 'mg N/(l.d)', 'eq ALK/(l.d)', 'mg N/(l.d)', 'mg COD/(l.d)'])

# collocation training
model, loss_list, grad_norm, lr_list = collocation_training(t, coll_y, coll_dy, model)
pred_y = validation(model, t)
plot_data(t, pred_y, true_y, y_labels=['Prediction','True'], show_RMSE=True)
plot_loss_grad(loss_list, grad_norm, lr_list)

# NODE training
model, loss_list, grad_norm, lr_list = NODE_training(t, true_y, model)
pred_y = validation(model, t)
plot_data(t, pred_y, true_y, y_labels=['Prediction','True'], show_RMSE=True)
plot_loss_grad(loss_list, grad_norm, lr_list)