import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

eff_data = pd.read_csv('ASM_CSTR/ASM1_NODE_package-upload/data/efficiency_test.csv',index_col=0)
# print(eff_data)

mean_std = eff_data.agg(['mean','std'])
# print(mean_std)

plt.style.use('seaborn-v0_8-paper')
fig = plt.figure(figsize=(12, 6))

for i in range(7):
    ax = fig.add_subplot(2, 4, i+1)
    sns.histplot(eff_data, x=eff_data.columns[i], stat='percent',kde=True)
    
    # plot the lines
    m = mean_std.iloc[0,i]
    s = mean_std.iloc[1,i]
    ax.margins(x=0,y=0)
    ax.axvline(x=m, c='b', ls='-', lw=2.5)
    ax.axvline(x=m + s, c='orange', ls='--', lw=2.5)
    ax.axvline(x=m - s, c='orange', ls='--', lw=2.5)
    ax.set_title('Mean=%.1f \n Std=%.1f' %(m,s))
    ax.tick_params(which="both", direction="in")

plt.subplots_adjust(left=0.05, bottom=0.1, right=0.97, top=0.9, wspace=0.5,hspace=0.5)
plt.show()