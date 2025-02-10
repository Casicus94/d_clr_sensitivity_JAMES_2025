import pandas as pd
import numpy as np
import seaborn as sns
import pdb
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pylab import *

mdir = '/home/tompkins-archive/acasallas/'
scr = '/home/netapp-clima/scratch/acasallas/'

columns = ['Tho-SM2-BLY_-10','Tho-SM2-BLY_0','GCE-SM2-BLY_-10','GCE-SM2-BLY_0','Tho-TKE-BLY_-10','Tho-TKE-BLY_0','GCE-TKE-BLY_-10','GCE-TKE-BLY_0']
##### dclr
df = pd.read_csv(mdir+'Data_to_plot/BLY_dclr.csv')
dat_dclr = pd.DataFrame(np.column_stack((df.ThoBLY[0:234],df.ThoBLY[240:],df.GCEBLY[0:234],df.GCEBLY[240:],df.ThoTBLY[0:234],df.ThoTBLY[240:],df.GCETBLY[0:234],df.GCETBLY[240:])), columns = columns)

##### Core Number
df_n = pd.read_csv(mdir+'Data_to_plot/BLY_core_num.csv')
dat_num = pd.DataFrame(np.column_stack((df_n.ThoBLY[0:234],df_n.ThoBLY[240:],df_n.GCEBLY[0:234],df_n.GCEBLY[240:],df_n.ThoTBLY[0:234],df_n.ThoTBLY[240:],df_n.GCETBLY[0:234],df_n.GCETBLY[240:])), columns = columns)

##### CP intensity
df_c = pd.read_csv(mdir+'Data_to_plot/BLY_cp_int.csv')
dat_cpi = pd.DataFrame(np.column_stack((df_c.ThoBLY[100:200],df_c.ThoBLY[240:340],df_c.GCEBLY[100:200],df_c.GCEBLY[240:340],df_c.ThoTBLY[100:200],df_c.ThoTBLY[240:340],df_c.GCETBLY[100:200],df_c.GCETBLY[240:340])), columns = columns)

##### CP convergence
df_w = pd.read_csv(mdir+'Data_to_plot/BLY_cp_ww.csv')
dat_ww = pd.DataFrame(np.column_stack((df_w.ThoBLY[0:234],df_w.ThoBLY[240:],df_w.GCEBLY[0:234],df_w.GCEBLY[240:],df_w.ThoTBLY[0:234],df_w.ThoTBLY[240:],df_w.GCETBLY[0:234],df_w.GCETBLY[240:])), columns = columns)

##### PBL
df_p = pd.read_csv(mdir+'Data_to_plot/BLY_PBLH.csv')
dat_pbl = pd.DataFrame(np.column_stack((df_p.ThoBLY[0:237],df_p.ThoBLY[240:],df_p.GCEBLY[0:237],df_p.GCEBLY[240:],df_p.ThoTBLY[0:237],df_p.ThoTBLY[240:],df_p.GCETBLY[0:237],df_p.GCETBLY[240:])), columns = columns)

############### Plotting ###############
my_palette = {'Tho-SM2-BLY_-10':'lightblue', 'Tho-SM2-BLY_0':'lightblue', 'GCE-SM2-BLY_-10':'azure', 'GCE-SM2-BLY_0':'azure', 'Tho-TKE-BLY_-10':'mediumpurple','Tho-TKE-BLY_0':'mediumpurple', 'GCE-TKE-BLY_-10':'darkcyan','GCE-TKE-BLY_0':'darkcyan'}

fig = plt.figure(figsize=(12,10))
gs = GridSpec(4,2,left = 0.09, right = 0.98, hspace=0.18, wspace=0.15, top = 0.95, bottom = 0.08, height_ratios = [1,1,0.45,1])

ax = plt.subplot(gs[0,0])
sns.boxplot(data=dat_dclr*2, showmeans=True, palette = my_palette, meanprops={"marker":"s","markerfacecolor":"white", "markeredgecolor":"k"})
plt.ylabel('Max. Free Conv. Dist. (km)')
plt.title('(a)', loc = 'left')
ax.xaxis.set_major_formatter(NullFormatter())

ax = plt.subplot(gs[0,1])
sns.boxplot(data=dat_num, showmeans=True, palette = my_palette, meanprops={"marker":"s","markerfacecolor":"white", "markeredgecolor":"k"})
plt.ylabel('Number of Cores')
plt.title('(b)', loc = 'left')
ax.xaxis.set_major_formatter(NullFormatter())

ax = plt.subplot(gs[1,0])
sns.boxplot(data=dat_cpi, showmeans=True, palette = my_palette, meanprops={"marker":"s","markerfacecolor":"white", "markeredgecolor":"k"})
plt.ylabel('Cold Pool Intensity (m s$^{-1}$)')
plt.title('(c)', loc = 'left')
plt.xticks(rotation=90)

ax = plt.subplot(gs[1,1])
sns.boxplot(data=dat_ww, showmeans=True, palette = my_palette, meanprops={"marker":"s","markerfacecolor":"white", "markeredgecolor":"k"})
plt.ylabel('$W_{100}$ at the gust fronts (m s$^{-1}$)')
plt.title('(d)', loc = 'left')
plt.xticks(rotation=90)

ax = plt.subplot(gs[3,:])
sns.boxplot(data=dat_pbl, showmeans=True, palette = my_palette, meanprops={"marker":"s","markerfacecolor":"white", "markeredgecolor":"k"})
plt.ylabel('PBLH (m)')
plt.title('(e)', loc = 'left')
plt.savefig(scr+'BLY_experiments_v2.jpg', bbox_inches = 'tight', dpi = 300)
plt.savefig(scr+'BLY_experiments_v2.pdf', bbox_inches = 'tight', dpi = 300)
plt.show()
plt.close()

fig = plt.figure(figsize=(12,8))
gs = GridSpec(3,2,left = 0.09, right = 0.98, hspace=0.18, wspace=0.15, top = 0.95, bottom = 0.08, height_ratios = [1,1,0.2])

ax = plt.subplot(gs[0,0])
sns.boxplot(data=dat_dclr*2, showmeans=True, palette = my_palette, meanprops={"marker":"s","markerfacecolor":"white", "markeredgecolor":"k"})
plt.ylabel('Max. Free Conv. Dist. (km)')
plt.title('(a)', loc = 'left')
ax.xaxis.set_major_formatter(NullFormatter())

ax = plt.subplot(gs[0,1])
sns.boxplot(data=dat_num, showmeans=True, palette = my_palette, meanprops={"marker":"s","markerfacecolor":"white", "markeredgecolor":"k"})
plt.ylabel('Number of Cores')
plt.title('(b)', loc = 'left')
ax.xaxis.set_major_formatter(NullFormatter())

ax = plt.subplot(gs[1,0])
sns.boxplot(data=dat_ww, showmeans=True, palette = my_palette, meanprops={"marker":"s","markerfacecolor":"white", "markeredgecolor":"k"})
plt.ylabel('$W_{100}$ at the gust fronts (m s$^{-1}$)')
plt.title('(c)', loc = 'left')
plt.xticks(rotation=90)

ax = plt.subplot(gs[1,1])
sns.boxplot(data=dat_pbl, showmeans=True, palette = my_palette, meanprops={"marker":"s","markerfacecolor":"white", "markeredgecolor":"k"})
plt.ylabel('PBLH (m)')
plt.title('(d)', loc = 'left')
plt.xticks(rotation=90)
plt.savefig(scr+'BLY_experiments_no_cp_int_v2.jpg', bbox_inches = 'tight', dpi = 300)
plt.savefig(scr+'BLY_experiments_no_cp_int_v2.pdf', bbox_inches = 'tight', dpi = 300)
plt.show()

