import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pdb
import matplotlib.gridspec as gridspec
from pylab import *
import xarray as xr
from Casicus import centroid, label_cluster
import pdb

med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/'
mdir = '/home/tompkins-archive/acasallas/'
scr = '/home/netapp-clima/scratch/acasallas/'

sta = 48
end = 24*44

#Convective distance data
dista = 'free' #conv or free
df3 = pd.read_csv(mdir+'Data_to_plot/Distances_'+dista+'_3D.csv')
dfT = pd.read_csv(mdir+'Data_to_plot/Distances_'+dista+'_TKE.csv')
df2 = pd.read_csv(mdir+'Data_to_plot/Distances_'+dista+'_2D.csv')
iqr = pd.read_csv(mdir+'Data_to_plot/IQR_data_tot.csv')
cp3 = pd.read_csv(mdir+'Data_to_plot/CP_intensity_3D.csv')
cpT = pd.read_csv(mdir+'Data_to_plot/CP_intensity_TKE.csv')
cp2 = pd.read_csv(mdir+'Data_to_plot/CP_intensity_2D.csv')

mps = ['Tho','GCE','WSM','Mor']
subgs = ['3','T','']
subgt = ['-SM3-','-TKE-','-SM2-']

print('File ready, lets plot')
names = []
titles = []
for i,subg in enumerate(subgs):
    for mp in mps:
        if mp == 'Tho' and subg == '':
            pbls = ['YSU','YDE','YHI','YTB','BL','BHE','BDI','BHB']
        else:
            pbls = ['YSU','BL']
        for pbl in pbls:
            names.append(mp+subg+pbl)
            titles.append(mp+subgt[i]+pbl)

iqr_max = []
dista = []
cps = []
dia = 5
dis_dia = 5
for i,name in enumerate(names): 
    iqr_max.append(iqr[name][-dia*24:].mean())
    dista.append(df2[names[i]][48:dis_dia*24].mean()*2)
    cps.append(cp2[names[i]][48:dis_dia*24].mean())

####### Plotting and so on
fig = plt.figure(figsize=(11,7))
gs = GridSpec(1,2,left = 0.08, right = 0.9, hspace=0.45, wspace=0.05, top = 0.92, bottom = 0.08)

color = ['darkblue','blue','royalblue','steelblue','deepskyblue','green','darkcyan','cyan'] 
marker = ['o']*8
ax = plt.subplot(gs[1])
for i,name in enumerate(names[0:8]):
    plt.scatter(dista[i],iqr_max[i], s=cps[i]*10, label = titles[i], color = color[i], marker = marker[i])

plt.axhline(7.5, linestyle = '--', color = 'k', linewidth = 0.7)
plt.legend(ncol=2, frameon=False, loc = 'lower center', bbox_to_anchor=(0.7,0.18))
plt.ylim(0,40)
ax.yaxis.set_major_formatter(NullFormatter())
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)
plt.xlabel('Max. Free Conv. Dist (km)')

color = ['darkblue','darkred','steelblue','red','salmon','orange','darkcyan','cyan',
         'darkblue','darkred','blue','magenta','pink','darkgreen','coral','royalblue','steelblue','red',
         'salmon','orange','darkcyan','cyan']
marker = ['*']*8 + ['s']*14 
ax = plt.subplot(gs[0])
for i,name in enumerate(names[8:]):
    plt.scatter(dista[i+8],iqr_max[i+8], s=cps[i+8]*25, label = titles[i+8], color = color[i], marker = marker[i])

plt.axhline(7.5, linestyle = '--', color = 'k', linewidth = 0.7)
plt.ylim(0,40)
plt.xlim(60,90)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.axvline(73, linestyle = '--', color = 'k', linewidth = 0.7)
plt.legend(ncol=1, frameon=False, bbox_to_anchor=(0.18,1.02), loc = 'upper center')
plt.xlabel('Max. Free Conv. Dist. (km)')
plt.ylabel('TCWV-NMTP (mm)')
plt.savefig(scr+'Scatter_dist_CP_intensity_v2.jpg', bbox_inches = 'tight', dpi = 300)
plt.savefig(scr+'Scatter_dist_CP_intensity_v2.pdf', bbox_inches = 'tight', dpi = 300)
plt.show()


