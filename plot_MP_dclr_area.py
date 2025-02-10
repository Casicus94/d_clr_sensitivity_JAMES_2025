import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pdb
import matplotlib.gridspec as gridspec
from pylab import *
import xarray as xr
from Casicus import centroid, label_cluster

mdir = '/home/tompkins-archive/acasallas/'
med = '/home/tompkins-archive/acasallas/Data_to_plot/Total_CP/'
scra = '/home/netapp-clima/scratch/acasallas/Idealized_plots/'
obdr = '/home/tompkins-archive2/acasallas/Ideal_RH/'
med1 = '/home/tompkins-archive2/acasallas/CP_characteristics/Total/'

dista = 'free' #conv or free
df3 = pd.read_csv(mdir+'Data_to_plot/Distances_'+dista+'_3D.csv')
dfT = pd.read_csv(mdir+'Data_to_plot/Distances_'+dista+'_TKE.csv')
df2 = pd.read_csv(mdir+'Data_to_plot/Distances_'+dista+'_2D.csv')
cp_2 = pd.read_csv(mdir+'Data_to_plot/CP_intensity_2D.csv')
cp_t = pd.read_csv(mdir+'Data_to_plot/CP_intensity_TKE.csv')
df_rh = pd.read_csv(obdr+'RH_df.csv')
df = pd.read_csv(mdir+'Data_to_plot/Size_and_number_coresall.csv')

sta = 48
end = 10*24
jump = 24

runs = []
mps = ['Tho','GCE','WSM','Mor']
subgs = [''] # 'T' or ''
for subg in subgs:
    for mp in mps:
        if mp == 'Tho' and subg == '':
            pbls = ['YSU','YDE', 'BL', 'BHE']
        else:
            pbls = ['YSU','BL']
        for pbl in pbls:
            runs.append(mp+subg+pbl)

dfs = pd.concat([df3,dfT,df2], axis = 1)
dfcp = pd.concat([cp_2,cp_t], axis = 1) 
dist = []
numm = []
rh = []
cps = []
ww_cp = [] 
cp_ww = []

for i,run in enumerate(runs):
    #### Core number
    numm.append(df['CN_'+run][48:240].mean())
    #### dclr
    dist.append(dfs[run][48:240].mean()*2)
    #### Relative Humidity
    rh.append(df_rh[run].mean())
    #### CP intensity
    cps.append(dfcp[run][48:10*24].mean())
    #### W
    df_char = pd.read_csv(med1+run+'_CPs_statistics.csv')
    dbins = np.arange(1,80.01)
    vari = df_char['W']
    obj,_ = np.histogram(df_char['Area'].where(vari>0), bins=dbins, weights=vari.where(vari>0))
    num,_ = np.histogram(df_char['Area'].where(vari>0), bins=dbins)
    mean_obj = obj*60/num
    cp_ww.append(np.mean(mean_obj))
    ww_cp.append(np.mean(vari.where(vari>0))*24)

####### Plotting and so on
fig = plt.figure(figsize=(14,10))
gs = GridSpec(4,6,left = 0.09, right = 0.98, hspace=0.5, wspace=0.5, top = 0.92, bottom = 0.08)
if subgs[0] == 'T':
    marker = ['*']*8
    color = ['darkblue','pink','steelblue','red','salmon','orange','darkcyan','cyan']
    linestyles = ['-','--','-','--','--','--','-','-']
else:
    marker = ['s']*10
    color = ['darkblue','darkred','pink','darkgreen','steelblue','red','salmon','orange','darkcyan','cyan']
    linestyles = ['-','--','--','-','-','--','--','--','-','-']

ax = plt.subplot(gs[0:2,0:2])
for i,run in enumerate(runs):
    plt.scatter(dist[i], cps[i]*3.6, label = run, color = color[i], marker = marker[i])

if subgs[0] == '':
    plt.legend(ncol=5, frameon=False, loc = 'upper center', bbox_to_anchor = (1.125,1.225))
else:
    plt.legend(ncol=4, frameon=False, loc = 'upper center', bbox_to_anchor = (1.05,1.225))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('Max. Free Conv. Dist. (km)')
plt.ylabel('Cold Pool Average Intensity (m s$^{-1}$)')
plt.title('(a)', loc = 'left', fontweight = 'bold')

ax = plt.subplot(gs[0:2,2:4])
for i,run in enumerate(runs):
    plt.scatter(numm[i], cps[i]*3.6, label = run, color = color[i], marker = marker[i])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('Number of conv. cores')
plt.ylabel('Cold Pool Average Intensity (m s$^{-1}$)')
plt.title('(b)', loc = 'left', fontweight = 'bold')

ax = plt.subplot(gs[2:4,0:2])
for i,run in enumerate(runs):
    plt.scatter(dist[i], ww_cp[i], label = run, color = color[i], marker = marker[i])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('Max. Free Conv. Dist. (km)')
plt.ylabel('100m Vertical velocity at the gust fronts (m s$^{-1}$)')
plt.title('(c)', loc = 'left', fontweight = 'bold')

ax = plt.subplot(gs[2:4,2:4])
for i,run in enumerate(runs):
    plt.scatter(rh[i], cps[i]*3.6, label = run, color = color[i], marker = marker[i])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('Relative Humidity (%) outside Conv. cores')
plt.ylabel('Cold Pool Average Intensity (m s$^{-1}$)')
plt.title('(d)', loc = 'left', fontweight = 'bold')

ax = plt.subplot(gs[0:4,4:6])
for i,run in enumerate(runs):
    mat_ra = pd.read_csv(med+run+'_CPs_area_radius.csv')
    y,x=np.histogram(mat_ra['Radius'].values,bins=np.linspace(0,8,10))
    bin_centers = 0.5*(x[1:]+x[:-1])
    plt.plot(bin_centers*bin_centers*np.pi, y, label = run, color = color[i], linestyles=linestyles[i])
plt.legend(frameon = False)
plt.yscale('log')
plt.ylabel('Frequency')
plt.xlabel(' ColdPool Area (km$^{2}$)')
plt.title('(e)', loc = 'left', fontweight = 'bold')

plt.savefig(scra+'MP_dclr_'+subgs[0]+'_area.jpg', bbox_inches = 'tight', dpi = 300)
plt.savefig(scra+'MP_dclr_'+subgs[0]+'_area.pdf', bbox_inches = 'tight', dpi = 300)
plt.show()


