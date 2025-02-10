import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pdb
import matplotlib.gridspec as gridspec
from pylab import *
import xarray as xr
from Casicus import centroid, label_cluster

mdir = '/home/tompkins-archive/acasallas/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/'
scra = '/home/netapp-clima/scratch/acasallas/Idealized_plots/'

dista = 'free' #conv or free
df3 = pd.read_csv(mdir+'Data_to_plot/Distances_'+dista+'_3D.csv')
dfT = pd.read_csv(mdir+'Data_to_plot/Distances_'+dista+'_TKE.csv')
df2 = pd.read_csv(mdir+'Data_to_plot/Distances_'+dista+'_2D.csv')
xh = xr.open_dataset(mdir+'Data_to_plot/Eddy_viscosity_tot_0-1056.nc')

sta = 48
end = 10*24
jump = 24

tipo = 'all'

if tipo == '':
    size_num = {'CS_Tho3YSU':[],'CS_Tho3BL':[],'CS_ThoTYSU':[],'CS_ThoTBL':[],'CS_ThoYSU':[],'CS_ThoBL':[],
                'CN_Tho3YSU':[],'CN_Tho3BL':[],'CN_ThoTYSU':[],'CN_ThoTBL':[],'CN_ThoYSU':[],'CN_ThoBL':[]}
    runs = ['Tho3YSU','Tho3BL','ThoTYSU','ThoTBL','ThoYSU','ThoBL']
elif tipo == 'MPs':
    size_num = {'CS_ThoYSU':[],'CS_ThoBL':[],'CS_GCEYSU':[],'CS_GCEBL':[],
                'CS_MorYSU':[],'CS_MorBL':[],'CS_WSMYSU':[],'CS_WSMBL':[],
                'CS_ThoYDE':[],'CS_ThoBHE':[],
                'CN_ThoYSU':[],'CN_ThoBL':[],'CN_GCEYSU':[],'CN_GCEBL':[],
                'CN_MorYSU':[],'CN_MorBL':[],'CN_WSMYSU':[],'CN_WSMBL':[],
                'CN_ThoYDE':[],'CN_ThoBHE':[]}
    runs = ['ThoYSU','ThoBL','GCEYSU','GCEBL','MorYSU','MorBL','WSMYSU','WSMBL','ThoYDE','ThoBHE']
elif tipo == 'all':
    runs = []
    titles = []
    size_num = {}
    mps = ['Tho','GCE','WSM','Mor']
    subgs = ['3','T','']
    pbls = ['YSU','BL']
    subgt = ['-SM3-','-TKE-','-SM2-']
    for i,subg in enumerate(subgs):
        for mp in mps:
            #if mp == 'Tho' and subg == '':
            #    pbls = ['YSU','YDE', 'BL', 'BHE']
            #else:
            #    pbls = ['YSU','BL']
            for pbl in pbls:
                runs.append(mp+subg+pbl)
                titles.append(mp+subgt[i]+pbl)
                size_num['CN_'+mp+subg+pbl] = [] 
                size_num['CS_'+mp+subg+pbl] = []       

stencil = np.ones((1,1))

try:
    df = pd.read_csv(mdir+'Data_to_plot/Size_and_number_cores'+tipo+'.csv')
except:
    for i,run in enumerate(runs):
        print('Processing: '+run)
        ds = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
        ww = ds.W[sta:end,19,:,:]
        del(ds)
        for t in range(0,end-sta):
            cluster_mask = np.where((ww[t,:,:]>1),1,0)
            num_pix = cluster_mask.sum()
            cluster_labels = cluster_labels = label_cluster(cluster_mask, stencil)
            centers = centroid(cluster_mask, cluster_labels).round(0).astype(int) 
            size_num['CS_'+run].append((num_pix/centers.shape[0])*4)
            size_num['CN_'+run].append(centers.shape[0])  
    df = pd.DataFrame.from_dict(size_num)
    df.to_csv(mdir+'Data_to_plot/Size_and_number_cores'+tipo+'.csv')

dfs = pd.concat([df3,dfT,df2], axis = 1)
sizem = []
numm = []
dist = []
kh = np.array([])
dia = 5
for i,run in enumerate(runs):
    #### Core size
    sizem.append(df['CS_'+run].mean())
    #### Core number
    numm.append(df['CN_'+run][48:dia*24].mean()) 
    #### dclr
    dist.append(dfs[run][48:dia*24].mean()*2)
    #### Entrainment
    kh = np.append(kh, xh[run+'_XKH'][48:dia*24].mean(dim='elev').max())

####### Plotting and so on
fig = plt.figure(figsize=(12,8))
gs = GridSpec(2,4,left = 0.08, right = 0.9, hspace=0.25, wspace=0.5, top = 0.92, bottom = 0.08)
marker = ['o']*8 + ['*']*8 + ['s']*8
color = ['darkblue','blue','royalblue','steelblue','deepskyblue','teal','darkcyan','cyan',
         'darkblue','pink','steelblue','red','salmon','orange','darkcyan','cyan',
         'darkblue','pink','steelblue','red','salmon','orange','darkcyan','cyan']
ax = plt.subplot(gs[0,0:2])
for i,run in enumerate(runs):
    plt.scatter(dist[i], kh[i], label = titles[i], color = color[i], marker = marker[i])

#plt.axhline(7.5, linestyle = '--', color = 'k', linewidth = 0.7)
#plt.ylim(40,250)
#plt.xlim(60,90)
plt.yscale('log')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#plt.axvline(71, linestyle = '--', color = 'k', linewidth = 0.7)
plt.xlabel('Max. Free Conv. Dist. (km)')
plt.ylabel('Sub-grid horizontal mixing ($K_{h}$)')
plt.title('(a)', loc = 'left', fontweight = 'bold')

ax = plt.subplot(gs[0,2:4])
for i,run in enumerate(runs):
    plt.scatter(dist[i], sizem[i], label = titles[i], color = color[i], marker = marker[i])

plt.yscale('log')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('Max. Free Conv. Dist. (km)')
plt.ylabel('Mean conv. core area (km$^{2}$)')
plt.title('(b)', loc = 'left', fontweight = 'bold')

ax = plt.subplot(gs[1,1:3])
for i,run in enumerate(runs):
    plt.scatter(numm[i], sizem[i], label = titles[i], color = color[i], marker = marker[i])

plt.legend(ncol=3, frameon=False, loc = 'upper center', bbox_to_anchor = (1.1, 0.9))
plt.yscale('log')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('Number of convective cores')
plt.ylabel('Mean convective core area (km$^{2}$)')
plt.title('(c)', loc = 'left', fontweight = 'bold')
plt.savefig(scra+'Sub-grid_scale_mixing_dclr_v2.jpg', bbox_inches = 'tight', dpi = 300)
plt.savefig(scra+'Sub-grid_scale_mixing_dclr_v2.pdf', bbox_inches = 'tight', dpi = 300)
plt.show()


