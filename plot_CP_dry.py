import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pylab import *
from Casicus import nan_out
import pdb

#--------------------------------------------
#Paths
mdir = '/home/tompkins-archive/acasallas/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/CP_characteristics/Dry/'
med1 = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Coldpools/'
scr = '/home/netapp-clima/scratch/acasallas/'

##Runs
runs = ['ThoYSU','ThoBHE','GCEYSU','MorYSU','MorBL','WSMYSU','WSMBL','GCEBL','ThoBL','ThoYDE']
colors = ['purple','k','blue','darkcyan','cyan','darkred','red','salmon','orange','yellow']
varis = ['Buoyancy','Tanom','Thvanom','Qvanom','W']

anom = sys.argv[1]
fig = plt.figure(figsize=(9,8)) #8,16
gs = GridSpec(2,3,left = 0.09, right = 0.95, hspace=0.15, wspace=0.3, top = 0.9, bottom = 0.08)
for i,run in enumerate(runs):
    print('-----> Starting with '+run)
    if anom == 'True':
        df = pd.read_csv(med+run+'_Dry_CP_charac_anom.csv')
        df_a = pd.read_csv(med+run+'_Dry_CP_area_anom.csv')
    else:
        df = pd.read_csv(med+run+'_Dry_CP_charac.csv')
        df_a = pd.read_csv(med+run+'_Dry_CP_area.csv')
    titles = ['Buoyancy (1e-2 m s$^{-2}$)',u"T$^{'} / \overline{T}$ (1e-3)",
              u"\u03F4$_{v}^{'} / \overline{\u03F4_{v}}$ (1e-3)",
              "(0.61q$_{v}^{'}) / (1+0.61\overline{q_{v}}$) (1e-3)",
              'W (1e-3 m s$^{-1}$)']
    area = df_a[df_a['Area'] != 0]
    cte = [1e2,1e3,1e3,1e3,1e3]
    if anom == 'True':
        lim_i = [-0.2,0.2,-0.15,-1,-20]
        lim_s = [0.6,1.2,0.5,0.1,2.5]
    else:
        lim_i = [-1,-1,-1,-1.2,-20]
        lim_s = [0.75,1.5,0.75,-0.2,2.5]
    for j,var in enumerate(varis):
        ax = plt.subplot(gs[j])
        dbins = np.arange(1,80) 
        vari = nan_out(df[var].values)
        area1 = nan_out(area['Area'].values)
        obj,_ = np.histogram(area1, bins=dbins, weights=vari*cte[j])
        num,_ = np.histogram(area1, bins=dbins)
        num = num*1.
        obj[obj==0.0] = np.nan
        num[num==0.0] = np.nan
        mean_obj = obj/num
        plt.plot(mean_obj, label = run, color = colors[i])
        plt.title(titles[j], fontsize = 10)
        #plt.ylim(lim_i[j], lim_s[j])
        if j < 3:
            ax.xaxis.set_major_formatter(NullFormatter())
        else:
            plt.xlabel('Area')
        if j == 4 and i == len(runs)-1:
            plt.legend(frameon=False,ncol=1,bbox_to_anchor=(1.75,0.8), loc = 'upper center') #1.18,1.3
if anom == 'True':
    plt.savefig(med1+'Dry_patch_CPs_statistics_anom.jpg', bbox_inches='tight')
else:
    plt.savefig(med1+'Dry_patch_CPs_statistics.jpg', bbox_inches='tight')
plt.show()
plt.close()

plt.figure(figsize=(8,8))
for i,run in enumerate(runs):
    if anom ==' True':
        mat_ar = pd.read_csv(med+run+'_Dry_CP_area_anom.csv')
    else:
        mat_ar = pd.read_csv(med+run+'_Dry_CP_area.csv')
    mat_ar = mat_ar['Area'].values
    plt.hist(mat_ar,histtype='step', bins=np.linspace(1,40,10), label = run, color = colors[i])
plt.legend(frameon = False)
plt.ylabel('Frequency', fontweight = 'bold')
plt.xlabel('Area', fontweight = 'bold')
if anom == True:
    plt.savefig(med1+'Dry_patch_CPs_area_anom.jpg', bbox_inches='tight')
else: 
    plt.savefig(med1+'Dry_patch_CPs_area.jpg', bbox_inches='tight')
#plt.show()
plt.close()


