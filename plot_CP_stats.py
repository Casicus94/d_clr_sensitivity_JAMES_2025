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
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/CP_characteristics/Total/'
med1 = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Coldpools/'
scr = '/home/netapp-clima/scratch/acasallas/'

##Runs
runs = ['ThoYSU','ThoBHE','GCEYSU','MorYSU','MorBL','WSMYSU','WSMBL','GCEBL','ThoBL','ThoYDE']
colors = ['purple','k','blue','darkcyan','cyan','darkred','red','salmon','orange','yellow']

varis = ['Buoyancy','Tanom','Thvanom','Qvanom','Qv','Tv','T2','MSE','W','U','V']

fig = plt.figure(figsize=(12,16)) #8,16
gs = GridSpec(4,3,left = 0.09, right = 0.95, hspace=0.1, wspace=0.3, top = 0.9, bottom = 0.08)

anom = sys.argv[1]

for i,run in enumerate(runs):
    print('-----> Starting with '+run)
    if anom == True:
        df = pd.read_csv(med+run+'_CPs_statistics.csv')
    else:
        df = pd.read_csv(med+run+'_CPs_statistics_no_anom.csv')
        
    titles = ['Buoyancy (1e-2 m s$^{-2}$)',u"T$^{'} / \overline{T}$ (1e-3)",
              u"\u03F4$_{v}^{'} / \overline{\u03F4_{v}}$ (1e-3)",
              "(0.61q$_{v}^{'}) / (1+0.61\overline{q_{v}}$) (1e-3)","QVAPOR (1e-3 g kg$^{-1}$)",
              'Tv (K)','T2 (K)','MSE/Cp (K)','W (1e-3 m s$^{-1}$)','U (m s$^{-1}$)','V (m s$^{-1}$)']
    cte = [1e2,1e3,1e3,1e3,1e3,1,1,1/1005,1e3,1,1]
    if anom == 'True':
        lim_i = [-0.25,-0.05,-0.2,-0.4,16.5,302.5,299.8,341,-5.75,-0.05,-0.05]
        lim_s = [0.55,0.5,0.7,0.15,17.7,302.9,300.25,344,1,0.05,0.05]
    else:
        lim_i = [-0.25,-0.05,-0.2,-0.4,16.5,302.5,299.8,341,-5.75,-0.05,-0.05]
        lim_s = [0.55,0.5,0.7,0.15,17.7,302.9,300.25,344,1,0.05,0.05]
    for j,var in enumerate(varis):
        ax = plt.subplot(gs[j])
        dbins = np.arange(1,80) 
        vari = df[var]
        obj,_ = np.histogram(df['Area'], bins=dbins, weights=vari*cte[j])
        num,_ = np.histogram(df['Area'], bins=dbins)
        mean_obj = obj/num
        plt.plot(mean_obj, label = run, color = colors[i])
        plt.ylabel(titles[j], fontsize = 10)
        if j == 8:
            plt.ylim(0, 13)
        if j < 8:
            ax.xaxis.set_major_formatter(NullFormatter())
        else:
            plt.xlabel('Area')
        if j == 10 and i == len(runs)-1:
            plt.legend(frameon=False,ncol=2,bbox_to_anchor=(1.75,0.8), loc = 'upper center') #1.18,1.3
if anom == 'True':
    plt.savefig(med1+'CPs_statistics_anom.jpg', bbox_inches='tight')
else:
    plt.savefig(med1+'CPs_statistics.jpg', bbox_inches='tight')
plt.show()
plt.close()

plt.figure(figsize=(8,8))
for i,run in enumerate(runs):
    mat_ar = pd.read_csv(med+run+'_CPs_area_radius.csv')
    mat_ar = mat_ar['Area'].values
    plt.hist(mat_ar,histtype='step', bins=np.linspace(1,40,10), label = run, color = colors[i])
plt.legend(frameon = False)
plt.ylabel('Frequency', fontweight = 'bold')
plt.xlabel('Area', fontweight = 'bold')
if anom == True:
    plt.savefig(med1+'CPs_area_anom.jpg', bbox_inches='tight')
else: 
    plt.savefig(med1+'CPs_area.jpg', bbox_inches='tight')
plt.show()
plt.close()


