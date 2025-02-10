import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pylab import *
import pdb

#--------------------------------------------
#Paths
mdir = '/home/tompkins-archive/acasallas/'
med = '/home/tompkins-archive/acasallas/Data_to_plot/Total_CP/'
med1 = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Coldpools/'
scr = '/home/netapp-clima/scratch/acasallas/'

#Runs
runs = []
mps = ['Tho','GCE','WSM','Mor']
subgs = ['T'] # 'T' or ''
for subg in subgs:
    for mp in mps:
        if mp == 'Tho' and subg == '':
            pbls = ['YSU','YDE', 'BL', 'BHE']
        else:
            pbls = ['YSU','BL']
        for pbl in pbls:
            runs.append(mp+subg+pbl)

if subgs[0] == '':
    colors = ['darkblue','darkred','pink','darkgreen','steelblue','red','salmon','orange','darkcyan','cyan'] 
else:
    colors = ['darkblue','pink','steelblue','red','salmon','orange','darkcyan','cyan']

plt.figure(figsize=(8,8))
for i,run in enumerate(runs):
    mat_ra = pd.read_csv(med+run+'_CPs_area_radius.csv')
    y,x=np.histogram(mat_ra['Radius'].values,bins=np.linspace(0,8,10))
    bin_centers = 0.5*(x[1:]+x[:-1])
    plt.plot(bin_centers, y, label = run, color = colors[i])
plt.legend(frameon = False)
plt.yscale('log')
plt.ylabel('Frequency', fontweight = 'bold')
plt.xlabel('Radius (km$^{2}$)', fontweight = 'bold')
plt.savefig(scr+'Coldpools_radius'+subgs[0]+'.jpg', bbox_inches='tight')
plt.show()
plt.close()

plt.figure(figsize=(8,8))
for i,run in enumerate(runs):
    mat_ra = pd.read_csv(med+run+'_CPs_area_radius.csv')
    y,x=np.histogram(mat_ra['Radius'].values,bins=np.linspace(0,8,10))
    bin_centers = 0.5*(x[1:]+x[:-1])
    plt.plot(bin_centers*np.pi, y, label = run, color = colors[i])
plt.legend(frameon = False)
plt.yscale('log')
plt.ylabel('Frequency', fontweight = 'bold')
plt.xlabel('Area (km$^{2}$)', fontweight = 'bold')
plt.savefig(scr+'Coldpools_area'+subgs[0]+'.jpg', bbox_inches='tight')
plt.show()
plt.close()


