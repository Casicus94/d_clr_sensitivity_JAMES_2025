import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pdb
import matplotlib.gridspec as gridspec
from pylab import *

tipo = '2D' #NE, DE, HE
runs = ['ThoYSU','ThoBHE','GCEYSU','MorYSU','MorBL','WSMYSU','WSMBL','GCEBL','ThoBL','ThoYDE']
colors = ['purple','k','blue','darkcyan','cyan','darkred','red','salmon','orange','yellow']

mdir = '/home/tompkins-archive/acasallas/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/'
scr = '/home/netapp-clima/scratch/acasallas/Cape/'

dista = 'free' #conv or free
dfd = pd.read_csv(mdir+'Data_to_plot/Distances_'+dista+'_'+tipo+'.csv')
dfcp = pd.read_csv(mdir+'Data_to_plot/CP_intensity_'+tipo+'.csv')
dfn = pd.read_csv(mdir+'Data_to_plot/Number_cores_'+tipo+'.csv')
dfs = pd.read_csv(mdir+'Data_to_plot/Size_and_number_coresMPs.csv')

sta = 2*24
end = 10*24

###Core size
plt.figure(figsize=(8,4))
for i,run in enumerate(runs):
    free = np.array(dfs['CS_'+run].rolling(24).mean())
    free = free[48+24:]
    x = np.linspace(sta/24,end/24,len(free[sta:end]))
    plt.plot(x, free[sta:end], label = run, color = colors[i], alpha=0.9)
    print('Core size '+run+': '+str(free[sta:end].mean()))
plt.legend(ncol=5, loc = 'upper center', bbox_to_anchor=(0.5,1), frameon = False)
plt.xticks(np.arange(2,(end/24.)+0.1,0.5))
plt.ylabel('Convective core area (km$^{2}$)')
#plt.show()
plt.close()

'''
import xarray as xr
###CAPE and CIN
jmp = 1
fig = plt.figure(figsize=(8,6)) #8,16
gs = GridSpec(2,1,left = 0.09, right = 0.95, hspace=0.1, wspace=0.3, top = 0.9, bottom = 0.08)
for i,run in enumerate(runs):
    print('CAPE and CIN, run: '+run)
    ds = xr.open_dataset(scr+'CAPE_CIN_'+run+'.nc')
    ds1 = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
    ww = ds1.W[48:,19,:,:]
    cape = np.array([])
    cin = np.array([])
    for t in range(0,15*24+1):
        cape = np.append(cape, ds.cape_2d[t,0,:,:].mean(dim=['south_north','west_east']))
        cin = np.append(cin, ds.cape_2d[t,1,:,:].mean(dim=['south_north','west_east']))
    dfca = pd.DataFrame(np.column_stack((cape,cin)), columns = ['CAPE','CIN'])
    cape = np.array(dfca['CAPE'].rolling(jmp).mean())
    cape = cape[48+jmp:]
    cin = np.array(dfca['CIN'].rolling(jmp).mean())
    cin = cin[48+jmp:]
    x = np.linspace(sta/24,end/24,len(cape[sta:end]))
    ax = plt.subplot(gs[0])
    plt.plot(x, cape[sta:end], label = run, color = colors[i], alpha=0.9)
    if i == 9:
        plt.legend(ncol=5, loc = 'upper center', bbox_to_anchor=(0.5,1.25), frameon = False)
        plt.title('(a) CAPE', x = 0.06, y = 0.88)
        plt.ylabel('J kg$^{-1}$')
    ax.xaxis.set_major_formatter(NullFormatter()) 
    ax = plt.subplot(gs[1])
    plt.plot(x, cin[sta:end], label = run, color = colors[i], alpha=0.9)
    if i == 9:
        plt.title('(b) CIN', x = 0.05, y = 0.88)
        plt.xlabel('Days')
        plt.ylabel('J kg$^{-1}$')
plt.savefig(med+'CAPE_CIN.jpg', bbox_inches = 'tight', dpi = 300)
plt.savefig(med+'CAPE_CIN.pdf', bbox_inches = 'tight', dpi = 300)
plt.show()
plt.close()
'''

######### Large plot here
fig = plt.figure(figsize=(14,8)) #8,16
gs = GridSpec(2,2,left = 0.10, right = 0.95, hspace=0.25, wspace=0.2, top = 0.9, bottom = 0.08)
##### Number of Cores first
roll = 18
alph = 0.7
ax = plt.subplot(gs[0])
for i,run in enumerate(runs):
    free = np.array(dfs['CN_'+run].rolling(roll).mean())
    free = free[48+roll:]
    x = np.linspace(sta/24,end/24,len(free[sta:end]))
    plt.plot(x, free[sta:end], label = run, color = colors[i], alpha=alph)
    print('Core number '+run+': '+str(free[sta:end].mean()))
plt.legend(ncol=5, loc = 'upper center', bbox_to_anchor=(1.1,1.2), frameon = False)
plt.xticks(np.arange(2,(end/24.)+0.1,0.5))
plt.title('(a)', x= 0.023, y = 0.90)
plt.ylabel('Number of convective cores', fontsize=10)
plt.xlabel('Days')
print('#########################')

##### Max. Free convective distance
ax = plt.subplot(gs[1])
for i,run in enumerate(runs):
    free = np.array(dfd[run].rolling(roll).mean())
    free = free[48+roll:]
    x = np.linspace(sta/24,end/24,len(free[sta:end]))
    plt.plot(x, free[sta:end]*2, label = run, color = colors[i], alpha=alph)
    print('Max distance '+run+': '+str(free[sta:end].mean()*2))
if dista == 'free':
    plt.ylabel('Max. Free Conv. Dist. (km)', fontsize=10)
else:
    plt.ylabel('Max. distance between conv. cores (km)')
plt.xticks(np.arange(2,(end/24.)+0.1,0.5))
plt.title('(b)', x= 0.023, y = 0.9)
plt.xlabel('Days')
print('#########################')

##### Coldpool intensity 
ax = plt.subplot(gs[2])
for i,run in enumerate(runs):
    free = np.array(dfcp[run].rolling(roll).mean())
    free = free[48+roll:]
    x = np.linspace(sta/24,end/24,len(free[sta:end]))
    plt.plot(x, free[sta:end]*3.6, label = run, color = colors[i], alpha=alph)
    print('Cold Pool int '+run+': '+str(free[sta:end].mean()*3.6))
plt.ylabel('Cold Pool Average Intensity (m s$^{-1}$)', fontsize=10)
plt.xticks(np.arange(2,(end/24.)+0.1,0.5))
plt.title('(c)', x= 0.023, y = 0.9)
plt.xlabel('Days')

med1 = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/CP_characteristics/Total/'

ax = plt.subplot(gs[3])
for i,run in enumerate(runs):
    df = pd.read_csv(med1+run+'_CPs_statistics.csv')
    dbins = np.arange(1,80.01)
    vari = df['W']
    obj,_ = np.histogram(df['Area'], bins=dbins, weights=vari)
    num,_ = np.histogram(df['Area'], bins=dbins)
    mean_obj = obj/num
    plt.plot(np.arange(2,160,2),mean_obj, label = run, color = colors[i], alpha = alph)
plt.ylabel('Vertical velocity (m s$^{-1}$)', fontsize = 10)
plt.xlabel('Area (km$^{2}$)')
plt.xticks(np.arange(0,161,10))
plt.title('(d)', x= 0.023, y = 0.9)
plt.savefig(med+'Dist_num_MPs_CPs_'+dista+'_'+str(int(sta/24))+'_to_'+str(int(end/24))+'.jpg', bbox_inches = 'tight')
plt.savefig(med+'Dist_num_MPs_CPs_'+dista+'_'+str(int(sta/24))+'_to_'+str(int(end/24))+'.pdf', bbox_inches = 'tight')
plt.show()
plt.close()
