import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import pdb
import matplotlib.gridspec as GridSpec
from pylab import *
from matplotlib import ticker

user = '/home/tompkins-archive2/acasallas/CP_RZ/' 
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/'
scra = '/home/netapp-clima/scratch/acasallas/'

###########################################################################################
######## First lets calculate and save the area of the dry zone from the added CPs ########
###########################################################################################

##Here the labels, in comment for each day
labels = [[202,113,2366,2274,1147,1066,1275,1507,1056,122], 
          [2174,2904,34,1259,1372,2147,135,1097,1026,48], 
          [2020,80,1016,857], 
          [1272],
          [1119]]

areas = {'Zone1':np.array([]), 'Zone2':np.array([]), 'Zone3':np.array([]), 
         'Zone4':np.array([]), 'Zone5':np.array([])}
times = np.arange(24,241,24)

for i,time in enumerate(times):
    filename='label_objects_'+str(time)+'.nc'
    obj = xr.open_dataset(user+filename)
    for label in range(np.shape(labels)[0]):
        for step in range(24):
            try:
                tmp = obj.labels[step,:,:].where(obj.labels[step,:,:] == labels[label][i]).count().values*4
                if tmp == 0:
                    areas['Zone'+str(label+1)] = np.append(areas['Zone'+str(label+1)], np.nan)
                else:
                    areas['Zone'+str(label+1)] = np.append(areas['Zone'+str(label+1)], tmp)
            except:
                areas['Zone'+str(label+1)] = np.append(areas['Zone'+str(label+1)], np.nan)

##############################################################
######## Second lets start with the dry patch growth! ########
##############################################################
run = 'MorBL'
mdir = '/home/tompkins-archive/acasallas/'
ds_dom = xr.open_dataset(mdir+run+'_run/Domain_buo_comp_'+run+'.nc')
ds_dry = xr.open_dataset(mdir+run+'_run/Dry_patches_Lagrangian_'+run+'.nc')
stas = [2*24, 5*24, 7*24, 10*24, 12*24, 15*24, 17*24, 19*24]
end = [3*24, 6*24, 8*24, 11*24, 13*24, 16*24, 18*24, 20*24-1]

dry_patch = {'THv_anom':[],'Qv_anom':[],'Temp_anom':[],'Qv-T':[],'RTHRATEN':[],'RTHRATLW':[],'RTHRATLWC':[]}
varis = ['THv_anom','Qv_anom','Temp_anom']
for i,sta in enumerate(stas):
    for var in varis:
        dry_patch[var].append((ds_dry[var][sta:end[i],:] - ds_dom[var][sta:end[i],:]).mean(dim = 'Time'))

for i,sta in enumerate(stas):
    dry_patch['Qv-T'].append(np.array(dry_patch['Qv_anom'])[i,:] + np.array(dry_patch['Temp_anom'])[i,:])

#varis = ['Qv-T','Qv_anom','Temp_anom']

ds_dom = xr.open_dataset(mdir+run+'_run/Domain_verticals_'+run+'.nc')
press_anom = (ds_dry.P[0:20*24-1,:] + ds_dry.PB[0:20*24-1,:]) - np.array(ds_dom.P[0:20*24-1,:,0,0] + ds_dom.PB[0:20*24-1,:,0,0])
net = ds_dry.RTHRATEN[0:20*24,:] - np.array(ds_dom.RTHRATEN[0:20*24-1,:,0,0]) 
lw = ds_dry.RTHRATLW[0:20*24,:] - np.array(ds_dom.RTHRATLW[0:20*24-1,:,0,0])
lwc = ds_dry.RTHRATLWC[0:20*24,:] - np.array(ds_dom.RTHRATLWC[0:20*24-1,:,0,0])
sw = ds_dry.RTHRATSW[0:20*24,:] - np.array(ds_dom.RTHRATSW[0:20*24-1,:,0,0])
swc = ds_dry.RTHRATSWC[0:20*24,:] - np.array(ds_dom.RTHRATSWC[0:20*24-1,:,0,0])
press = ds_dom.P[0:20*24,:] + ds_dom.PB[0:20*24,:]
varis1 = ['RTHRATEN','RTHRATLW','RTHRATLWC']
for i,sta in enumerate(stas):
    for var in varis1:
        fluxx = ds_dry[var][sta:end[i],0:59].mean(dim='Time')- ds_dom[var][sta:end[i],0:59,0,0].mean(dim='XTIME')
        dry_patch[var].append(fluxx)

###############################################
######## Finally lets plot our results ########
###############################################
#z = [40,87.5,141,233,333,436,542,646,752,858,961,1060,1213,1421,1623,1825,2031,2241,2450,2657,2863,3068,3271,3472,3678,3887,4094,4297,4500,4705,4910]

z = [40,87.5,141,233,333,436,542,646,752,858,961, 1060,1213,1421,1623,1825,2031,2241,2450,2657,2863,3068,3271,3472,3678,3887,4094,4297,4500,4705,4910,5116,5471,5978,6488,6997,7500,8005,8513,9013,9509,10009,10503,10988,11468,11949,12425,12896,13165,13356,13829,14294,14754,15447,16374,17302,18242,19190,20148]

ylim = 14
days = [2,5,7,10,12,15,17,20]
#colores = [['plum','violet','mediumorchid','mediumpurple','blueviolet','darkviolet','purple','indigo'],
#            ['lightblue','cyan','darkcyan','royalblue','blue','mediumblue','darkblue','midnightblue'],
#            ['yellow','orange','sandybrown','lightsalmon','salmon','tomato','red','darkred'] ]

colores = [['lightblue','cyan','darkcyan','royalblue','plum','mediumorchid','blueviolet','purple'],['lightblue','cyan','darkcyan','royalblue','plum','mediumorchid','blueviolet','purple'],['lightblue','cyan','darkcyan','royalblue','plum','mediumorchid','blueviolet','purple']]

title = [u"(a) \u03F4$_{v}^{'} / \overline{\u03F4_{v}}$","(b) (0.61q$_{v}^{'}) / (1+0.61\overline{q_{v}}$)",u"(c) T$^{'} / \overline{T}$"]

fig = plt.figure(figsize=(12,8))
gs = GridSpec(1,6,left = 0.1, right = 0.95, hspace=0.45, wspace=0.15, top = 0.92, bottom = 0.08, width_ratios = [1,1,1,0.15,1,0.1])

for i,var in enumerate(varis):
    ax = plt.subplot(gs[i])
    for time in range(len(stas)):
        plt.plot(dry_patch[var][time][0:len(z)]*1e3,np.array(z)/1000, label='Day '+str(days[time]),color=colores[i][time])
    plt.axvline(0, linestyle=':', linewidth=0.5, color = 'k') 
    if i == 1: 
        plt.legend(frameon=False, fontsize = 8)
    if i == 0:
        plt.ylabel('Height (km)')
    else:
        ax.yaxis.set_major_formatter(NullFormatter())
    plt.ylim(0,ylim)
    plt.yticks(np.arange(0,12.1,2))
    plt.title(title[i], loc='left')
    plt.xlabel('*10$^{-3}$', loc = 'right') 

levels = np.arange(-40,-9.9,10)
levels = np.append(levels,np.arange(-9,-0.9,1))
levels = np.append(levels,np.arange(1,9.1,1))
levels = np.append(levels,np.arange(10,41,10))
ax = plt.subplot(gs[4])
im = plt.contourf(np.array(z)/1000,np.linspace(0,20,20*24-1), press_anom[:,0:len(z)], cmap = 'seismic',extend='both',levels = levels)
plt.ylabel('Days')
plt.xlabel('Height (km)')
plt.title('(d) Pressure anomaly', loc = 'left')
plt.xlim(0,ylim)
plt.xticks(np.arange(0,12.1,2))

ax = plt.subplot(gs[5])
cbar = plt.colorbar(im, cax = ax, ticks = levels)
cbar.set_label('hPa')
#plt.savefig(scra+'Dry_patch_growth_buo_'+run+'.jpg', bbox_inches='tight', dpi = 300)
#plt.savefig(scra+'Dry_patch_growth_buo_'+run+'.pdf', bbox_inches='tight', dpi = 300)
#plt.show()
plt.close()

######### Horizontal pressure hov
fig = plt.figure(figsize=(8,12))
gs = GridSpec(6,3,left = 0.1, right = 0.95, hspace=0.1, wspace=0.2, top = 0.92, bottom = 0.08, height_ratios = [1,1,0.2,0.7,0.1,0.1])

for i,var in enumerate(varis):
    ax = plt.subplot(gs[0:2,i])
    for time in range(len(stas)):
        plt.plot(dry_patch[var][time][0:len(z)]*1e3,np.array(z)/1000, label='Day '+str(days[time]),color=colores[i][time])
    plt.axvline(0, linestyle=':', linewidth=0.5, color = 'k')
    if i == 1:
        plt.legend(frameon=False, fontsize = 8)
    if i == 0:
        plt.ylabel('Height (km)')
    else:
        ax.yaxis.set_major_formatter(NullFormatter())
    plt.ylim(0,ylim)
    plt.yticks(np.arange(0,12.1,2))
    plt.title(title[i], loc='left')
    plt.xlabel('*10$^{-3}$', loc = 'right')

levels = np.arange(-40,-9.9,10)
levels = np.append(levels,np.arange(-9,-0.9,1))
levels = np.append(levels,np.arange(1,9.1,1))
levels = np.append(levels,np.arange(10,41,10))

ax = plt.subplot(gs[3,:])
im = plt.contourf(np.linspace(0,20,20*24-1), np.array(z)/1000, press_anom[:,0:len(z)].T, cmap = 'seismic',extend='both',levels = levels)
plt.xlabel('Days')
plt.ylabel('Height (km)')
plt.title('(d) Pressure anomaly', loc = 'left')
plt.ylim(0,ylim)
plt.yticks(np.arange(0,12.1,2))

ax = plt.subplot(gs[5,:])
cbar = plt.colorbar(im, cax = ax, ticks = levels, orientation = 'horizontal')
cbar.set_label('hPa')
plt.savefig(scra+'Dry_patch_growth_buo_hor_'+run+'.jpg', bbox_inches='tight', dpi = 300)
plt.savefig(scra+'Dry_patch_growth_buo_hor_'+run+'.pdf', bbox_inches='tight', dpi = 300)
plt.show()
plt.close()

pdb.set_trace()

######
#run = 'ThoYSU-full'
if run == 'ThoYSU-full':
    fig = plt.figure(figsize=(8,12))
    gs = GridSpec(4,4,left = 0.08, right = 0.9, hspace=0.45, wspace=0.15, top = 0.92, bottom = 0.08, width_ratios = [1,1,1,0.1])
    plot_var = [lw, lw-lwc, lwc]
    levs = [-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1] 
    titles = ['(a) Q$_{rad}$ Net','(b) Q$_{rad}$ LW CRE',' (c) Q$_{rad}$ LW Clear-Sky']
    for i,var in enumerate(plot_var):
        ax = plt.subplot(gs[i,0:3])
        im = plt.contourf(np.linspace(0,20,20*24-1), np.array(z)/1000, var[:,0:len(z)].T*86400, cmap = 'seismic', levels = levs, extend = 'both')
        #if i == 0:
        #    plt.ylabel('Days')
        #else:
        #    ax.yaxis.set_major_formatter(NullFormatter())
        plt.xlabel('Days')
        plt.ylabel('Height (km)')
        plt.title(titles[i], loc = 'left')
        plt.ylim(0,ylim)
        plt.yticks(np.arange(0,12.1,2))
    ax = plt.subplot(gs[0:3,3])
    cbar = plt.colorbar(im, cax = ax, ticks = levs)
    cbar.set_label('K day$^{-1}$')

    colors = ['purple', 'blue', 'red', 'salmon', 'orange']
    ax = plt.subplot(gs[3,:])
    for i in np.arange(1,6): 
        plt.plot(np.linspace(0,10,240), areas['Zone'+str(i)], label = 'Zone '+str(i), color = colors[i-1])

    plt.yscale('log')
    plt.ylabel('Dry patch Area (km$^{2}$)')
    plt.xlabel('Days')
    plt.title('(d)', loc ='left') 
    plt.legend(frameon=False)
    plt.xlim(0,ylim)
    plt.xticks(np.arange(0,12.1,2))
    #plt.savefig(scra+'Dry_patch_growth_rad_CP_init_hov_'+run+'.jpg', bbox_inches='tight', dpi = 300)
    #plt.savefig(scra+'Dry_patch_growth_rad_CP_init_hov_'+run+'.pdf', bbox_inches='tight', dpi = 300)
    plt.show()
    plt.close()
else:
    fig = plt.figure(figsize=(12,12))
    gs = GridSpec(4,7,left = 0.08, right = 0.9, hspace=0.45, wspace=0.3, top = 0.92, bottom = 0.08, width_ratios = [1,1,1,1,1,1,0.1])
    plot_var = [[lw, lw-lwc, lwc],[sw, sw-swc, swc]]
    levs = [-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1]
    titles = ['(a) Q$_{rad}$ LW Net','(b) Q$_{rad}$ LW CRE',' (c) Q$_{rad}$ LW Clear-Sky', '(d) Q$_{rad}$ SW Net','(e) Q$_{rad}$ SW CRE','(f) Q$_{rad}$ SW Clear-Sky']
    count = 0
    for i,pvar in enumerate(plot_var):
        for j,var in enumerate(pvar):
            if i == 0:
                ax = plt.subplot(gs[j,0:3])
            else: 
                ax = plt.subplot(gs[j,3:6])
            #im = plt.contourf(np.linspace(0,20,20*24-1), np.array(z)/1000, var[:,0:len(z)].T*86400, cmap = 'seismic', levels = levs, extend = 'both')
            im = plt.contourf(np.linspace(0,19,19*24-1), np.array(z)/1000, var.rolling(Time=24).mean()[24:,0:len(z)].T*86400, cmap = 'seismic', levels = levs, extend = 'both')
            if i == 0:
                plt.ylabel('Height (km)') 
            else:
                ax.yaxis.set_major_formatter(NullFormatter())
            plt.xlabel('Days')
            plt.title(titles[count], loc = 'left')
            plt.ylim(0,ylim)
            plt.yticks(np.arange(0,12.1,2))
            count +=1

    ax = plt.subplot(gs[3,0:6])
    #im = plt.contourf(np.linspace(0,20,20*24-1), np.array(z)/1000, net[:,0:len(z)].T*86400, cmap = 'seismic', levels = levs, extend = 'both')
    im = plt.contourf(np.linspace(0,19,19*24-1), np.array(z)/1000, net.rolling(Time=24).mean()[24:,0:len(z)].T*86400, cmap = 'seismic', levels = levs, extend = 'both')
    plt.ylabel('Height (km)')
    plt.xlabel('Days')
    plt.title('(g) Q$_{rad}$ Net', loc = 'left')
    plt.ylim(0,ylim)
    plt.yticks(np.arange(0,12.1,2))
    
    ax = plt.subplot(gs[0:4,6])
    cbar = plt.colorbar(im, cax = ax, ticks = levs)
    cbar.set_label('K day$^{-1}$')    
    #plt.savefig(scra+'Dry_patch_growth_rad_CP_init_hov_day_'+run+'_full.jpg', bbox_inches='tight', dpi = 300)
    #plt.savefig(scra+'Dry_patch_growth_rad_CP_init_hov_day_'+run+'_full.pdf', bbox_inches='tight', dpi = 300)
    #pdb.set_trace()
    plt.show()
    plt.close()

######
if run == 'ThoYSU-full':
    fig = plt.figure(figsize=(10,12))
    gs = GridSpec(4,3,left = 0.08, right = 0.9, hspace=0.45, wspace=0.1, top = 0.92, bottom = 0.08, width_ratios = [1,1,1])
else:
    fig = plt.figure(figsize=(10,9))
    gs = GridSpec(3,3,left = 0.08, right = 0.9, hspace=0.45, wspace=0.1, top = 0.92, bottom = 0.08, width_ratios = [1,1,1])
colores = ['lightblue','cyan','darkcyan','royalblue','plum','mediumorchid','blueviolet','purple']
for i,var in enumerate(varis1):
    ax = plt.subplot(gs[0:3,i])
    for time in range(len(stas)):
        if i == 1:
            plt.plot((dry_patch[var][time][0:len(z)]-dry_patch['RTHRATLWC'][time][0:len(z)])*86400,np.array(z)/1000, label='Day '+str(days[time]),color=colores[time])
        else:
            plt.plot(dry_patch[var][time][0:len(z)]*86400,np.array(z)/1000, label='Day '+str(days[time]),color=colores[time])
    plt.axvline(0, linestyle=':', linewidth=0.5, color = 'k')
    if i == 1:
        plt.legend(frameon=False, fontsize = 8)
    if i == 0:
        plt.ylabel('Height (km)')
    else:
        ax.yaxis.set_major_formatter(NullFormatter())
    plt.ylim(0,ylim)
    plt.yticks(np.arange(0,12.1,2))
    plt.title(titles[i], loc='left')
    plt.xlabel('K day$^{-1}$')
    
if run == 'ThoYSU-full':
    colors = ['purple', 'blue', 'red', 'salmon', 'orange']
    ax = plt.subplot(gs[3,:])
    for i in np.arange(1,6):
        plt.plot(np.linspace(0,10,240), areas['Zone'+str(i)], label = 'Zone '+str(i), color = colors[i-1])

    plt.yscale('log')
    plt.ylabel('Dry patch Area (km$^{2}$)')
    plt.xlabel('Days')
    plt.title('(d)', loc ='left')
    plt.legend(frameon=False)
    plt.xticks(np.arange(0,10.01))
    plt.xlim(0,10)

#plt.savefig(scra+'Dry_patch_growth_vert_rad_CP_init_'+run+'.jpg', bbox_inches='tight', dpi = 300)
#plt.savefig(scra+'Dry_patch_growth_vert_rad_CP_init_'+run+'.pdf', bbox_inches='tight', dpi = 300)
#plt.show()
plt.close()

plt.figure(figsize=(10,4))
colors = ['purple', 'blue', 'red', 'salmon', 'orange']
for i in np.arange(1,6):
    plt.plot(np.linspace(0,10,240), areas['Zone'+str(i)], label = 'Zone '+str(i), color = colors[i-1])

plt.yscale('log')
plt.ylabel('Dry patch Area (km$^{2}$)')
plt.xlabel('Days')
plt.legend(frameon=False)
plt.xticks(np.arange(0,10.01))
plt.xlim(0,10)
#plt.savefig(scra+'CP_init_experiments.jpg', bbox_inches='tight', dpi = 300)
#plt.savefig(scra+'CP_init_experiments.pdf', bbox_inches='tight', dpi = 300)
#plt.show()

