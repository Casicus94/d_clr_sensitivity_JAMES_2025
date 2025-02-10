import numpy as np
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.gridspec as gridspec
from pylab import *
import Casicus as casi
import pdb


#Path
mdir = '/home/tompkins-archive/acasallas/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Dry_composites/Verticals/'

#Run
runs = ['ThoYSU','GCEYSU','MorYSU','MorBL']

#Constants
g = 9.81
rho = 1.225
Rd=287.05
Cp=1005.0
Lv=2.5e6
Ls = 2.83e6
Lf = Ls - Lv

colors = ['purple','blue','darkcyan','cyan']
names = ['Qv (1e3 g kg$^{-1}$)','Ql (1e4 g kg$^{-1}$)','Qi (1e5 g kg$^{-1}$)',
         'Qs (1e4 g kg$^{-1}$)','Qr (1e4 g kg$^{-1}$)','MSE/Cp (K)','Buoyancy (m s$^{-2}$)']
nam_rad = ['RTHRATEN (K day$^{-1}$)','RTHRATLW (K day$^{-1}$)','RTHRATLWC (K day$^{-1}$)',
           'LWCRE (K day$^{-1}$)','RTHRATSW (K day$^{-1}$)','RTHRATSWC (K day$^{-1}$)']

z = [40,87.5,141,233,333,436,542,646,752,858,961,
       1060,1213,1421,1623,1825,2031,2241,2450,2657,
       2863,3068,3271,3472,3678,3887,4094,4297,4500,
       4705,4910,5116,5471,5978,6488,6997,7500,8005,
       8513,9013,9509,10009,10503,10988,11468,11949,
       12425,12896,13165,13356,13829,14294,14754,15447,
       16374,17302,18242.6,19190.13,20148.62,21084.01,
       21990.93,22890.83]

stas = [2*24, 2*24, 5*24, 10*24, 15*24]
end = [4*24, 3*24, 10*24, 15*24, 20*24] 

for i,sta in enumerate(stas):
    fig,(ax)=plt.subplots(nrows=1,ncols=7,figsize=(22,8))
    plt.subplots_adjust(hspace=0.0,wspace=0.25,top=0.9,right=0.98,left=0.05)
    for j,run in enumerate(runs):
        print('-----> Starting with '+run)
        #### Water first
        print('Plotting water species')
        ds_dry = xr.open_dataset(mdir+run+'_run/Dry_patches_Lagrangian_'+run+'.nc')
        ### First the variables
        print('Water species')
        #Water
        qv_d = ds_dry.QVAPOR[sta:end[i],:]
        ql_d = ds_dry.QCLOUD[sta:end[i],:]
        qi_d = ds_dry.QICE[sta:end[i],:]
        qs_d = ds_dry.QSNOW[sta:end[i],:]
        qr_d = ds_dry.QRAIN[sta:end[i],:]
        #MSE
        press_d = ds_dry.P[sta:end[i],:] + ds_dry.PB[sta:end[i],:]
        theta_d = ds_dry.T[sta:end[i],:] + 300
        temp_d = theta_d*np.power(press_d/1.e5,Rd/Cp)
        tempv_d=(1.0+0.61*qv_d)*temp_d # tempv is virtual T
        thv_d = tempv_d*np.power(1.e5/press_d*100,Rd/Cp)
        thv_mean_d = thv_d.mean(dim=['bottom_top'])
        thv_pert_d = thv_d - thv_mean_d
        mse_d = Cp*temp_d + g*np.array(z) + Lv*qv_d - Lf*qi_d
        buo_d = g*(thv_pert_d/thv_mean_d)
        del(theta_d);del(temp_d);del(tempv_d)
        del(thv_d);del(thv_mean_d);del(thv_pert_d)
        ########################
        #          Dry         #
        ########################
        to_plot = [qv_d*1e3,ql_d*1e4,qi_d*1e5,qs_d*1e4,qr_d*1e4,mse_d/Cp,buo_d]
        fs = 12
        for iplot,plot in enumerate(to_plot):
            ax[iplot].plot(plot.mean(dim='Time'),z, color = colors[j], label = run)
            if iplot == 0:
                ax[iplot].set_ylabel('Pressure (hPa)', fontsize = fs)
                if j == len(runs)-1:
                    ax[iplot].legend(frameon=False)
            else:
                ax[iplot].yaxis.set_major_formatter(NullFormatter())
            ax[iplot].set_title(names[iplot], fontsize = fs)
    plt.savefig(med+'Dry_compos_hum_vert_'+str(int(sta/24))+'_'+str(int(end[i]/24))+'.jpg', bbox_inches = 'tight')
    #plt.show()
    plt.close('all')

    fig,(ax)=plt.subplots(nrows=1,ncols=6,figsize=(22,8))
    plt.subplots_adjust(hspace=0.0,wspace=0.25,top=0.9,right=0.98,left=0.05)
    for j,run in enumerate(runs):
        print('-----> Starting with '+run)
        ds_dry = xr.open_dataset(mdir+run+'_run/Dry_patches_Lagrangian_'+run+'.nc')    
        ##### Rad Cooling
        print('Rad Cooling')
        #Rad
        press_d = ds_dry.P[sta:end[i],:] + ds_dry.PB[sta:end[i],:]
        net_d = ds_dry.RTHRATEN[:,:]
        LW_d = ds_dry.RTHRATLW[:,:]
        LWC_d = ds_dry.RTHRATLWC[:,:]
        LWCRE_d = LW_d - LWC_d
        SW_d = ds_dry.RTHRATSW[:,:]
        SWC_d = ds_dry.RTHRATSWC[:,:]
        ########################
        #          Dry         #
        ########################
        to_plot = [net_d,LW_d,LWC_d,LWCRE_d,SW_d,SWC_d]
        for iplot,plot in enumerate(to_plot):
            ax[iplot].plot(plot.mean(dim='Time')*86400,z, color = colors[j], label = run)
            if iplot == 0:
                ax[iplot].set_ylabel('Pressure (hPa)', fontsize = fs)
                if j == len(runs)-1:
                    ax[iplot].legend(frameon=False)
            else:
                ax[iplot].yaxis.set_major_formatter(NullFormatter())
            ax[iplot].set_title(nam_rad[iplot], fontsize = fs)
    plt.savefig(med+'Dry_compos_rad_vert_'+str(int(sta/24))+'_'+str(int(end[i]/24))+'.jpg', bbox_inches = 'tight')
    #plt.show()
    plt.close('all')

