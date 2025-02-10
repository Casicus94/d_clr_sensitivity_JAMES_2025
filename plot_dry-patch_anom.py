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
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Dry_composites/Hovm/'
user = '/home/netapp-clima/users/acasallas/'

#Run
runs = ['ThoYSU','GCEYSU','MorYSU','MorBL','ThoNDYSU','ThoNDBL']
tipos = ['Dry'] ## Moist, Dry
#Constants
g = 9.81
rho = 1.225
Rd=287.05
Cp=1005.0
Lv=2.5e6
Ls = 2.83e6
Lf = Ls - Lv

z = [40,87.5,141,233,333,436,542,646,752,858,961,
       1060,1213,1421,1623,1825,2031,2241,2450,2657,
       2863,3068,3271,3472,3678,3887,4094,4297,4500,
       4705,4910,5116,5471,5978,6488,6997,7500,8005,
       8513,9013,9509,10009,10503,10988,11468,11949,
       12425,12896,13165,13356,13829,14294,14754,15447,
       16374,17302,18242.6,19190.13,20148.62,21084.01,
       21990.93,22890.83]

for j,run in enumerate(runs):
    print('-----> Starting with '+run)
    #### Water first
    print('Plotting water species')
    fig,(ax)=plt.subplots(nrows=2,ncols=6,figsize=(21,18)) #width,height
    plt.subplots_adjust(hspace=0.2,wspace=0.35,top=0.9,right=0.98,left=0.05)
    for it,tipo in enumerate(tipos):
        print('##### '+tipo+' #####')
        if tipo == 'Dry':
            selec = '_Lagrangian'
        else:
            selec = ''
        if run == 'ThoNDYSU' or run == 'ThoNDBL':
            ds = xr.open_dataset(user+run+'_run/'+tipo+'_patches'+selec+'_'+run+'.nc')
        else: 
            ds = xr.open_dataset(mdir+run+'_run/'+tipo+'_patches'+selec+'_'+run+'.nc')
        ### First the variables
        #Water
        qv = ds.QVAPOR[:,:]
        ql = ds.QCLOUD[:,:]
        qi = ds.QICE[:,:]
        qs = ds.QSNOW[:,:]
        #Press and temperature
        press = ds.P[:,:] + ds.PB[:,:]
        theta = ds.T[:,:] + 300
        temp = theta*np.power(press/1.e5,Rd/Cp)
        tempv=(1.0+0.61*qv)*temp # tempv is virtual T
        thv = tempv*np.power(1.e5/press*100,Rd/Cp)
        thv_mean = thv.mean(dim=['bottom_top'])
        thv_pert = thv - thv_mean
        #Extra diagnostics
        buo = g*(thv_pert/thv_mean)
        mse = Cp*temp + g*np.array(z) + Lv*qv - Lf*qi
        del(theta);del(tempv);del(thv);del(thv_mean);del(thv_pert)
         
        ############ Plot
        to_plot = [qv*1e2,ql*1e4,qi*1e5,qs*1e4,mse/Cp,buo]
        colors = ['Spectral','Blues','Blues','Blues','Spectral_r', 'Purples_r']
        names = ['Qv (1e2 kg kg$_{-1}$)','Ql (1e4 kg kg$_{-1}$)','Qi (1e5 kg kg$_{-1}$)',
                 'Qs (1e4 kg kg$^{-1}$)','MSE/Cp (K)','Buoyancy (m s$^{-2}$)']
        levels = [np.arange(0,1.76,0.25), np.arange(0,0.76,0.15), np.arange(0,0.76,0.15),
                  np.arange(0,0.76,0.15), np.arange(320,381,5), np.arange(-1.5,0.1,0.25)]
        ylabel = ['Days (Dry area)','Days (Random area)']
        fs = 12
        for iplot,plot in enumerate(to_plot):         
            im = ax[it,iplot].contourf(press[0,:]/100,np.linspace(0,20,480),plot,
                                       cmap = colors[iplot],levels = levels[iplot],extend = 'max')
            fig.colorbar(im, ax=ax[it,iplot], shrink = 0.6)
            if iplot == 0 or iplot == 6:
                ax[it,iplot].set_ylabel(ylabel[it], fontsize = fs)
            ax[it,iplot].set_xlim(1000,0)
            ax[it,iplot].set_xlabel('Pressure (hPa)', fontsize = fs)
            ax[it,iplot].set_xticks([1000,750,500,250,0])
            ax[it,iplot].set_title(names[iplot], fontsize = fs)
    plt.savefig(med+run+'/Dry_composites'+selec+'_water_'+run+'.jpg', bbox_inches = 'tight')
    #plt.show()
    plt.close('all')
    del(qv);del(ql);del(qi);del(qs);del(buo);del(mse)
    
    ### Wind   
    print('Plotting wind')
    fig,(ax)=plt.subplots(nrows=2,ncols=3,figsize=(12,18)) #width,height
    plt.subplots_adjust(hspace=0.20,wspace=0.35,top=0.9,right=0.98,left=0.1)
    for it,tipo in enumerate(tipos):
        print('##### '+tipo+' #####')
        if tipo == 'Dry':
            selec = '_Lagrangian'
        else:
            selec = ''
        if run == 'ThoNDYSU' or run == 'ThoNDBL':
            ds = xr.open_dataset(user+run+'_run/'+tipo+'_patches'+selec+'_'+run+'.nc')
        else:
            ds = xr.open_dataset(mdir+run+'_run/'+tipo+'_patches'+selec+'_'+run+'.nc')
        press = ds.P[:,:] + ds.PB[:,:]
        if selec != '_Lagrangian':
            w = ds.W[:,:]
            u = ds.U[:,:]
            v = ds.V[:,:]
            ############ Plot
            to_plot = [u,v,w]
            colors = ['seismic','seismic', 'seismic']
            names = ['Zonal (m s$^{-1}$)','Meridional (m s$^{-1}$)','Vertical velocity (m s$^{-1}$)']
            levels = [np.arange(-0.5,0.51,0.1), np.arange(-0.5,0.51,0.1), np.arange(-0.5,0.51,0.1)] 
            ylabel = ['Days (Dry area)','Days (Random area)']
            for iplot,plot in enumerate(to_plot):
                im = ax[it,iplot].contourf(press[0,:]/100,np.linspace(0,20,480),plot,
                                           cmap = colors[iplot], levels = levels[iplot])
                fig.colorbar(im, ax=ax[it,iplot], shrink = 0.6)
                if iplot == 0 or iplot == 3:
                    ax[it,iplot].set_ylabel(ylabel[it], fontsize = fs)
                ax[it,iplot].set_xlim(1000,0)
                ax[it,iplot].set_xlabel('Pressure (hPa)', fontsize = fs)
                ax[it,iplot].set_xticks([1000,750,500,250,0])
                ax[it,iplot].set_title(names[iplot], fontsize = fs)
            plt.savefig(med+run+'/Dry_composites'+selec+'_wind_'+run+'.jpg', bbox_inches = 'tight')
            #plt.show()
            plt.close('all')
            del(u);del(v);del(w) 
    ### Rad cooling     
    print('Plotting Radiative Cooling')
    fig,(ax)=plt.subplots(nrows=2,ncols=6,figsize=(21,18)) #width,height
    plt.subplots_adjust(hspace=0.2,wspace=0.35,top=0.9,right=0.98,left=0.05)  
    for it,tipo in enumerate(tipos): 
        print('##### '+tipo+' #####') 
        if tipo == 'Dry':
            selec = '_Lagrangian'
        else:
            selec = ''  
        if run == 'ThoNDYSU' or run == 'ThoNDBL':
            ds = xr.open_dataset(user+run+'_run/'+tipo+'_patches'+selec+'_'+run+'.nc')
        else:
            ds = xr.open_dataset(mdir+run+'_run/'+tipo+'_patches'+selec+'_'+run+'.nc')
        press = ds.P[:,:] + ds.PB[:,:]
        #Radiative cooling
        net = ds.RTHRATEN[:,:]
        LW = ds.RTHRATLW[:,:]
        LWC = ds.RTHRATLWC[:,:]
        LWCRE = LW - LWC
        SW = ds.RTHRATSW[:,:]
        SWC = ds.RTHRATSWC[:,:]
        ############ Plot
        to_plot = [net,LW,LWC,LWCRE,SW,SWC]
        colors = ['Spectral_r','Spectral_r','Spectral_r','Reds','Reds','Reds']
        names = ['Net (K day$^{-1}$)','LW all-sky (K day$^{-1}$)',
                 'LW clear-sky (K day$^{-1}$)','LW CRE (K day$^{-1}$)',
                 'SW all-sky (K day$^{-1}$)','SW clear-sky (K day$^{-1}$)']
        levels = [np.arange(-4,4.1,1), np.arange(-3.0,3.01,0.6), np.arange(-3.5,3.51,0.7),
                  np.arange(0,4.01,0.5), np.arange(0,7.51,1.5),np.arange(0,7.51,1.5)]
        ylabel = ['Days (Dry area)','Days (Random area)']
        for iplot,plot in enumerate(to_plot):
            im = ax[it,iplot].contourf(press[0,:]/100,np.linspace(0,20,480),plot*86400, 
                                       cmap = colors[iplot],levels = levels[iplot], extend = 'both')
            fig.colorbar(im, ax=ax[it,iplot], shrink = 0.6)
            if iplot == 0 or iplot == 6:
                ax[it,iplot].set_ylabel(ylabel[it], fontsize = fs)
            ax[it,iplot].set_xlim(1000,0)
            ax[it,iplot].set_xlabel('Pressure (hPa)', fontsize = fs)
            ax[it,iplot].set_xticks([1000,750,500,250,0])
            ax[it,iplot].set_title(names[iplot], fontsize = fs)
    plt.savefig(med+run+'/Dry_composites'+selec+'_Rad_cool_'+run+'.jpg', bbox_inches = 'tight')
    #plt.show()
    #pdb.set_trace()
    plt.close('all')
    del(net);del(LW);del(LWC);del(LWCRE);del(SW);del(SWC)

    #
    ##### Starting the anomalies
    #

    print('########## Starting anomalies!!! ##########')
    if run == 'ThoNDYSU' or run == 'ThoNDBL':
        ds = xr.open_dataset(user+run+'_run/Domain_verticals_'+run+'.nc')
        ds_dry = xr.open_dataset(user+run+'_run/Dry_patches'+selec+'_'+run+'.nc')
    else:
        ds = xr.open_dataset(mdir+run+'_run/Domain_verticals_'+run+'.nc')
        ds_dry = xr.open_dataset(mdir+run+'_run/Dry_patches'+selec+'_'+run+'.nc')
    #ds_moi = xr.open_dataset(mdir+run+'_run/Moist_patches_'+run+'.nc')
    ### First the variables
    #Water
    print('Water species')
    qv_d = ds_dry.QVAPOR[:,:]
    ql_d = ds_dry.QCLOUD[:,:]
    qi_d = ds_dry.QICE[:,:]
    qs_d = ds_dry.QSNOW[:,:]
    qr_d = ds_dry.QRAIN[:,:]
    #Domain
    qv_t = ds.QVAPOR[0:480,:]
    ql_t = ds.QCLOUD[0:480,:]
    qi_t = ds.QICE[0:480,:]
    qs_t = ds.QSNOW[0:480,:]
    qr_t = ds.QRAIN[0:480,:]
    #Anomalies
    #Dry - domain
    qv_ta = qv_d - np.array(qv_t[:,:,0,0])
    ql_ta = ql_d - np.array(ql_t[:,:,0,0])
    qi_ta = qi_d - np.array(qi_t[:,:,0,0])
    qs_ta = qs_d - np.array(qs_t[:,:,0,0])
    qr_ta = qr_d - np.array(qr_t[:,:,0,0])
    #Press and temperature
    #Dry
    mse_d = ds_dry.MSE[0:480,:]
    tp_d = ds_dry.Temp_anom[0:480,:]
    thvp_d = ds_dry.THv_anom[0:480,:]
    qvp_d = ds_dry.Qv_anom[0:480,:] 
    buo_d = ds_dry.Buoyancy[0:480,:]
    #Domain
    if run == 'ThoNDYSU' or run == 'ThoNDBL':
        ds11 = xr.open_dataset(user+run+'_run/Domain_buo_comp_'+run+'.nc')
    else:
        ds11 = xr.open_dataset(mdir+run+'_run/Domain_buo_comp_'+run+'.nc')
    mse_t = ds11.MSE[0:480,:]
    tp_t = ds11.Temp_anom[0:480,:]
    thvp_t = ds11.THv_anom[0:480,:]
    qvp_t = ds11.Qv_anom[0:480,:]
    buo_t = ds11.Buoyancy[0:480,:]
    del(ds11)
    #Anomaly
    #MSE
    mse_ta = mse_d - np.array(mse_t)
    press_ta = (ds_dry.P[0:480,:] + ds_dry.PB[0:480,:]) - np.array(ds.P[0:480,:,0,0] + ds.PB[0:480,:,0,0])
    del(mse_d)
    del(mse_t)
    #Buoyancy
    qvp_ta = qvp_d - np.array(qvp_t)
    tp_ta = tp_d - np.array(tp_t)
    thvp_ta = thvp_d - np.array(thvp_t)
    buo_ta = buo_d - np.array(buo_t)
    del(buo_d);del(tp_d);del(qvp_d);del(qvp_t)
    del(buo_t);del(tp_t);del(thvp_d);del(thvp_t)
    
    #
    ##### Plot Anomalies #####
    #
    
    ########################
    #     Dry - Domain     #
    ########################
    to_plot = [qv_ta*1e3,ql_ta*1e4,qi_ta*1e5,qs_ta*1e4,qr_ta*1e4,mse_ta/Cp,press_ta]
    colors = ['seismic_r','seismic_r','seismic_r','seismic_r','seismic_r','seismic_r','seismic']
    names = ['Qv (1e3 g kg$^{-1}$)','Ql (1e4 g kg$^{-1}$)','Qi (1e5 g kg$^{-1}$)',
             'Qs (1e4 g kg$^{-1}$)','Qr (1e4 g kg$^{-1}$)','MSE/Cp (K)','Pressure (hPa)']
    levels = [np.arange(-1,1.1,0.2),np.arange(-0.03,0.031,0.005),np.arange(-0.1,0.101,0.01),np.arange(-0.07,0.0701,0.005),np.arange(-0.015,0.016,0.0025),np.arange(-2,2.1,0.25),np.arange(-8,8.1,1)]

    fig,(ax)=plt.subplots(nrows=1,ncols=7,figsize=(22,8))
    plt.subplots_adjust(hspace=0.0,wspace=0.35,top=0.9,right=0.98,left=0.05)
    for iplot,plot in enumerate(to_plot):
        im = ax[iplot].contourf(press[0,:]/100,np.linspace(0,20,480),plot,
                                   cmap = colors[iplot],levels = levels[iplot],extend = 'both')
        fig.colorbar(im, ax=ax[iplot], shrink = 0.6)
        ax[iplot].set_xlim(1000,0)
        ax[iplot].set_xlabel('Pressure (hPa)', fontsize = fs)
        ax[iplot].set_xticks([1000,750,500,250,0])
        ax[iplot].set_title(names[iplot], fontsize = fs)
    plt.savefig(med+run+'/Anom_composites'+selec+'_hum_dom_'+run+'.jpg', bbox_inches = 'tight')
    #plt.show()
    plt.close('all')
    
    ### Buoyancy
    to_plot = [buo_ta*1e3,tp_ta*1e3,thvp_ta*1e4,qvp_ta*1e3]
    colors = ['seismic','seismic','seismic','seismic']
    names = ['Buoyancy (1e3 m s$^{-2}$)',u"T$^{'} / \overline{T}$ (1e3)",
             u"\u03F4$_{v}^{'} / \overline{\u03F4_{v}}$ (1e4)",
             "(0.61q$_{v}^{'}) / (1+0.61\overline{q_{v}}$) (1e3)"]
    levels = [np.arange(-2.,2.01,0.25), np.arange(-1.,1.01,0.1), 
              np.arange(-2.,2.01,0.25), np.arange(-1.,1.01,0.1)]

    fig,(ax)=plt.subplots(nrows=1,ncols=4,figsize=(16,8))
    plt.subplots_adjust(hspace=0.0,wspace=0.35,top=0.9,right=0.98,left=0.05)
    for iplot,plot in enumerate(to_plot):
        im = ax[iplot].contourf(press[0,:]/100,np.linspace(0,20,480),plot,
                                cmap = colors[iplot],levels = levels[iplot],extend = 'both')
        fig.colorbar(im, ax=ax[iplot], shrink = 0.6)
        ax[iplot].set_xlim(1000,0)
        ax[iplot].set_xlabel('Pressure (hPa)', fontsize = fs)
        ax[iplot].set_xticks([1000,750,500,250,0])
        ax[iplot].set_title(names[iplot], fontsize = fs)
    plt.savefig(med+run+'/Anom_composites'+selec+'_buo_dom_'+run+'.jpg', bbox_inches = 'tight')
    #plt.show()
    plt.close('all')
    
    if selec != '_Lagrangian':
        ##### Wind
        print('Wind')
        #Dry
        w_d = ds_dry.W[:,:]
        u_d = ds_dry.U[:,:]
        v_d = ds_dry.V[:,:]
        #Random
        w_m = ds_moi.W[:,:]
        u_m = ds_moi.U[:,:]
        v_m = ds_moi.V[:,:]
        #Domain
        w_t = ds.W[0:480,0:62]
        u_t = ds.U[0:480,:]
        v_t = ds.V[0:480,:]
        #Anomaly
        w_a = w_d - w_m
        u_a = u_d - u_m
        v_a = v_d - v_m
        w_ta = w_d - np.array(w_t[:,:,0,0])
        u_ta = u_d - np.array(u_t[:,:,0,0])
        v_ta = v_d - np.array(v_t[:,:,0,0])
        del(w_d);del(u_d);del(v_d)
        del(w_m);del(u_m);del(v_m)
        del(w_t);del(u_t);del(v_t)
    
        #
        ##### Plot Anomalies #####
        #

        ########################
        #     Dry - Random     #
        ######################## 
        to_plot = [u_a,v_a,w_a]
        colors = ['seismic','seismic', 'seismic']
        names = ['Zonal (m s$^{-1}$)','Meridional (m s$^{-1}$)','Vertical velocity (m s$^{-1}$)']
        levels = [np.arange(-1,1.1,0.2), np.arange(-1,1.1,0.2), np.arange(-1,1.1,0.2)]
        fig,(ax)=plt.subplots(nrows=1,ncols=3,figsize=(12,8))
        plt.subplots_adjust(hspace=0.0,wspace=0.35,top=0.9,right=0.98,left=0.05)
        for iplot,plot in enumerate(to_plot):
            im = ax[iplot].contourf(press[0,:]/100,np.linspace(0,20,480),plot,
                                   cmap = colors[iplot],levels = levels[iplot])
            fig.colorbar(im, ax=ax[iplot], shrink = 0.6)
            ax[iplot].set_xlim(1000,0)
            ax[iplot].set_xlabel('Pressure (hPa)', fontsize = fs)
            ax[iplot].set_xticks([1000,750,500,250,0])
            ax[iplot].set_title(names[iplot], fontsize = fs)
        plt.savefig(med+run+'/Anom_composites'+selec+'_wind_'+run+'.jpg', bbox_inches = 'tight')
        #plt.show()
        plt.close('all')
   
        ########################
        #     Dry - Domain     #
        ########################
        to_plot = [u_ta,v_ta,w_ta]
        colors = ['seismic','seismic', 'seismic']
        names = ['Zonal (m s$^{-1}$)','Meridional (m s$^{-1}$)','Vertical velocity (m s$^{-1}$)']
        levels = [np.arange(-1,1.1,0.2), np.arange(-1,1.1,0.2), np.arange(-1,1.1,0.2)]
        fig,(ax)=plt.subplots(nrows=1,ncols=3,figsize=(12,8))
        plt.subplots_adjust(hspace=0.0,wspace=0.35,top=0.9,right=0.98,left=0.05)
        for iplot,plot in enumerate(to_plot):
            im = ax[iplot].contourf(press[0,:]/100,np.linspace(0,20,480),plot,
                                   cmap = colors[iplot],levels = levels[iplot])
            fig.colorbar(im, ax=ax[iplot], shrink = 0.6)
            ax[iplot].set_xlim(1000,0)
            ax[iplot].set_xlabel('Pressure (hPa)', fontsize = fs)
            ax[iplot].set_xticks([1000,750,500,250,0])
            ax[iplot].set_title(names[iplot], fontsize = fs)
        plt.savefig(med+run+'/Anom_composites'+selec+'_wind_dom_'+run+'.jpg', bbox_inches = 'tight')
        #plt.show()
        plt.close('all')
    elif selec == '_Lagrangian':
        w_d = ds_dry.W[:,:]
        w_t = ds.W[0:480,0:62]
        w_ta = w_d - np.array(w_t[:,:,0,0])
        to_plot = [w_ta*10]
        colors = ['seismic']
        names = ['Vertical velocity (m s$^{-1}$)']
        levels = [np.arange(-0.25,0.26,0.05)]
        fig,(ax)=plt.subplots(nrows=1,ncols=1,figsize=(5,8))
        for iplot,plot in enumerate(to_plot):
            im = ax.contourf(press[0,:]/100,np.linspace(0,20,480),plot,
                                   cmap = colors[iplot],levels = levels[iplot])
            fig.colorbar(im, ax=ax, shrink = 0.6)
            ax.set_xlim(1000,0)
            ax.set_ylabel('Days')
            ax.set_xlabel('Pressure (hPa)', fontsize = fs)
            ax.set_xticks([1000,750,500,250,0])
            ax.set_title(names[iplot], fontsize = fs)
        plt.savefig(med+run+'/Anom_composites'+selec+'_wind_dom_'+run+'.jpg', bbox_inches = 'tight')
        #plt.show()
        plt.close('all')
  
    ##### Rad Cooling
    print('Rad Cooling')
    #Dry
    net_d = ds_dry.RTHRATEN[:,:]
    LW_d = ds_dry.RTHRATLW[:,:]
    LWC_d = ds_dry.RTHRATLWC[:,:]
    LWCRE_d = LW_d - LWC_d
    SW_d = ds_dry.RTHRATSW[:,:]
    SWC_d = ds_dry.RTHRATSWC[:,:]
    #Random
    #net_m = ds_moi.RTHRATEN[:,:]
    #LW_m = ds_moi.RTHRATLW[:,:]
    #LWC_m = ds_moi.RTHRATLWC[:,:]
    #LWCRE_m = LW_m - LWC_m
    #SW_m = ds_moi.RTHRATSW[:,:]
    #SWC_m = ds_moi.RTHRATSWC[:,:]
    #Domain
    net_t = ds.RTHRATEN[0:480,:]
    LW_t = ds.RTHRATLW[0:480,:]
    LWC_t = ds.RTHRATLWC[0:480,:]
    LWCRE_t = LW_t - LWC_t
    SW_t = ds.RTHRATSW[0:480,:]
    SWC_t = ds.RTHRATSWC[0:480,:]
    #Anomalies
    #net_a = net_d - net_m
    #LW_a = LW_d - LW_m    
    #LWC_a = LWC_d - LWC_m
    #LWCRE_a = LWCRE_d - LWCRE_m
    #SW_a = SW_d - SW_m
    #SWC_a = SWC_d - SW_m
    net_ta = net_d - np.array(net_t[:,:,0,0])
    LW_ta = LW_d - np.array(LW_t[:,:,0,0])
    LWC_ta = LWC_d - np.array(LWC_t[:,:,0,0])
    LWCRE_ta = LWCRE_d - np.array(LWCRE_t[:,:,0,0])
    SW_ta = SW_d - np.array(SW_t[:,:,0,0])
    SWC_ta = SWC_d - np.array(SW_t[:,:,0,0])
    del(net_d);del(LW_d);del(LWC_d);del(LWCRE_d);del(SW_d);del(SWC_d)
    #del(net_m);del(LW_m);del(LWC_m);del(LWCRE_m);del(SW_m);del(SWC_m)
    del(net_t);del(LW_t);del(LWC_t);del(LWCRE_t);del(SW_t);del(SWC_t)
    
    #
    ##### Plot Anomalies #####
    #

    ########################
    #     Dry - Random     #
    ######################## 
    #fig,(ax)=plt.subplots(nrows=1,ncols=6,figsize=(18,8))
    #plt.subplots_adjust(hspace=0.0,wspace=0.35,top=0.9,right=0.98,left=0.05)
    #to_plot = [net_a,LW_a,LWC_a,LWCRE_a,SW_a,SWC_a]
    #colors = ['seismic','seismic','seismic','seismic','seismic','seismic']
    #names = ['Net (K day$^{-1}$)','LW all-sky (K day$^{-1}$)',
    #         'LW clear-sky (K day$^{-1}$)','LW CRE (K day$^{-1}$)',
    #         'SW all-sky (K day$^{-1}$)','SW clear-sky (K day$^{-1}$)']
    #levels = [np.arange(-3,3.1,0.5), np.arange(-3.0,3.1,0.5), np.arange(-1.0,1.01,0.2),
    #          np.arange(-2.5,2.51,0.25), np.arange(-2,2.01,0.2),np.arange(-2.0,2.01,0.2)]
    #for iplot,plot in enumerate(to_plot):
    #    im = ax[iplot].contourf(press[0,:]/100,np.linspace(0,20,480),plot*86400,
    #                               cmap = colors[iplot],levels = levels[iplot],extend = 'both')
    #    fig.colorbar(im, ax=ax[iplot], shrink = 0.6)
    #    ax[iplot].set_xlim(1000,0)
    #    ax[iplot].set_xlabel('Pressure (hPa)', fontsize = fs)
    #    ax[iplot].set_xticks([1000,750,500,250,0])
    #    ax[iplot].set_title(names[iplot], fontsize = fs)
    #plt.savefig(med+run+'/Anom_composites'+selec+'_Rad_Cooling_'+run+'.jpg', bbox_inches = 'tight')
    #plt.show() 
    #plt.close('all')
    
    ########################
    #     Dry - Domain     #
    ########################
    fig,(ax)=plt.subplots(nrows=1,ncols=6,figsize=(18,8))
    plt.subplots_adjust(hspace=0.0,wspace=0.35,top=0.9,right=0.98,left=0.05)
    to_plot = [net_ta,LW_ta,LWC_ta,LWCRE_ta,SW_ta,SWC_ta]
    colors = ['seismic','seismic','seismic','seismic','seismic','seismic']
    names = ['Net (K day$^{-1}$)','LW all-sky (K day$^{-1}$)',
             'LW clear-sky (K day$^{-1}$)','LW CRE (K day$^{-1}$)',
             'SW all-sky (K day$^{-1}$)','SW clear-sky (K day$^{-1}$)']
    levels = [np.arange(-1,1.1,0.25), np.arange(-1.0,1.1,0.25), np.arange(-0.5,0.51,0.1),
              np.arange(-1.5,1.51,0.25), np.arange(-1,1.01,0.2),np.arange(-1.0,1.01,0.2)]
    for iplot,plot in enumerate(to_plot):
        im = ax[iplot].contourf(press[0,:]/100,np.linspace(0,20,480),plot*86400,
                                   cmap = colors[iplot],levels = levels[iplot],extend = 'both')
        fig.colorbar(im, ax=ax[iplot], shrink = 0.6)
        ax[iplot].set_xlim(1000,0)
        ax[iplot].set_xlabel('Pressure (hPa)', fontsize = fs)
        ax[iplot].set_xticks([1000,750,500,250,0])
        ax[iplot].set_title(names[iplot], fontsize = fs)
    plt.savefig(med+run+'/Anom_composites'+selec+'_Rad_Cool_dom_'+run+'.jpg', bbox_inches = 'tight')
    #plt.show()
    plt.close('all')
print('Complete!!!') 
