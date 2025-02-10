import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pdb
import matplotlib.gridspec as gridspec
from pylab import *

scra = '/home/netapp-clima/scratch/acasallas/'
mdir = '/home/tompkins-archive/acasallas/'

sta = 10 
end = 15
nx=256
ny=256
nz=62
npts=nx*ny # hardwired for moment
hovchunk=1024  #1024 # averaging length for chunking, power of 2, larger=smoother.
nhov=int(npts/hovchunk)
hovx=100*(np.arange(nhov)+0.5)/nhov

runs = ['ThoYSU']#,'ThoBL','GCEYSU','GCEBL','MorYSU','MorBL','WSMYSU','WSMBL']
fluxes = ['lw','lwc','lwcl','sw','swc','swcl']
sflevs = [-36,-24,-12,-8,-4,0,4,8,12,24,36] 
title = ['(a) Q$_{rad}$ LW net','(b) Q$_{rad}$ LW Clear-sky','(c) Q$_{rad}$ LW CRE','(d) Q$_{rad}$ SW net','(e) Q$_{rad}$ SW Clear-sky','(f) Q$_{rad}$ SW CRE']

for run in runs:
    ds = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
    znu=  ds.ZNU[end*24,:]*ds.PSFC[end*24,0,0]/100.
    fig = plt.figure(figsize=(10,11))
    gs = GridSpec(3,3,left = 0.09, right = 0.98, hspace=0.4, wspace=0.1, top = 0.92, bottom = 0.08, height_ratios = [1,1,0.1])
    for i,flux in enumerate(fluxes): 
        df = pd.read_csv(mdir+'Data_to_plot/circulation/'+run+'_'+flux+'_'+str(sta)+'_'+str(end)+'.csv')
        df = df.drop('Unnamed: 0', axis = 1)
        cdf = pd.read_csv(mdir+'Data_to_plot/circulation/'+run+'_c'+flux+'_'+str(sta)+'_'+str(end)+'.csv')
        cdf = cdf.drop('Unnamed: 0', axis = 1)
        pdb.set_trace()
        ax = plt.subplot(gs[i]) 
        im = plt.contourf(hovx, np.array(znu), np.array(df)*86400, cmap = 'bwr', levels = np.arange(-0.4,0.41,0.04), extend = 'both')
        #cssf = plt.contour(hovx,np.array(znu),np.array(cdf)*3600*100, colors = 'k', levels = np.array(sflevs), linewidths = 0.7) 
        #plt.clabel(cssf, fontsize=8, inline=1,fmt = '%1.2f', levels = sflevs)
        plt.ylim(1000,0)
        plt.title(title[i], loc = 'left')
        if i == 0 or i == 3:
            plt.ylabel('Pressure (hPa)')
        else:
            ax.yaxis.set_major_formatter(NullFormatter())
        plt.xlabel('TCWV %-tile')
    ax = plt.subplot(gs[2,0:])
    cbar = plt.colorbar(im, cax = ax, orientation = 'horizontal')
    cbar.set_label('K day$^{-1}$')
    plt.savefig(scra+'circulation/plots/'+run+'_NS_rad_cooling_'+str(sta)+'_'+str(end)+'.jpg', bbox_inches = 'tight', dpi = 300)
    plt.savefig(scra+'circulation/plots/'+run+'_NS_rad_cooling_'+str(sta)+'_'+str(end)+'.pdf', bbox_inches = 'tight', dpi = 300)
    #plt.show()
    plt.close()

for run in runs:
    ds = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
    znu=  ds.ZNU[end*24,:]*ds.PSFC[end*24,0,0]/100.
    fig = plt.figure(figsize=(6,6))
    gs = GridSpec(2,1,left = 0.09, right = 0.98, hspace=0.4, wspace=0.1, top = 0.92, bottom = 0.08, height_ratios = [1,0.1])
    df = pd.read_csv(mdir+'Data_to_plot/circulation/'+run+'_net_'+str(sta)+'_'+str(end)+'.csv')
    df = df.drop('Unnamed: 0', axis = 1)
    cdf = pd.read_csv(mdir+'Data_to_plot/circulation/'+run+'_cnet_'+str(sta)+'_'+str(end)+'.csv')
    cdf = cdf.drop('Unnamed: 0', axis = 1)
    ax = plt.subplot(gs[0])
    im = plt.contourf(hovx, np.array(znu), np.array(df)*86400, cmap = 'bwr', levels = np.arange(-0.4,0.41,0.04), extend = 'both')
    #cssf = plt.contour(hovx,np.array(znu),np.array(cdf)*3600*100, colors = 'k', levels = np.array(sflevs), linewidths = 0.7)
    #plt.clabel(cssf, fontsize=8, inline=1,fmt = '%1.2f', levels = sflevs)
    plt.ylim(1000,0)
    plt.title('Q$_{rad}$ Net', loc = 'left')
    plt.xlabel('TCWV %-tile')
    plt.ylabel('Pressure (hPa)')
    ax = plt.subplot(gs[1])
    cbar = plt.colorbar(im, cax = ax, orientation = 'horizontal')
    cbar.set_label('K day$^{-1}$')
    plt.savefig(scra+'circulation/plots/'+run+'_NS_rad_cooling_net_'+str(sta)+'_'+str(end)+'.jpg', bbox_inches = 'tight', dpi = 300)
    plt.savefig(scra+'circulation/plots/'+run+'_NS_rad_cooling_net_'+str(sta)+'_'+str(end)+'.pdf', bbox_inches = 'tight', dpi = 300)
    plt.close() 

