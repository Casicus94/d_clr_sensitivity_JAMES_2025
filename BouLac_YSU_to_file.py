import numpy as np
import pandas as pd
import xarray as xr
import pdb
import scipy.spatial as spatial
from netCDF4 import Dataset
import Casicus as casi
from wrf import tvirtual
from Casicus import CPs_object, labels_cp, cp_geometry

rho = 1.225
Rd=287.05
Cp=1005.0
mdir = '/home/tompkins-archive/acasallas/'
scra = '/home/netapp-clima/scratch/acasallas/'

nx = 256
ny = 256
elv = 18

### pcentile for edge detection in gradient 
pcen=60

def dictionary_list(runs):
    dict_list = {}
    for run in runs:
        dict_list[run] = []
    return(dict_list)

def dictionary_array(runs):
    dict_list = {}
    for run in runs:
        dict_list[run] = np.array([])
    return(dict_list)

def d_clr_num_cores(wt,nx,ny):
    ### Calculate d_clr and the number of cores
    nt = nx*ny
    cnv_coords = np.argwhere(wt > 1)
    numm_conv = np.shape(cnv_coords)[0]
    allidx=np.argwhere(wt<0)
    for xoff in [0,nx,-nx]:
        for yoff in [0,-ny,ny]:
            if xoff==0 and yoff==0:
                j9=cnv_coords.copy()
            else:
                jo=cnv_coords.copy()
                jo[:,0]+=xoff
                jo[:,1]+=yoff
                j9=np.vstack((j9,jo))
    f_cnvdst = np.array([])
    tree=spatial.cKDTree(j9)
    d2a,minidx=tree.query(allidx)
    f_cnvdst = np.hstack((f_cnvdst, d2a)) 
    if len(f_cnvdst) == 0:
        dist_all=0
    else:
        dist_all=np.nanmax(f_cnvdst)
    cnvdst = np.array([])
    for k in range(len(cnv_coords)):
        a = np.delete(j9, list(range(k, j9.shape[0], len(cnv_coords))), axis=0)
        tree=spatial.cKDTree(a)
        d2c,minidx=tree.query(cnv_coords[k,:])
        cnvdst = np.hstack((cnvdst, d2c))
    if len(cnvdst) == 0:
        dist_conv=0
    else:
        dist_conv=np.nanmax(cnvdst)
    return(dist_all,dist_conv,numm_conv)

def cp_intensity(ds,t):
    press = (ds.P[t,:,:,:] + ds.PB[t,:,:,:])
    T = (ds.T[t,:,:,:]+300)*np.power(press*100/1.e5,Rd/Cp)
    Tv = tvirtual(T,ds.QVAPOR[t,:,:,:])
    thv = Tv*np.power(1.e5/press*100,Rd/Cp)
    thv_mean = thv.mean(dim=['south_north','west_east'])
    thv_pert = thv - thv_mean
    B = 9.81*(thv_pert/thv_mean)
    del(Tv); del(thv); del(thv_mean); del(thv_pert); del(press); del(T)
    colp_b = B[0:elv,:,:].where(B[0:elv,:,:]<-0.005)
    inte = 2*casi.verint_xarray(ds,colp_b*(-1),t,elv)
    return(inte)

runs = ['ThoBLY','ThoTBLY','GCEBLY','GCETBLY']*2

dist = dictionary_list(runs)
numm = dictionary_list(runs)
gust_conv = dictionary_list(runs)
cp_int = dictionary_array(runs)
pblh = dictionary_array(runs)

#### Here we calculate d_clr and the number of cores
try:
    df = pd.read_csv(mdir+'Data_to_plot/BLY_dclr.csv')
    df1 = pd.read_csv(mdir+'Data_to_plot/BLY_core_num.csv')
except:
    for j,run in enumerate(runs):
        print('----------------> Calculating d_clr and core counts')
        if j < len(runs)/2:
            print('-----> Starting with '+run[0:-1])
            ds = Dataset(mdir+run[0:-1]+'_run/Variables_'+run[0:-1]+'_all.nc')
        else:
            print('-----> Starting with '+run)
            ds = Dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
        count = 0 
        for i in np.arange(0,237):
            print('Step: '+str(i)+' of 240')
            if j < len(runs)/2:
                wt = ds.variables['W'][-240+i,19,:,:]
            else:
                wt = ds.variables['W'][i,19,:,:]
            d_all,d_conv,n_conv = d_clr_num_cores(wt,nx,ny)
            dist[run].append(d_all)
            numm[run].append(n_conv)
    data_dist = pd.DataFrame.from_dict(dist)
    data_dist.to_csv(mdir+'Data_to_plot/BLY_dclr.csv', index = False)
    data_numm = pd.DataFrame.from_dict(numm)
    data_numm.to_csv(mdir+'Data_to_plot/BLY_core_num.csv', index = False)

### Here we calculate the PBLH
try:
    df_pbl = pd.read_csv(mdir+'Data_to_plot/BLY_PBLH.csv')
except: 
    for j,run in enumerate(runs):
        print('----------------> Calculating PBLH')
        if j < len(runs)/2:
            print('-----> Starting with '+run[0:-1])
            ds = xr.open_dataset(mdir+run[0:-1]+'_run/Variables_'+run[0:-1]+'_all.nc')
            time = np.arange(-240,0)
        else:
            print('-----> Starting with '+run)
            ds = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
            time = np.arange(0,237)
        for t in time:
            print('Step: '+str(t)+' of 240')
            pbl = ds.PBLH[t,:,:]
            pblh[run] = np.append(pblh[run], pbl.mean())
    data_pblh = pd.DataFrame.from_dict(pblh)
    data_pblh.to_csv(mdir+'Data_to_plot/BLY_PBLH.csv', index = False)

### Here we calculate the CP intensity
try:
    df_cp = pd.read_csv(mdir+'Data_to_plot/BLY_cp_int.csv')
except:
    for j,run in enumerate(runs): 
        print('----------------> Calculating CP intensity')
        if j < len(runs)/2:
            print('-----> Starting with '+run[0:-1])
            ds = xr.open_dataset(mdir+run[0:-1]+'_run/Variables_'+run[0:-1]+'_all.nc')
            time = np.arange(-240,0)
        else:
            print('-----> Starting with '+run)
            ds = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
            time = np.arange(0,237)
        for t in time:
            print('Step: '+str(t)+' of 240')
            inte = cp_intensity(ds,t)
            cp_int[run] = np.append(cp_int[run], inte.mean())
    data_cp_int = pd.DataFrame.from_dict(cp_int)
    data_cp_int.to_csv(mdir+'Data_to_plot/BLY_cp_int.csv', index = False)

#### Here we calculate the gust fronts convergence
try:
    df_ww = pd.read_csv(mdir+'Data_to_plot/BLY_cp_ww.csv')
except:
    for j,run in enumerate(runs):
        print('----------------> Calculating CP convergence')
        if j < len(runs)/2:
            print('-----> Starting with '+run[0:-1])
            ds = xr.open_dataset(mdir+run[0:-1]+'_run/Variables_'+run[0:-1]+'_all.nc')
            time = np.arange(-240,0)
        else:
            print('-----> Starting with '+run)
            ds = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
            time = np.arange(0,238)
        t2=ds.T2[time[0]:time[-1],:,:]
        t2=(t2.differentiate('south_north')*t2.differentiate('south_north')+t2.differentiate('west_east')*t2.differentiate('west_east'))**0.5
        ww = ds.W[time[0]:time[-1],1,:,:]
        #Select coldpools
        print('Selecting coldpools')
        coldpools, distance = casi.CPs_object(t2,pcen)
        for t in np.arange(0,237):
            print('Step '+str(t)+' from 240')
            #Coldpools labels
            mat = casi.labels_cp(distance[t,:,:],coldpools[t,:,:])
            ww_tmp = ww[t,:,:].where(mat>0)
            gust_conv[run].append(ww_tmp.mean().values*24)
    data_cp_ww = pd.DataFrame.from_dict(gust_conv)
    data_cp_ww.to_csv(mdir+'Data_to_plot/BLY_cp_ww.csv', index = False)

