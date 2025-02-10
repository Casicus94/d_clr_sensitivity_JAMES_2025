import numpy as np
import pandas as pd
from netCDF4 import Dataset,num2date,date2num
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.gridspec as gridspec
from pylab import *
import Casicus as casi
import pdb

#Run
runs = ['ThoYSU','GCEYSU','MorYSU','MorBL','ThoNDBL','ThoNDYSU']
#runs = ['ThoNDBL']

#Some constants and dimensions
times = 480
date_ini = datetime.datetime(2007,6,1)
elv = 62
elev = np.arange(0,elv)
#Path
mdir = '/home/tompkins-archive/acasallas/'
user = '/home/netapp-clima/users/acasallas/'
 
### All variables
varis = ['Buoyancy','THv_anom', 'Temp_anom', 'MSE', 'Qv_anom']

sta = 0
end = 480

z = [40,87.5,141,233,333,436,542,646,752,858,961,
       1060,1213,1421,1623,1825,2031,2241,2450,2657,
       2863,3068,3271,3472,3678,3887,4094,4297,4500,
       4705,4910,5116,5471,5978,6488,6997,7500,8005,
       8513,9013,9509,10009,10503,10988,11468,11949,
       12425,12896,13165,13356,13829,14294,14754,15447,
       16374,17302,18242.6,19190.13,20148.62,21084.01,
       21990.93,22890.83]
#Constants
Rd=287.05
Cp=1005.0
Lv=2.5e6
Ls = 2.83e6
Lf = Ls - Lv
g = 9.81

for j,run in enumerate(runs):
    print('-----> Starting with '+run)
    if run == 'ThoNDYSU' or run == 'ThoNDBL':
        ds = Dataset(user+run+'_run/Variables_'+run+'_all.nc')
    else:
        ds = Dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
    unout = ds.variables['XTIME'].description
    dates = np.array([])
    #Create dates 
    for t in range(times):
        dates =np.append(dates, date_ini+datetime.timedelta(minutes=int(ds.variables['XTIME'][t])))
    del(ds)
    if run == 'ThoNDYSU' or run == 'ThoNDBL':
        ds = xr.open_dataset(user+run+'_run/Variables_'+run+'_all.nc')
    else:
        ds = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
    ### Add more complicated variables!
    print('Reading variables from file')
    print('##### Qv #####')
    qv = ds.QVAPOR[sta:end,:,:,:]
    print('##### Qi #####')
    qi = ds.QICE[sta:end,:,:,:]
    print('##### Press #####')
    press = ds.P[sta:end,:,:,:] + ds.PB[sta:end,:,:,:]
    print('##### Theta #####')
    theta = ds.T[sta:end,:,:,:] + 300
    temp = theta*np.power(press/1.e5,Rd/Cp)
    del(theta)
    print('Starting with Diagnostics')
    data_dry = []
    tempv=(1.0+0.61*qv)*temp # tempv is virtual T
    thv = tempv*np.power(1.e5/press,Rd/Cp)
    del(tempv);del(press)
    thv_mean = thv.mean(dim=['south_north','west_east'])
    thv_pert = thv - thv_mean
    del(thv)
    print('##### Buoyancy #####')
    buo = g*(thv_pert/thv_mean)
    data_dry.append(buo.mean(dim=['south_north','west_east']))
    del(buo)
    print('##### Buo Components #####')
    thvanom = thv_pert/thv_mean
    data_dry.append(thvanom.mean(dim=['south_north','west_east']))
    del(thv_pert);del(thv_mean)
    del(thvanom)
    tanom = (temp-temp.mean(dim=['south_north','west_east']))/(temp.mean(dim=['south_north','west_east']))
    data_dry.append(tanom.mean(dim=['south_north','west_east']))
    del(tanom)
    print('##### MSE #####')
    mse = Cp*temp + g*np.array(z)[None,:,None,None] + Lv*qv - Lf*qi
    data_dry.append(mse.mean(dim=['south_north','west_east']))
    del(temp);del(qi);del(mse)
    print('##### Qv Anom #####')
    qvanom = (0.61*(qv-qv.mean(dim=['south_north','west_east'])))/(1+0.61*qv.mean(dim=['south_north','west_east']))
    data_dry.append(qvanom.mean(dim=['south_north','west_east']))
    del(ds);del(qv)

    ###Creating files
    print('Adding to file')
    ### File
    if run == 'ThoNDYSU' or run == 'ThoNDBL':
        ncout = Dataset(user+run+'_run/Domain_buo_comp_'+run+'.nc','w','NETCDF4')
    else:
        ncout = Dataset(mdir+run+'_run/Domain_buo_comp_'+run+'.nc','w','NETCDF4')
    ncout.createDimension('bottom_top',elv)
    ncout.createDimension('Time',times);
    elvvar = ncout.createVariable('bottom_top', 'float32', ('bottom_top'));elvvar[:] = elev;
    timevar = ncout.createVariable('Time','float32',('Time'));timevar.setncattr('units',unout);timevar[:]=date2num(dates);
    for i,v in enumerate(varis):
        print('Including: '+v)
        myvar = ncout.createVariable(v,'float32',('Time','bottom_top'));myvar[:] = data_dry[i][:,0:62];
    ncout.close();


