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
runs = ['ThoYSU','GCEYSU','MorYSU','MorBL']

#Some constants and dimensions
times = 479
date_ini = datetime.datetime(2007,6,1)
elv = 62
elev = np.arange(0,elv)
#Path
mdir = '/home/tompkins-archive/acasallas/'
user = '/home/netapp-clima/users/acasallas/MTD/MTD_output/DryP_QVAPORTH44p5/'
user1 = '/home/netapp-clima/users/acasallas/'

### Variables from file
varis = ['QVAPOR','QCLOUD','QICE','QSNOW','P','PB','PH','PHB','RTHRATEN',
         'RTHRATLW','RTHRATLWC','RTHRATSW','RTHRATSWC']
### All variables
varis1 = ['QVAPOR','QCLOUD','QICE','QSNOW','P','PB','PH','PHB','RTHRATEN',
          'RTHRATLW','RTHRATLWC','RTHRATSW','RTHRATSWC','Theta']

tiempo = [224,150,298,367,224,224]
sta = 0
end = 479

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
ny,nx = (256,256)
lon = np.arange(0,nx,1)
lat = np.arange(0,ny,1)

for j,run in enumerate(runs):
    print('-----> Starting with '+run)
    if run == 'ThoNDYSU' or run == 'ThoNDBL':
        ds = Dataset(user1+run+'_run/Variables_'+run+'_all.nc')
    else:
        ds = Dataset(mdir+run+'_run/Variables_'+run+'_all.nc') 
    unout = ds.variables['XTIME'].description
    dates = np.array([])
    #Create dates 
    for t in range(times):
        dates =np.append(dates, date_ini+datetime.timedelta(minutes=int(ds.variables['XTIME'][t])))
    del(ds)
    if run == 'ThoNDYSU' or run == 'ThoNDBL':
        ds = xr.open_dataset(user1+run+'_run/Variables_'+run+'_all.nc')
    else:
        ds = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
    ds1 = Dataset(user+'mtd_'+run+'_TCWV44p5_minvol10___20091201_003000V_obj.nc')
    lon,lat = np.meshgrid(np.arange(0,256),np.arange(0,256))
    #pdb.set_trace()
    obj = xr.Dataset(coords={'XLONG': (['south_north', 'west_east'], lon),
                     'XLAT': (['south_north', 'west_east'], lat),
                     'Time': pd.date_range('2007-06-01', periods=end, freq='H')})
    obj["Objects"] = (['Time','south_north','west_east'],
                       np.array(ds1.variables['fcst_object_id'][sta:end,:,:]))
    ### Add easy variables
    data_dry = []
    for count,var in enumerate(varis):
        print('##### '+var+' #####')
        value = ds.variables[var][sta:end,:,:,:].where(obj.Objects[tiempo[j],:,:]>0)
        #value = value.mean(dim=['south_north','west_east'])
        data_dry.append(value)
    ### Add more complicated variables!
    theta = ds.T[sta:end,:,:,:] + 300
    data_dry.append(theta.where(obj.Objects[tiempo[j],:,:]>0))
    #del(theta)
 
    ###Creating files
    print('Adding to file: Dry')
    ### Dry file
    lon = np.arange(0,nx,1)
    lat = np.arange(0,ny,1)
    ncout = Dataset(mdir+run+'_run/Dry_3D_Lagrangian_'+run+'.nc','w','NETCDF4')
    ncout.createDimension('lon',nx)
    ncout.createDimension('lat',ny)
    ncout.createDimension('bottom_top',elv)
    ncout.createDimension('Time',times);
    lonvar = ncout.createVariable('lon', 'float32', ('lon'));lonvar[:] = lon;
    latvar = ncout.createVariable('lat', 'float32', ('lat'));latvar[:] = lat;    
    elvvar = ncout.createVariable('bottom_top', 'float32', ('bottom_top'));elvvar[:] = elev;
    timevar = ncout.createVariable('Time','float32',('Time'));timevar.setncattr('units',unout);timevar[:]=date2num(dates);
    for i,v in enumerate(varis1):
        print('Including: '+v)
        myvar = ncout.createVariable(v,'float32',('Time','bottom_top','lat','lon'));myvar[:] = data_dry[i][:,0:62,:,:];
    ncout.close();
    

