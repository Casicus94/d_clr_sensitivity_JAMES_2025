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
runs = ['ThoYSU']

#Some constants and dimensions
times = 480
date_ini = datetime.datetime(2007,6,1)
elv = 62
elev = np.arange(0,elv)
#Path
mdir = '/home/tompkins-archive/acasallas/'
#Variables to file
varis = ['QVAPOR','QCLOUD','QICE','QSNOW','QRAIN','U','V',
         'T','W','P','PB','PH','PHB','RTHRATEN',
         'RTHRATLW','RTHRATLWC','RTHRATSW','RTHRATSWC']
#Dry and Moist patches coordinates
dries_x = [29,2,135,82,51,185,74,142,250]
dries_y = [243,185,217,161,91,201,92,98,79]
moist_x = [143,206,74,227,36,237,64,173,17]
moist_y = [159,19,227,134,33,238,6,166,143]

for j,run in enumerate(runs):
    print('-----> Starting with '+run)
    ds = Dataset(mdir+run+'_run/Variables_'+run+'_all.nc') 
    unout = ds.variables['XTIME'].description
    dates = np.array([])
    #Create dates 
    for t in range(times):
        dates =np.append(dates, date_ini+datetime.timedelta(minutes=int(ds.variables['XTIME'][t])))
    #Create list for each variable!
    data_dry = []
    data_moi = [] 
    for count,var in enumerate(varis):
        print('##### '+var+' #####')
        for coord,dry_x in enumerate(dries_x):
            print('Coordinate: '+str(coord+1)+' of '+str(len(dries_x)))
            #Select coordinate, 0 since its the beggining!
            if coord == 0:
                ds_dry = ds.variables[var][0:times,:,dries_y[coord],dry_x]
                ds_moi = ds.variables[var][0:times,:,moist_y[coord],moist_x[coord]]
            else:
                ds_dry = ds_dry + ds.variables[var][0:times,:,dries_y[coord],dry_x]
                ds_moi = ds_moi + ds.variables[var][0:times,:,moist_y[coord],moist_x[coord]]
        print('Calculating composite')
        #Composite, so this is the mean!
        ds_dry = ds_dry/len(dries_x)
        ds_moi = ds_moi/len(moist_x) 
        #Append variable data to the list
        data_dry.append(ds_dry)
        data_moi.append(ds_moi)
    del(ds)
    
    ###Creating files
    print('Adding to file: Dry')
    ### Dry file
    ncout = Dataset(mdir+run+'_run/Dry_patches_'+run+'.nc','w','NETCDF4')
    ncout.createDimension('bottom_top',elv)
    ncout.createDimension('Time',times);
    elvvar = ncout.createVariable('bottom_top', 'float32', ('bottom_top'));elvvar[:] = elev;
    timevar = ncout.createVariable('Time','float32',('Time'));timevar.setncattr('units',unout);timevar[:]=date2num(dates);
    for i,v in enumerate(varis):
        print('Including: '+v)
        myvar = ncout.createVariable(v,'float32',('Time','bottom_top'));myvar[:] = data_dry[i][:,0:62];
    ncout.close();
    
    print('Adding to file: Moist')
    ### Moist file
    ncout = Dataset(mdir+run+'_run/Moist_patches_'+run+'.nc','w','NETCDF4')
    ncout.createDimension('bottom_top',elv)
    ncout.createDimension('Time',times);
    elvvar = ncout.createVariable('bottom_top', 'float32', ('bottom_top'));elvvar[:] = elev;
    timevar = ncout.createVariable('Time','float32',('Time'));timevar.setncattr('units',unout);timevar[:]=date2num(dates);
    for i,v in enumerate(varis):
        print('Including: '+v)
        myvar = ncout.createVariable(v,'float32',('Time','bottom_top'));myvar[:] = data_moi[i][:,0:62];
    ncout.close();



