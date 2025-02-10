from netCDF4 import Dataset,num2date,date2num
import numpy as np
import xarray as xr
from Casicus import CPs_object, labels_cp, cp_geometry
import datetime
import pdb

#--------------------------------------------
#Run
tipo = ''
if tipo == '':
    runs = ['Tho'+tipo+'YSU','GCE'+tipo+'YSU','Mor'+tipo+'YSU','Mor'+tipo+'BL','WSM'+tipo+'YSU','WSM'+tipo+'BL','GCE'+tipo+'BL','Tho'+tipo+'BL', 'ThoYDE', 'ThoBHE', 'ThoYTB', 'ThoBHB']
else:
    runs = ['Tho'+tipo+'YSU','GCE'+tipo+'YSU','Mor'+tipo+'YSU','Mor'+tipo+'BL','WSM'+tipo+'YSU','WSM'+tipo+'BL','GCE'+tipo+'BL','Tho'+tipo+'BL']

#Path
mdir = '/home/tompkins-archive/acasallas/'
user = '/home/netapp-clima/users/acasallas/'

# Constants
rho = 1.225
Rd=287.05
Cp=1005.0
Lv=2.5e6

date_ini = datetime.datetime(2007,6,1)
nx = 256
ny = 256
sta = 0*24
end = 20*24
lon = np.arange(0,nx,1)
lat = np.arange(0,ny,1)

### pcentile for edge detection in gradient 
pcen=60

for i,run in enumerate(runs):
    print('-----> Starting with '+run)
    if run == 'ThoNDYSU' or run == 'ThoNDBL':
        ds = Dataset(user+run+'_run/Variables_'+run+'_all.nc')
    else:
        ds = Dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
    unout = ds.variables['XTIME'].description
    dates = np.array([])
    #Create dates 
    for t1 in range(end-sta):
        dates = np.append(dates, date_ini+datetime.timedelta(minutes=int(ds.variables['XTIME'][t1])))
    del(ds)
    #Read the data!
    if run == 'ThoNDYSU' or run == 'ThoNDBL':
        ds=xr.open_dataset(user+run+'_run/Variables_'+run+'_all.nc')
    else:
        ds=xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
    t2=ds.T2[sta:end,:,:]
    t2=(t2.differentiate('south_north')*t2.differentiate('south_north')+t2.differentiate('west_east')*t2.differentiate('west_east'))**0.5
    #Select coldpools
    print('Selecting coldpools')
    coldpools, distance = CPs_object(t2,pcen)
    mat = np.ones((len(np.arange(sta,end)),nx,ny))
    print('Calculating the labels')
    for t in range(end-sta):
        print('Step '+str(t)+' from '+str(end))
        #Coldpools labels
        mat[t,:,:] = labels_cp(distance[t,:,:],coldpools[t,:,:])
    ###Creating file
    print('Adding Cp objects to file')
    ### Dry file
    if run == 'ThoNDYSU' or run == 'ThoNDBL':
        ncout = Dataset(user+run+'_run/CP_objects_'+run+'_per'+str(pcen)+'.nc','w','NETCDF4')
    else:
        ncout = Dataset(mdir+run+'_run/CP_objects_'+run+'_per'+str(pcen)+'.nc','w','NETCDF4')
    ncout.createDimension('west_east',nx)
    ncout.createDimension('south_north',ny)
    ncout.createDimension('Time',end);
    timevar = ncout.createVariable('Time','float32',('Time'));timevar.setncattr('units',unout);timevar[:]=date2num(dates,unout);
    lonvar = ncout.createVariable('west_east', 'float32', ('west_east'));lonvar[:] = lon;
    latvar = ncout.createVariable('south_north', 'float32', ('south_north'));latvar[:] = lat;
    myvar = ncout.createVariable('CP_object','float32',('Time','south_north','west_east'));myvar[:] = mat[:,:,:];
    ncout.close();
