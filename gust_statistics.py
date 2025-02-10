import numpy as np
import xarray as xr
from Casicus import labels_cp, cp_geometry, CPs_object_thr, nan_out
import sys
import pdb
import datetime
import time
import pandas as pd

def easy_mean(data,labels):
    var = np.array([])
    for obj in range(1,np.max(labels.flatten())):
        var = np.append(var,data.where(labels == obj).mean())
    return(var)

def sel_variables(ds,sta,end):
    qv=ds.QVAPOR[sta:end,0,:,:]
    u = ds.U[sta:end,0,:,:]
    v = ds.V[sta:end,0,:,:]
    w = ds.W[sta:end,1,:,:]
    qi = ds.QICE[sta:end,0,:,:]
    t2=ds.T2[sta:end,:,:]
    press = (ds.P[sta:end,0,:,:] + ds.PB[sta:end,0,:,:])
    theta_t = ds.T[sta:end,0,:,:] + 300
    temp_t = theta_t*np.power(press/1.e5,Rd/Cp)
    Tv =(1.0+0.61*qv)*temp_t # tempv is virtual T
    thv = Tv*np.power(1.e5/press,Rd/Cp)
    thv_mean = thv.mean()
    thvanom = (thv - thv_mean)/thv_mean
    temp_mean = temp_t.mean()
    tanom = (temp_t - temp_mean)/temp_mean
    qv_mean = qv.mean() 
    qvanom = (0.61*(qv-qv_mean))/(1+0.61*qv_mean)
    #Bouyancy and MSE
    Tv_mean = Tv.mean() 
    buo = 9.81*(Tv - Tv_mean)/Tv_mean
    mse = Cp*temp_t + g*z + Lv*qv - Lf*qi
    del(qi);del(theta_t)
    return(qv,u,v,w,t2,press,temp_t,Tv,thv,thvanom,tanom,qvanom,buo,mse)

def cal_lab_charac(qv,buo,mse,Tv,t2,u,v,w,thvanom,tanom,qvanom,index,t,labels):
    qv_cp = qv[t,:,:].where(labels == index).mean()
    buo_cp = buo[t,:,:].where(labels == index).mean()
    mse_cp = mse[t,:,:].where(labels == index).mean()
    tv_cp = Tv[t,:,:].where(labels == index).mean()
    t2_cp = t2[t,:,:].where(labels == index).mean()
    u_cp = u[t,:,0:256].where(labels == index).mean()
    v_cp = v[t,0:256,:].where(labels == index).mean()
    w_cp = w[t,:,:].where(labels == index).mean()
    thvanom_cp = thvanom[t,:,:].where(labels == index).mean()
    tanom_cp = tanom[t,:,:].where(labels == index).mean()
    qvanom_cp = qvanom[t,:,:].where(labels == index).mean()
    return(qv_cp,buo_cp,mse_cp,tv_cp,t2_cp,u_cp,v_cp,w_cp,thvanom_cp,tanom_cp,qvanom_cp)

#--------------------------------------------
#Paths
mdir = '/home/tompkins-archive/acasallas/'
med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/Coldpools/' 
user = '/home/netapp-clima/users/acasallas/' 
scr = '/home/netapp-clima/scratch/acasallas/'

# Constants
rho = 1.225
Rd=287.05
Cp=1005.0
Lv=2.5e6
Ls = 2.83e6
Lf = Ls - Lv
tipo = 'T'

#Run
if sys.argv == 'All':
    runs = ['Tho'+tipo+'YSU','GCE'+tipo+'YSU','Mor'+tipo+'YSU','Mor'+tipo+'BL','WSM'+tipo+'YSU','WSM'+tipo+'BL','GCE'+tipo+'BL','Tho'+tipo+'BL']
else:
    runs = [sys.argv[1]]

date_ini = datetime.datetime(2007,6,1)
nx = 256
ny = 256
sta = int(sys.argv[2])*24
end = int(sys.argv[3])*24 #10
lon = np.arange(0,nx,1)
lat = np.arange(0,ny,1)
g = 9.81

pcen = 60
z = 40

anom = sys.argv[4]

varis = ['Area','Qv','Buoyancy','MSE','Tv','T2','U','V','W','Tanom','Thvanom','Qvanom']
varis1 = ['Qv','Buoyancy','MSE','Tv','T2','U','V','W','Tanom','Thvanom','Qvanom']
mat_ar = []
mat_ra = []

for i,run in enumerate(runs):
    print('-----> Starting with '+run)
    #Read the data!
    if run == 'ThoNDYSU' or run == 'ThoNDBL':
        ds=xr.open_dataset(user+run+'_run/Variables_'+run+'_all.nc')
        ds1 = xr.open_dataset(user+run+'_run/CP_objects_'+run+'_per'+str(pcen)+'.nc')
    else:
        ds=xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
        ds1 = xr.open_dataset(mdir+run+'_run/CP_objects_'+run+'_per'+str(pcen)+'.nc')
    #Variables!!!
    print('Selecting variables')
    qv,u,v,w,t2,press,temp_t,Tv,thv,thvanom,tanom,qvanom,buo,mse = sel_variables(ds,sta,end) 
    del(ds)
    print('Selecting variables in coldpools')
    area = np.array([])
    radius = np.array([])
    CPs = {'Area':np.array([]),'Qv':np.array([]),'Buoyancy':np.array([]),'MSE':np.array([]),'Tv':np.array([]),'T2':np.array([]),'U':np.array([]),'V':np.array([]),'W':np.array([]),'Tanom':np.array([]),'Thvanom':np.array([]),'Qvanom':np.array([])}
    coldpools, distance = CPs_object_thr(ds1['CP_object'][sta:end,:,:],0)
    del(ds1)
    print('Calculate area and radius')
    for t in range(end-sta):
        start_time = time.time()
        print('Time: '+str(t+sta)+' of '+str(end))
        #Coldpools labels
        labels = labels_cp(distance[t,:,:],coldpools[t,:,:])
        za,zr = cp_geometry(labels)
        area = np.append(area,za)
        radius = np.append(radius,zr)
        if anom == 'True':
            mqv = qv[t,:,:].mean() 
            mbuo = buo[t,:,:].mean()
            mmse = mse[t,:,:].mean()
            mtv = Tv[t,:,:].mean() 
            mt2 = t2[t,:,:].mean()
            mu=u[t,:,0:256].mean()
            mv = v[t,0:256,:].mean()
            mw = w[t,:,:].mean()
            mthvanom = thvanom[t,:,:].mean()
            mtanom = tanom[t,:,:].mean()
            mqvanom = qvanom[t,:,:].mean() 
        for index in range(1,labels.max()):
            CPs['Area'] = np.append(CPs['Area'], qv[t,:,:].where(labels == index).count())
            #CP values in vars!
            if anom == 'True':
                qv_cp,buo_cp,mse_cp,tv_cp,t2_cp,u_cp,v_cp,w_cp,thvanom_cp,tanom_cp,qvanom_cp = cal_lab_charac(qv,buo,mse,Tv,t2,u,v,w,thvanom,tanom,qvanom,index,t,labels)
                means = [mqv,mbuo,mmse,mtv,mt2,mu,mv,mw,mthvanom,mtanom,mqvanom]
                colds = [qv_cp,buo_cp,mse_cp,tv_cp,t2_cp,u_cp,v_cp,w_cp,thvanom_cp,tanom_cp,qvanom_cp]
                for c,var in enumerate(varis1):
                    tmp = colds[c]-means[c]
                    CPs[var] = np.append(CPs[var],tmp)
                del(colds);del(means)
            elif anom == 'False':
                qv_cp,buo_cp,mse_cp,tv_cp,t2_cp,u_cp,v_cp,w_cp,thvanom_cp,tanom_cp,qvanom_cp = cal_lab_charac(qv,buo,mse,Tv,t2,u,v,w,thvanom,tanom,qvanom,index,t,labels)
                colds = [qv_cp,buo_cp,mse_cp,tv_cp,t2_cp,u_cp,v_cp,w_cp,thvanom_cp,tanom_cp,qvanom_cp]
                for c,var in enumerate(varis1):
                    CPs[var] = np.append(CPs[var],colds[c])
                del(colds)
        print("--- %s seconds ---" % (time.time() - start_time))
    mat_ra.append(radius)
    mat_ar.append(area) 
    data = pd.DataFrame.from_dict(CPs)
    if anom == 'True':
        data.to_csv(scr+run+'_CPs_statistics.csv', index=False)
    elif anom == 'False':
        data.to_csv(scr+run+'_CPs_statistics_no_anom.csv', index=False)
        areass = pd.DataFrame(np.column_stack((radius,area)), columns=['Radius','Area'])
        areass.to_csv(scr+run+'_CPs_area_radius.csv', index=False) 

print('Complete')
