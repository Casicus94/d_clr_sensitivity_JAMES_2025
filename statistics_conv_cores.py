import xarray as xr
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import pdb
import Casicus as casi
import scipy.spatial as spatial
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
from wrf import pw,rh,tvirtual

#Read the run and the elevation
runs = ['WSMYSU', 'WSMBL',
        'GCEYSU', 'GCEBL',
        'ThoYSU', 'ThoBL',
        'MorYSU', 'MorBL']

# Constants
g = 9.81
rho = 1.225
Rd=287.05
Cp=1005.0
Lv=2.5e6
Ls = 2.83e6
Lf = Ls - Lv

z =np.array([40,87.5,141,233,333,436,542,646,752,858,961,
       1060,1213,1421,1623,1825,2031,2241,2450,2657,
       2863,3068,3271,3472,3678,3887,4094,4297,4500,
       4705,4910,5116,5471,5978,6488,6997,7500,8005,
       8513,9013,9509,10009,10503,10988,11468,11949,
       12425,12896,13165,13356,13829,14294,14754,15447,
       16374,17302,18242.6,19190.13,20148.62,21084.01,
       21990.93,22890.83])

mdir = '/home/tompkins-archive/acasallas/'
scra = '/home/netapp-clima/scratch/acasallas/wrf/'

### Treat the time
times = [2*24, 2*24, 5*24, 10*24, 15*24, 20*24, 30*24, 40*24]
time_e = [4*24, 3*24, 10*24, 15*24, 20*24, 25*24, 35*24, 44*24]
level = 8
elv = 18
nx = 256
ny = 256
nt = nx*ny

varis = ['W','QVAPOR','QCLOUD','QICE','QSNOW','QRAIN','MSE','Buoyancy',
         'RTHRATEN','RTHRATLW','RTHRATLWC','RTHRATSW','LWCRE',
         'RTHRATSWC','Cloud Mass Flux','Conv Mass Flux',
         'Downdraft Mass Flux', 'RTHBLTEN']

for i,run in enumerate(runs):
    print('-----> Starting with '+run)
    for it,time in enumerate(times):
        print('----> Time: '+str(int(time/24)))
        ##### Create matrix
        dist_conv = np.ones((len(varis),62)) #Number of variables, number of levels
        dist_doma = np.ones((len(varis)-3,62))
        ########## The idea is to check some variables in the places where W > 1
        # Lets try to calculate the vertical profiles!
        ########## Read data for the other variables!
        ds = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
        ##### Vertical velocity in the updrafts
        wt = ds.W[time:time_e[it],19,:,:]
        print('--> W')
        w = ds.W[time:time_e[it],:,:,:]
        wc = w.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        w = w.mean(dim=['Time','south_north','west_east'])
        dist_doma[0,:] = w[0:-1]
        dist_conv[0,:] = wc[0:-1] 
        del(w);del(wc)
        print('--> QVAPOR')
        qv = ds.QVAPOR[time:time_e[it],:,:,:]
        qvc = qv.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        qv = qv.mean(dim=['Time','south_north','west_east'])
        dist_doma[1,:] = qv
        dist_conv[1,:] = qvc
        del(qv);del(qvc)
        print('--> QCLOUD')
        ql = ds.QCLOUD[time:time_e[it],:,:,:]
        qlc = ql.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        ql = ql.mean(dim=['Time','south_north','west_east'])
        dist_doma[2,:] = ql
        dist_conv[2,:] = qlc
        del(ql);del(qlc)
        print('--> QICE')
        qi = ds.QICE[time:time_e[it],:,:,:]
        qic = qi.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        qi = qi.mean(dim=['Time','south_north','west_east'])
        dist_doma[3,:] = qi
        dist_conv[3,:] = qic
        del(qi);del(qic)
        print('--> QSNOW')
        qs = ds.QSNOW[time:time_e[it],:,:,:]
        qsc = qs.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        qs = qs.mean(dim=['Time','south_north','west_east'])
        dist_doma[4,:] = qs
        dist_conv[4,:] = qsc
        del(qs);del(qsc)
        print('--> QRAIN')
        qr = ds.QRAIN[time:time_e[it],:,:,:]
        qrc = qr.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        qr = qr.mean(dim=['Time','south_north','west_east'])
        dist_doma[5,:] = qr
        dist_conv[5,:] = qrc
        del(qr);del(qrc)
        print('--> MSE')
        press = ds.P[time:time_e[it],:,:,:] + ds.PB[time:time_e[it],:,:,:]
        theta = ds.T[time:time_e[it],:,:,:] + 300
        qv1 = ds.QVAPOR[time:time_e[it],:,:,:]
        temp = theta*np.power(press/1.e5,Rd/Cp)
        qi1 = ds.QICE[time:time_e[it],:,:,:]
        mse = Cp*temp + g*z[None,:,None,None] + Lv*qv1 - Lf*qi1
        del(qi1);del(theta) 
        tempv =(1.0+0.61*qv1)*temp # tempv is virtual T
        del(qv1);del(temp)
        msec = mse.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        mse = mse.mean(dim=['Time','south_north','west_east'])
        dist_doma[6,:] = mse
        dist_conv[6,:] = msec
        del(mse);del(msec)
        print('--> Buoyancy')
        thv = tempv*np.power(1.e5/press*100,Rd/Cp)
        del(tempv)
        press = press.mean(dim=['Time','south_north','west_east'])
        thv_mean = thv.mean(dim=['south_north','west_east'])
        thv_pert = thv - thv_mean
        del(thv)
        buo = g*(thv_pert/thv_mean)
        del(thv_mean);del(thv_pert)
        buoc = buo.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        buo = buo.mean(dim=['Time','south_north','west_east'])
        dist_doma[7,:] = buo
        dist_conv[7,:] = buoc
        del(buo);del(buoc)
        print('--> Rad Cooling')
        print('--> Net')
        net = ds.RTHRATEN[time:time_e[it],:,:,:]
        netc = net.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        net = net.mean(dim=['Time','south_north','west_east'])
        dist_doma[8,:] = net
        dist_conv[8,:] = netc 
        del(net);del(netc)
        print('--> LW')
        LW = ds.RTHRATLW[time:time_e[it],:,:,:]
        LCW = LW.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        LWC = ds.RTHRATLWC[time:time_e[it],:,:,:]
        LCWC = LWC.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        LWCRE = LW - LWC
        LW = LW.mean(dim=['Time','south_north','west_east'])
        dist_doma[9,:] = LW
        dist_conv[9,:] = LCW
        del(LW);del(LCW)
        LWC = LWC.mean(dim=['Time','south_north','west_east'])
        dist_doma[10,:] = LWC
        dist_conv[10,:] = LCWC
        del(LWC);del(LCWC)
        LWCREc = LWCRE.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        LWCRE = LWCRE.mean(dim=['Time','south_north','west_east'])
        dist_doma[11,:] = LWCRE
        dist_conv[11,:] = LWCREc
        del(LWCRE);del(LWCREc)
        print('--> SW')
        SW = ds.RTHRATSW[time:time_e[it],:,:,:]
        SCW = SW.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        SW = SW.mean(dim=['Time','south_north','west_east'])
        dist_doma[12,:] = SW
        dist_conv[12,:] = SCW
        del(SW);del(SCW)
        SWC = ds.RTHRATSWC[time:time_e[it],:,:,:]
        SCWC = SWC.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        SWC = SWC.mean(dim=['Time','south_north','west_east'])
        dist_doma[13,:] = SWC
        dist_conv[13,:] = SCWC
        del(SWC);del(SCWC)
        print('--> Mass Flux')
        ww = ds.W[time:time_e[it],:,:,:]
        del(ds)
        w = (ww.where(wt > 0))
        wwc = (ww.where(wt > 1))
        wd = (ww.where(wt < -0.5))
        n = casi.number_pixels(w).sum()
        nd = casi.number_pixels(wd).sum()
        nc = casi.number_pixels(wwc).sum()
        sig = n/nt
        sigd = nd/nt
        sigc = nc/nt
        m = rho*w*sig
        md = rho*wd*sigd
        mc = rho*wwc*sigc
        dist_conv[14,:] = m.mean(dim=['Time','south_north','west_east'])[0:-1]
        dist_conv[15,:] = mc.mean(dim=['Time','south_north','west_east'])[0:-1]
        dist_conv[16,:] = md.mean(dim=['Time','south_north','west_east'])[0:-1]
        del(m);del(md);del(mc);del(sig);del(sigd);del(sigc)
        del(n);del(nd);del(nc);del(ww);del(w);del(wwc);del(wd)
        print('Theta tendency PBL')
        ds = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
        pbl = ds.RTHBLTEN[time:time_e[it],:,:,:]
        pblc = pbl.where(wt > 1).mean(dim=['Time','south_north','west_east'])
        pbl = pbl.mean(dim=['Time','south_north','west_east'])
        dist_doma[14,:] = pbl
        dist_conv[17,:] = pblc
        del(pbl);del(pblc);del(ds)
        print('--> Creating file')
        print('--> Conv cores')
        df_conv = pd.DataFrame(np.column_stack((dist_conv)),columns=varis)
        df_conv.to_csv(mdir+'Data_to_plot/Conv_cores_stats_'+run+'_day_'+str(int(time/24))+'_to_'+str(int(time_e[it]/24))+'.csv', index = False)
        print('--> Domain')
        df_doma = pd.DataFrame(np.column_stack((dist_doma)),columns=varis[0:len(varis)-3])
        df_doma.to_csv(mdir+'Data_to_plot/Conv_domain_stats_'+run+'_day_'+str(int(time/24))+'_to_'+str(int(time_e[it]/24))+'.csv', index = False)
        del(df_conv);del(df_doma)

print('##### Complete #####')

