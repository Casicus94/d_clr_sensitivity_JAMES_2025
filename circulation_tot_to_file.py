import xarray as xr
import numpy as np
import pandas as pd
import Casicus as casi
import pdb
import sys

scra = '/home/netapp-clima/scratch/acasallas/'
mdir = '/home/tompkins-archive/acasallas/'

#Time
sta = int(sys.argv[1])*24
end = int(sys.argv[2])*24
time1 = sta
time2 = end
ntime=time2-time1+1
dtime=3600 #This depends on the output frequency
difftime=str(dtime)
fstr=difftime+'-'+str(time1)+'-'+str(time2)

#File
runlist = ['ThoYSU','ThoBL','GCEYSU','GCEBL','ThoTYSU','ThoTBL','GCETYSU','GCETBL']

#Domain and splits
nx=256
ny=256
nz=62
npts=nx*ny # hardwired for moment
hovchunk=1024  #1024 # averaging length for chunking, power of 2, larger=smoother.
nhov=int(npts/hovchunk)
hovx=100*(np.arange(nhov)+0.5)/nhov
ilev1=0
ilev2=8
pbl = 20

#Constants
Rd=287.05
Cp=1005.0
Lv=2.5e6
Lf=334000
Cpdsinday=Cp/86400.
qv_press=1100.

#Start matrixes 
times = np.arange(sta,end)
nrun = len(runlist)

#Variables
heights=np.array(range(nz))
lw_hov=np.zeros([nrun,nz,ntime,nhov])
lwc_hov=np.zeros([nrun,nz,ntime,nhov])
lwcl_hov=np.zeros([nrun,nz,ntime,nhov])
sw_hov=np.zeros([nrun,nz,ntime,nhov])
swc_hov=np.zeros([nrun,nz,ntime,nhov])
swcl_hov=np.zeros([nrun,nz,ntime,nhov])
net_hov=np.zeros([nrun,nz,ntime,nhov])
sflw_hov=np.zeros([nrun,nz,ntime,nhov])
sflwc_hov=np.zeros([nrun,nz,ntime,nhov])
sflwcl_hov=np.zeros([nrun,nz,ntime,nhov])
sfsw_hov=np.zeros([nrun,nz,ntime,nhov])
sfswc_hov=np.zeros([nrun,nz,ntime,nhov])
sfswcl_hov=np.zeros([nrun,nz,ntime,nhov])
sfnet_hov=np.zeros([nrun,nz,ntime,nhov])
tcwv_hov=np.zeros([nrun,ntime,nhov])
qv_hov=np.zeros([nrun,nz,ntime,nhov])

for irun,run in enumerate(runlist):
    print('--------------> Reading: '+run)
    tcwv = xr.open_dataset(mdir+run+'_run/'+run+'_vert_dataQVAPOR.nc')
    ds = xr.open_dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
    ps = ds.PSFC[100,:,:]
    znu = ds.ZNU[100,:]*ps[0,0]/100.
    for itime,time in enumerate(np.arange(int(sta),int(end))):
        print('Step: '+str(itime)+' from '+str(end-sta))
        tcwvt = tcwv.QVAPOR[time,:,:]
        isort=np.argsort(tcwvt.values.flatten())
        pr = ds.P[time,:,:,:] + ds.PB[time,:,:,:]
        qvt = ds.QVAPOR[time,:,:,:]
        rnet = ds.RTHRATEN[time,:,:]
        rlw = ds.RTHRATLW[time,:,:]
        rlwc =  ds.RTHRATLWC[time,:,:]
        try:
            rlwcl = rlw-rlwc
            tlw = np.array(ds.RTHRATLW[time,0:61,:,:]/theta.diff(dim='bottom_top'))
        except:
            rlw = (ds.RTHRATLW[time-1,:,:] + ds.RTHRATLW[time+1,:,:])/2
            rlwcl = rlw-rlwc
            tlw = np.array((ds.RTHRATLW[time-1,0:61,:,:]+ds.RTHRATLW[time+1,0:61,:,:])/(2*theta.diff(dim='bottom_top')))
        rsw = ds.RTHRATSW[time,:,:]
        rswc = ds.RTHRATSWC[time,:,:]
        rswcl = rsw-rswc
        ### Calculate circulation
        theta = ds.T[time,:,:,:] + 300
        tnet = np.array(ds.RTHRATEN[time,0:61,:,:]/theta.diff(dim='bottom_top'))
        tnet[tnet==np.inf] = np.nan
        tnet[tnet==-np.inf] = np.nan
        tlw = np.array(ds.RTHRATLW[time,0:61,:,:]/theta.diff(dim='bottom_top'))
        tlw[tlw==np.inf] = np.nan
        tlw[tlw==-np.inf] = np.nan
        tlwc = np.array(ds.RTHRATLWC[time,0:61,:,:]/theta.diff(dim='bottom_top'))
        tlwc[tlwc==np.inf] = np.nan
        tlwc[tlwc==-np.inf] = np.nan
        tlwcl = tlw - tlwc
        tsw = np.array(ds.RTHRATSW[time,0:61,:,:]/theta.diff(dim='bottom_top'))
        tsw[tsw==np.inf] = np.nan
        tsw[tsw==-np.inf] = np.nan
        tswc = np.array(ds.RTHRATSWC[time,0:61,:,:]/theta.diff(dim='bottom_top'))
        tswc[tswc==np.inf] = np.nan
        tswc[tswc==-np.inf] = np.nan
        tswcl = tsw - tswc
        #Averages
        temp = theta[:,:,:]*np.power(pr/1.e5,Rd/Cp)
        tempv=(1.0+0.61*qvt)*temp
        rhohgt=pr/(Rd*tempv)
        rhohgt=rhohgt.mean(dim=['south_north','west_east'])
        qvmean = qvt.mean(dim=['south_north','west_east'])
        for iz in range(nz-1):
            qvpert=100.*(qvt[iz,:,:]-qvmean[iz])/qvmean[iz]
            qvpert=np.array(qvpert).flatten()           
            qv_hov[irun,iz,itime,:]=np.mean(qvpert[isort].reshape(-1,hovchunk),axis=1) 
            #Net
            netper = np.array(rnet[iz,:,:]).flatten()
            net_hov[irun,iz,itime,:] = np.mean(netper[isort].reshape(-1,hovchunk),axis=1) 
            #LW
            lwper = np.array(rlw[iz,:,:]).flatten()
            lw_hov[irun,iz,itime,:] = np.mean(lwper[isort].reshape(-1,hovchunk),axis=1)
            lwcper = np.array(rlwc[iz,:,:]).flatten() 
            lwc_hov[irun,iz,itime,:] = np.mean(lwcper[isort].reshape(-1,hovchunk),axis=1)
            lwclper = np.array(rlwcl[iz,:,:]).flatten()
            lwcl_hov[irun,iz,itime,:] = np.mean(lwclper[isort].reshape(-1,hovchunk),axis=1)
            #SW
            swper = np.array(rsw[iz,:,:]).flatten()
            sw_hov[irun,iz,itime,:] = np.mean(swper[isort].reshape(-1,hovchunk),axis=1)
            swcper = np.array(rswc[iz,:,:]).flatten()
            swc_hov[irun,iz,itime,:] = np.mean(swcper[isort].reshape(-1,hovchunk),axis=1)
            swclper = np.array(rswcl[iz,:,:]).flatten()
            swcl_hov[irun,iz,itime,:] = np.mean(swclper[isort].reshape(-1,hovchunk),axis=1)
            ### Circulation
            #Net
            wnetf = tnet[iz,:,:].flatten()
            sfnet_hov[irun,iz,itime,:] = np.nanmean(wnetf[isort].reshape(-1,hovchunk),axis=1)*np.array(rhohgt[iz])
            sfnet_hov[irun,iz,itime,:]=np.cumsum(sfnet_hov[irun,iz,itime,:]) 
            #LW
            wlwf = tlw[iz,:,:].flatten()
            sflw_hov[irun,iz,itime,:] = np.nanmean(wlwf[isort].reshape(-1,hovchunk),axis=1)*np.array(rhohgt[iz])
            sflw_hov[irun,iz,itime,:]=np.cumsum(sflw_hov[irun,iz,itime,:])
            wlwcf = tlwc[iz,:,:].flatten() 
            sflwc_hov[irun,iz,itime,:] = np.nanmean(wlwcf[isort].reshape(-1,hovchunk),axis=1)*np.array(rhohgt[iz])
            sflwc_hov[irun,iz,itime,:]=np.cumsum(sflwc_hov[irun,iz,itime,:])
            wlwclf = tlwcl[iz,:,:].flatten()
            sflwcl_hov[irun,iz,itime,:] = np.nanmean(wlwclf[isort].reshape(-1,hovchunk),axis=1)*np.array(rhohgt[iz])
            sflwcl_hov[irun,iz,itime,:]=np.cumsum(sflwcl_hov[irun,iz,itime,:])
            #SW
            wswf = tsw[iz,:,:].flatten() 
            sfsw_hov[irun,iz,itime,:] = np.nanmean(wswf[isort].reshape(-1,hovchunk),axis=1)*np.array(rhohgt[iz])
            sfsw_hov[irun,iz,itime,:]=np.cumsum(sfsw_hov[irun,iz,itime,:])
            wswcf = tswc[iz,:,:].flatten()
            sfswc_hov[irun,iz,itime,:] = np.nanmean(wswcf[isort].reshape(-1,hovchunk),axis=1)*np.array(rhohgt[iz])
            sfswc_hov[irun,iz,itime,:]=np.cumsum(sfswc_hov[irun,iz,itime,:])
            wswclf = tswcl[iz,:,:].flatten()
            sfswcl_hov[irun,iz,itime,:] = np.nanmean(wswclf[isort].reshape(-1,hovchunk),axis=1)*np.array(rhohgt[iz])
            sfswcl_hov[irun,iz,itime,:]=np.cumsum(sfswcl_hov[irun,iz,itime,:])

lista = [lw_hov, lwc_hov, lwcl_hov, sw_hov, swc_hov, swcl_hov, net_hov, sflw_hov, sflwc_hov, sflwcl_hov, sfsw_hov, sfswc_hov, sfswcl_hov, sfnet_hov]
names = ['lw','lwc','lwcl','sw','swc','swcl','net', 'clw','clwc','clwcl','csw','cswc','cswcl','cnet']

for irun,run in enumerate(runlist):
    for il,lis in enumerate(lista):
        tmp = pd.DataFrame(np.mean(lis[irun,:,:,:],1))
        tmp.to_csv(scra+'circulation_Qrad/'+run+'_'+names[il]+'_Qrad_'+str((int(sta/24)))+'_'+str((int(end/24)))+'.csv')
        tmp.to_csv(mdir+'Data_to_plot/circulation_Qrad/'+run+'_Qrad_'+names[il]+'_'+str((int(sta/24)))+'_'+str((int(end/24)))+'.csv')

