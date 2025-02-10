import pandas as pd
import numpy as np
import xarray as xr
import pdb
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import scipy.spatial as spatial

#Read the run and the elevation
tipo = '2D' #2D, 3D, TKE

if tipo == '2D':
    runs = ['ThoYSU','ThoBHE','GCEYSU','MorYSU','MorBL','WSMYSU','WSMBL','GCEBL','ThoBL','ThoYDE','ThoYHI','ThoBDI','ThoBHB','ThoYTB']
elif tipo == '3D':
    runs = ['Tho3YSU','Tho3BL','GCE3YSU','Mor3YSU','Mor3BL','WSM3YSU','WSM3BL','GCE3BL']
elif tipo == 'TKE':
    runs = ['ThoTYSU','ThoTBL','GCETYSU','MorTYSU','MorTBL','WSMTYSU','WSMTBL','GCETBL']

elevs = np.arange(0,62)
rho = 1.225
mdir = '/home/tompkins-archive/acasallas/'
scra = '/home/netapp-clima/scratch/acasallas/wrf/'
nx = 256
ny = 256
nt = nx*ny

### Treat the time
time = 0*24
time_e = 44*24

#Distances matrix
dist_conv = np.ones((len(runs),time_e-time))
dist_all = np.ones((len(runs),time_e-time))


for j,run in enumerate(runs):
    print('-----> Starting with '+run)
    ds = Dataset(mdir+run+'_run/Variables_'+run+'_all.nc')
    count = 0
    for i in np.arange(time,time_e):
        print('Step: '+str(i)+' of '+str(time_e)+' in '+run)
        wt = ds.variables['W'][i,19,:,:]
        cnv_coords = np.argwhere(wt > 1)
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
            dist_all[j,i]=0
        else:
            dist_all[j,i]=np.nanmax(f_cnvdst)
        #f_cnvdst*=dxkm
        cnvdst = np.array([])
        for k in range(len(cnv_coords)):
            a = np.delete(j9, list(range(k, j9.shape[0], len(cnv_coords))), axis=0)
            tree=spatial.cKDTree(a)
            d2c,minidx=tree.query(cnv_coords[k,:])
            cnvdst = np.hstack((cnvdst, d2c))
        if len(cnvdst) == 0:
            dist_conv[j,i]=0
        else:
            dist_conv[j,i]=np.nanmax(cnvdst)
        #print(cnvdst)
        #cnvdst*=dxkm

dist_all[dist_all == np.inf] = 0
dist_conv[dist_conv == np.inf] = 0

if tipo == '2D':
    free_cnv = pd.DataFrame(np.column_stack((dist_all[0,:],dist_all[1,:],dist_all[2,:],
                                             dist_all[3,:],dist_all[4,:],dist_all[5,:],
                                             dist_all[6,:],dist_all[7,:],dist_all[8,:],
                                             dist_all[9,:],dist_all[10,:],dist_all[11,:],
                                             dist_all[12,:],dist_all[13,:])), columns = runs)

    cnv_tot = pd.DataFrame(np.column_stack((dist_conv[0,:],dist_conv[1,:],dist_conv[2,:],
                                            dist_conv[3,:],dist_conv[4,:],dist_conv[5,:],
                                            dist_conv[6,:],dist_conv[7,:],dist_conv[8,:],
                                            dist_conv[9,:],dist_conv[10,:],dist_conv[11,:],
                                            dist_conv[12,:],dist_conv[13,:])), columns = runs)
else:
    free_cnv = pd.DataFrame(np.column_stack((dist_all[0,:],dist_all[1,:],dist_all[2,:],
                                             dist_all[3,:],dist_all[4,:],dist_all[5,:],
                                             dist_all[6,:],dist_all[7,:])), columns = runs)

    cnv_tot = pd.DataFrame(np.column_stack((dist_conv[0,:],dist_conv[1,:],dist_conv[2,:],
                                            dist_conv[3,:],dist_conv[4,:],dist_conv[5,:],
                                            dist_conv[6,:],dist_conv[7,:])), columns = runs)

free_cnv.to_csv(mdir+'Data_to_plot/Distances_free_'+tipo+'.csv', index = False)
cnv_tot.to_csv(mdir+'Data_to_plot/Distances_conv_'+tipo+'.csv', index = False)

