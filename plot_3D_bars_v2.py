import xarray as xr
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pdb

med = '/media/acasallas/ALEJO_HD/Mac_backup/Documents/PhD/Plots/'
scr = '/home/netapp-clima/scratch/acasallas/'
mdir = '/home/tompkins-archive/acasallas/'
horz = ['Smag2D','TKE','Smag3D']
pbls = ['BL', 'YSU']
micro = ['WSM', 'GCE', 'Tho', 'Mor']*len(pbls)

htipo = ['','T','3']
runs = ['WSM','GCE','Tho','Mor']
kk = '_XK'
IQR = np.ones((len(micro),len(htipo)))*np.nan
kh = np.ones((len(micro),len(htipo)))*np.nan

xh = xr.open_dataset(mdir+'Data_to_plot/Eddy_viscosity_tot_0-1056.nc')
iq = pd.read_csv(mdir+'Data_to_plot/IQR_data_tot.csv')
dia = 5

for j,hor in enumerate(htipo):
    count = 0
    for pbl in pbls:
        for run in runs:
            IQR[count,j] = iq[run+hor+pbl][-dia*24:].mean()
            xh_mean = xh[run+hor+pbl+kk][:,:].mean(dim='elev') 
            kh[count,j] = xh_mean.mean()
            count+=1

#pdb.set_trace()
### Matrix dimensions
lx= len(horz)
ly= len(micro)
xpos = np.arange(lx)    # Set up a mesh of positions
ypos = np.arange(ly)
xpos, ypos = np.meshgrid(xpos, ypos)
### 1D arrays
xpos = xpos.flatten()   
ypos = ypos.flatten()
zpos = np.zeros(lx*ly)
# Size
dx = 0.5*np.ones_like(zpos)
dy = dx.copy()
dz = IQR.flatten()

### Ticks
ticksy = np.arange(0.45, ly, 1)
ticksx = np.arange(0.35, lx, 1)
### This is to create the cmap and the colorbar
top = cm.get_cmap('gist_ncar', 128)
bottom = cm.get_cmap('Spectral_r', 128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)), bottom(np.linspace(0, 1, 128))))
del(top);del(bottom)
#cmap = ListedColormap(newcolors, name='Double')
cmap = cm.get_cmap('gist_ncar') # Get desired colormap
max_b = np.max(kh.flatten())   # get range of colorbars
min_b = np.min(kh.flatten())
rgba = [cmap((k-min_b)/max_b) for k in kh.flatten()]
norm = Normalize(vmin=min(kh.flatten()), vmax=max(kh.flatten()))
col = cmap(norm(kh.flatten()))
### Y ticks colors
colors = ['k','k', 'k','k', 'purple','purple', 'purple','purple']

### Plot
fig = plt.figure(figsize=(12,6))
ax = Axes3D(fig)
for i in range(ly*lx):
    cs = ax.bar3d(xpos[i],ypos[i],zpos[i], dx[i], dy[i], dz[i], color=rgba[i], alpha = 0.6, shade = True, edgecolor='black')
plt.xticks(ticksx, horz)
plt.yticks(ticksy, micro)
ax.text(-1,0.9,0,'BL', color = 'k', fontweight = 'bold')
ax.text(-1.07,4.5,0,'YSU', color = 'purple', fontweight = 'bold')
for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), colors):
    ticklabel.set_color(tickcolor)
#ax.set_xlabel('Eddy Viscosities schemes')
#ax.set_ylabel('Microphysics schemes')
ax.set_zlabel('TCWV-NMTP (mm)')
# Colorbar
sc = cm.ScalarMappable(cmap=cmap,norm=norm)
sc.set_array([])
cbar = plt.colorbar(sc, shrink=0.7)
cbar.ax.set_ylabel('K$_{h}$')
ax.view_init(19,230)
plt.savefig(scr+'3D_bars_v2.jpg',dpi=300,bbox_inches='tight')
plt.savefig(scr+'3D_bars_v2.pdf',dpi=300,bbox_inches='tight')
plt.show()

