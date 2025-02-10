import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops, perimeter
import pandas as pd
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1 import make_axes_locatable
from bisect import bisect_left
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology

### Hov for WRF output
def hov_wrf(sst_ds,tcwv_ds):
    time1 = 0*24
    time2 = len(tcwv_ds.QVAPOR) - 1
    ntime = time2-time1+1
    dtime = 3600 #This depends on the ouput frequency
    difftime = str(dtime)
    fstr = difftime+'-'+str(time1)+'-'+str(time2)
    #Domain and splits
    nx = 552; ny = 399; nz = 32; npts = nx*ny # hardwired for moment
    hovchunk=1656  #1024 # averaging length for chunking, power of 2, larger=smoother.
    nhov=int(npts/hovchunk)
    hovx=100*(np.arange(nhov)+0.5)/nhov
    ilev1=0
    ilev2=8
    #constants
    Rd=287.05
    Cp=1005.0
    Lv=2.5e6
    Cpdsinday=Cp/86400.
    #Start matrixes 
    times=np.array(range(ntime)) # start from day 1 for spin up
    hovtimes=(times+time1)*dtime/86400.
    tcwv_hov=np.zeros([1,ntime,nhov])
    sst_hov=np.zeros([1,ntime,nhov])
    for itime in times:
        print('Time: '+str(itime))
        tcwv = np.array(tcwv_ds.variables['QVAPOR'][itime,:,:])
        sst =  np.array(sst_ds.variables['TSK'][itime,:,:])
        # sort the tcwv
        isort=np.argsort(tcwv.flatten()) # sort tcwv
        sstpert=sst.flatten()-np.mean(sst)
        sst_hov[0,itime,:]=np.mean(sstpert[isort].reshape(-1,hovchunk),axis=1)
        tcwv_hov[0,itime,:]=np.mean(tcwv.flatten()[isort].reshape(-1,hovchunk),axis=1)
    return(sst_hov, tcwv_hov, hovx, ntime)

### Autocorrelation!
def autocorr(x,lags,mode):
    '''numpy.correlate, non partial'''
    mean=x.mean()
    var=np.var(x)
    xp=x-mean
    corr=np.correlate(xp,xp,mode)[len(x)-1:]/var/len(x)
    return corr[:len(lags)]

### The next 3 are to detect the local minima of the coldpools
def detect_local_minima(arr,minmax):
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    if minmax == 0:
        local_min = (filters.minimum_filter(arr, footprint=neighborhood)==arr)
    elif minmax == 1:
        local_min = (filters.maximum_filter(arr, footprint=neighborhood)==arr)
    background = (arr==0)
    eroded_background = morphology.binary_erosion(background, structure=neighborhood, border_value=1)
    detected_minima = local_min ^ eroded_background
    #return np.where(detected_minima)       
    return detected_minima

def number_pixels(tv_fil):
    #### Lets calculate the area of pixels that follow a condition
    cp_n_pix = (tv_fil.notnull()).astype('int')
    cp_n_pix = cp_n_pix.where(cp_n_pix>0)
    return cp_n_pix

def n_cp_f(tv,tv_fil,minmax):
    #### Lets count the number of coldpools, or w or any core
    l_coord = detect_local_minima(tv_fil,minmax)
    lmin = tv.where(tv*l_coord>0)
    n_cp = (lmin.notnull()).astype('int')
    n_cp = (n_cp.where(n_cp>0))
    return n_cp

### Plot the area by the local minima
def make_plot_area(rad,std,mean,run,j):
    plt.plot(rad, label = run, color = 'b')
    plt.plot(rad + std, color = 'b', linestyle = ':', linewidth = 0.7)
    plt.plot(rad - std, color = 'b', linestyle = ':', linewidth = 0.7)
    plt.grid(linestyle = ':', linewidth = 0.5)
    plt.xticks(np.arange(0,len(rad)-1,len(rad)/6.))
    plt.xlim(0,len(rad))
    plt.legend(loc = 'lower left')
    if j == 4:
        plt.xlabel('Hours', fontweight = 'bold')
    else:
        ax.xaxis.set_major_formatter(NullFormatter())

### Vertical int but for xarray
def verint_xarray(ds,var,time,elv):
    #### var have to be the variable for each time
    Rd=287.05
    Cp=1005.0

    press = (ds.P[time,0:elv,:,:] + ds.PB[time,0:elv,:,:])
    T = (ds.T[time,0:elv,:,:]+300)*np.power(press/1.e5,Rd/Cp)
    hgt = (ds.PH[time,0:elv,:,:] + ds.PHB[time,0:elv,:,:]).mean(dim=['south_north','west_east'])/9.81
    diffhgt=np.diff(hgt)    
    vari = (press[:-1,:,:]/(Rd*T[:-1,:,:]))*var[:-1,:,:]*diffhgt[:,None,None]
    del(press); del(T); del(diffhgt)
    varint = vari.sum(dim='bottom_top')
    return(varint)
    
### Usefull to calculate the vertical integration of any wrf variable!
# Input, dataset, name of the variable, start day and for var or string
def vertint(ds,varname,sta_day,tipo):
    Rd=287.05
    Cp=1005.0

    pr= np.array(ds.variables["P"][sta_day:sta_day+1])
    pb= np.array(ds.variables["PB"][sta_day:sta_day+1])
    press=pr+pb
    del(pr)
    del(pb)

    ph= np.array(ds.variables["PH"][sta_day:sta_day+1])
    phb=np.array(ds.variables["PHB"][sta_day:sta_day+1])
    hgt=np.mean(ph+phb,axis=(2,3))/9.81
    del(ph)
    del(phb)

    theta=np.array(ds.variables["T"][sta_day:sta_day+1])

    diffhgt=np.diff(hgt,axis=1)
    del(hgt)
    theta=theta+300. # pot temp  
    #press=np.array(ds.variables["P"][sta_day:sta_day+1])

    # temperature
    temp=theta*np.power(press/1.e5,Rd/Cp)
    del(theta)
    if tipo == 0:
        var=np.array(ds.variables[varname][sta_day:sta_day+1])
        varint=np.sum(press/(Rd*temp)*var*diffhgt[:,:,None,None],axis=1)
    else: 
        var = varname 
        varint=np.sum(press[0,:-1,:,:]/(Rd*temp[0,:-1,:,:])*var[:-1,:,:]*diffhgt[:,:-1,None,None],axis=1)
    #varint=np.sum(press[:,:-1,:,:]/(Rd*temp[:,:-1,:,:])*var[:,:-1,:,:]*diffhgt[:,:,None,None],axis=1)
    #varint[varint<1e-6] = 0
    del(press)
    del(temp)
    del(var)
    del(diffhgt)
    return(varint)

### Usefull to calculate the moving average
# Input array and number of data points for running
def moving_average(a,n):
    N=len(a)
    return np.array([np.mean(a[i:i+n]) for i in np.arange(0,N-n+1)])

##Glitch restart point
def glitch(x): #replace the entire step for glitches in the data
    newx=x.copy()
    for i in range(1,len(x)-1):
        if np.any((x[i]-x[i-1])*(x[i]-x[i+1])>1.0):
            newx[i]=0.5*(x[i-1]+x[i+1])
    return (newx)

def glitch2(x): #replace the entire step for glitches in the data
    newx=x.copy()
    for i in range(1,len(x)-1):
        if np.any((x[i])<0.0):
            newx[i]=x[i+2]
    return (newx)

### Make hovmuller plots!
def make_fig(fig,ax,hovx,times,hov,*args,**kwargs):
    lab=kwargs.get('lab','')
    xlab=kwargs.get('xlab','')
    ylab=kwargs.get('ylab','')
    units=kwargs.get('units','')
    xoff=kwargs.get('xoff',False)
    yoff=kwargs.get('yoff',False)
    kfac=kwargs.get('kfac',1.0)
    xticks=kwargs.get('xticks',[])
    yticks=kwargs.get('yticks',0)
    cbar_pos=kwargs.get('cbar_pos',1) # 0=off, 1=default side, 2=top
    yticki=kwargs.get('yticki',0.25)
    title=kwargs.get('title',"")
    fsize=kwargs.get('fsize',12)
    cmap=kwargs.get('cmap','bwr')
    debug=kwargs.get('debug',False)
    zhov=hov*kfac # may need to flip
    vmax=max(zhov.max(),-zhov.min())
    cont_levs=[0] #np.arange(-110,111,20)
    ymin=np.min(times)
    ymax=np.max(times)

    vmax=kwargs.get('vmax',vmax)
    cont_levs=kwargs.get('conts',cont_levs)
    ymin=kwargs.get('ymin',ymin)
    ymax=kwargs.get('ymax',ymax)
    if yticks==0:
        yticks=np.arange(min(times), max(times)+yticki,yticki)


    if debug:
       print("vmax",vmax)
       exit()

    # image
    img=ax.pcolormesh(hovx,times,zhov,cmap=cmap,
        rasterized=False,vmax=vmax,vmin=-vmax)

    ax.tick_params(axis='both', which='major', labelsize=fsize)
    ax.tick_params(axis='both', which='minor', labelsize=fsize)
    plt.xticks(fontsize=fsize)
    ax.set_title(lab+" "+title,loc="center",fontsize=fsize)
    ax.set_xlabel(xlab,fontsize=fsize)
    ax.set_ylabel(ylab,fontsize=fsize)
    ax.set_yticks(yticks)
    ax.set_ylim(ymin,ymax)
    if len(xticks)>0:
        ax.set_xticks(xticks)

#    plt.tight_layout()

    if yoff==True:
       ax.set_yticklabels([])
    if xoff==True:
       ax.set_xticklabels([])
    if cbar_pos==2:
        cbar_ax=fig.add_axes([0.2, 0.94, 0.6, 0.015])
        cbar=plt.colorbar(img,cax=cbar_ax,orientation='horizontal')
        cbar.set_label(units,rotation=270)
        cbar.ax.tick_params(labelsize=fsize)
    if cbar_pos==1:
        cbar=plt.colorbar(img,ax=ax)
        cbar.set_label(units)
        cbar.ax.tick_params(labelsize=fsize)

    if times.shape[0]>zhov.shape[0]:
       ctimes=times[:-1]+0.5
    else:
       ctimes=times
    # contours, switch off with -999
    #if cont_levs!=-999:
    #    cs=ax.contour(hovx,ctimes,zhov,cont_levs,colors='k',linewidth=0.5)
    #    ax.clabel(cs,inline=1,fontsize=fsize,fmt='%3.0f')
    zhov=None

### Based on Beucler et al 2020, calculate BLW
def calc_MMLi(CRH, c, perc, dx = 2.0e3):
    ### Input:Column Water Vapour (x, y)
    ### Returns: aggregation index calculated as 4*area*pi/(perimeter length)**2
    ###          which compares the length of the boundary between the moist and 
    ###          the dry region with the circumference of the moist region
    ### Parameters: dx grid spacing of input field, perc 
    ### c is True if circle, else is a horizontal band

    binary = np.where(CRH > np.percentile(CRH, perc, interpolation='linear',  axis=(0,1))[np.newaxis, np.newaxis], 1, 0)
    if c == True:
        A_tot = np.sum(binary)*(dx**2) #Area equals to the number of grid points times the gridspacing
        BLW = 2*np.sqrt(np.pi*A_tot)/(dx*perimeter(binary))
    else:
        # 399 is m domain height and 552 its width
        A_tot = 552*dx
        BLW = 2*A_tot/(dx*perimeter(binary))
    return BLW

### Iorg based on Schulz
def nearest_neighbour(points,output_nn=False):
    """
    Calculates for a given list of points (x,y) in a 2D space
    the nearest point.
    
    Returns the nearest neighbour distances and
    optionally two lists containing the original point
    and its neighbour.
    """
    from scipy.spatial import cKDTree
    import numpy as np

    tree = cKDTree(points)
    NN_distance = np.empty((len(points)))
    if output_nn: NN_points = np.empty((len(points),2))

    for i,point in enumerate(points):
        if any(np.isnan(point)): continue
        distance, ind = tree.query(point, k=2)
        NN_distance[i] = distance[1]
        try:
            closestPoint = points[ind[1]]
        except IndexError:
            continue
        if output_nn: NN_points[i] = closestPoint

    if output_nn:
        return NN_distance, NN_points
    else:
        return NN_distance

def NNCDF_ran(l,r=None, rmax=1,resolution=2000):
    """
    l: lambda = number of points / area
    """
    import numpy as np
    if r is None:
        r = np.linspace(0,rmax,resolution)
    return 1-np.exp(-l*np.pi*r**2)

def label_cluster(cluster_mask, stencil=np.ones((1,1))):
    """
    Find coherent clusters.
    
    Method to find and label coherent clusters.
    The stencil can be used to count also masked points
    which are not direct neighbours to the cluster.
    """
    import numpy as np
    from scipy.ndimage import label, binary_dilation
    labels = label(binary_dilation(cluster_mask, structure=stencil))[0]
    return labels

def NNCDF_distances(data, thresholds,dshape=None,stencil=np.ones((1,1))):
    if dshape is not None:
        data = data.reshape(dshape)
    domain_size = np.sqrt(data.shape[0]*data.shape[1])//1
    cluster_mask = np.where((data > thresholds[0]) &\
                            (data < thresholds[1]),1,0)
    if ~np.any(cluster_mask): return np.nan
    cluster_labels = label_cluster(cluster_mask, stencil)
    points = centroid(cluster_mask, cluster_labels).round(0).astype(int)
    distances = nearest_neighbour(points)
    return distances,domain_size

def centroid(cluster_mask, cluster_labels):
    """
    Find centroid of cluster
    """
    import numpy as np
    from scipy.ndimage.measurements import center_of_mass as center
    labels2query = np.unique(cluster_labels[cluster_labels!=0]) # label=0 is the background
    return np.asarray(center(cluster_mask, cluster_labels, labels2query))

def calc_Iorg(distances, domain_size):
    """
    Calculating the Iorg index
    """
    import numpy as np
    n = len(distances)
    return 1-np.sum(NNCDF_ran(n/domain_size**2,r=np.sort(distances)))/n

def calc_Iorg_complete(data, thresholds,dshape=None,mode="normal",stencil=np.ones((1,1))):
    """
    Calculating the Iorg index
    
    Parameters
    ----------
    data : array-like
        input data array, e.g. brightness temperature,
        outgoing longwave radiation, column moisture
        if data is flattened, the dshape parameter has to be defined
    thresholds : tuple, (minimum, maximum) 
        thresholds to mask the data to extract cloud clusters. 
        Data between minimum and maximum is regarded as clouds.
    dshape : tuple, optional
        in case data is flattened (this might be the case
        by applying the function with ndimage.generic_filter)
    mode : string, {"normal","GCM"}, optional
        define mode of Iorg calculation. In case of GCMs
        the option "GCM" might be chosen to handle each gridcell as
    Returns
    -------
    Iorg : scalar
        Organization index
    
    Note
    ----
    It might be necessary to adapt this function in case parts of the
    domain should be excluded e.g. when deep convection is present,
    although one is interested in shallow convection only. The threshold
    may mask it, but left an unnatrual space with no clouds at all.
    """
    import numpy as np
    if dshape is not None:
        data = data.reshape(dshape)
    domain_size = np.sqrt(data.shape[0]*data.shape[1])//1
    cluster_mask = np.where((data > thresholds[0]) &\
                            (data < thresholds[1]),1,0)
    if ~np.any(cluster_mask): return np.nan
    if mode == "normal":
        cluster_labels = label_cluster(cluster_mask, stencil)
    elif mode == "GCM":
        # Give every non-zero mask entry a new number/label
        cluster_labels = np.zeros_like(cluster_mask,dtype="int")
        cluster_labels[np.nonzero(cluster_mask)] = np.arange(1,len(np.nonzero(cluster_mask)[0])+1)
    points = centroid(cluster_mask, cluster_labels).round(0).astype(int)
    distances = nearest_neighbour(points)
    return calc_Iorg(distances, domain_size)

#### Here the watershed
def CPs_object(var,pcen):
    from scipy import ndimage as ndi
    coldpools=var>np.percentile(var,pcen)
    distance = ndi.distance_transform_edt(coldpools)
    return(coldpools,distance)

def CPs_object_thr(var,pcen): ### This depend on a threshold not a percentile!!
    from scipy import ndimage as ndi
    coldpools=var>pcen ### Is > for the coldpools when they come from a file!
    distance = ndi.distance_transform_edt(coldpools)
    return(coldpools,distance)

def labels_cp(distance, coldpools):
    from scipy import ndimage as ndi
    from skimage.feature import canny,peak_local_max
    from skimage.segmentation import active_contour,watershed
    coords = peak_local_max(np.array(distance), footprint=np.ones((3, 3)), labels=np.array(coldpools))
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers, _ = ndi.label(mask)
    labels = watershed(-distance, markers, mask=coldpools)
    return(labels)

def cp_geometry(labels):
    radius = np.array([])
    area = np.array([])
    for obj in range(1,np.max(labels.flatten())):
        z=len(np.argwhere(labels==obj))
        area = np.append(area,z)
        radius = np.append(radius,np.sqrt(z/np.pi))
    return(area,radius)

def nan_out(data):
    data = np.array(data)
    data = data[~np.isnan(data)]
    return(data)

### Bouyancy components calculation, vertical average
def buo_vert_comp(ds, min_times, max_time, min_lons, max_lons, min_lats, max_lats):
    import xarray as xr
    import numpy as np
    import pdb
    #Constants
    Rd=287.05
    Cp=1005.0
    g = 9.81
    bu_rst_mat = np.ones((len(min_times),len(min_lons),62))
    t_rst_mat = np.ones((len(min_times),len(min_lons),62))
    th_rst_mat = np.ones((len(min_times),len(min_lons),62))
    qv_rst_mat = np.ones((len(min_times),len(min_lons),62))
    for j,min_time in enumerate(min_times):
        for i,min_lon in enumerate(min_lons):
            ##### Domain Mean
            qv_mean = ds.QVAPOR[min_time:max_time[j],:,:,:]
            pr_mean = ds.P[min_time:max_time[j],:,:,:] + ds.PB[min_time:max_time[j],:,:,:]
            th_mean = ds.T[min_time:max_time[j],:,:,:] + 300
            tem_mean = th_mean*np.power(pr_mean/1.e5,Rd/Cp)
            tv_mean = (1.0+0.61*qv_mean)*tem_mean
            thv_mean = tv_mean*np.power(1.e5/pr_mean,Rd/Cp)
            del(pr_mean);del(th_mean);del(tv_mean)
            qv_mean = qv_mean.mean(dim=['Time','south_north','west_east'])
            tem_mean = tem_mean.mean(dim=['Time','south_north','west_east'])
            thv_mean = thv_mean.mean(dim=['Time','south_north','west_east'])
            ##### CP zone
            qv = ds.QVAPOR[min_time:max_time[j],:,min_lon:max_lons[i],min_lats[i]:max_lats[i]]
            press = ds.P[min_time:max_time[j],:,min_lon:max_lons[i],min_lats[i]:max_lats[i]] + ds.PB[min_time:max_time[j],:,min_lon:max_lons[i],min_lats[i]:max_lats[i]]
            theta = ds.T[min_time:max_time[j],:,min_lon:max_lons[i],min_lats[i]:max_lats[i]]+300
            temp = theta*np.power(press/1.e5,Rd/Cp)
            tempv = (1.0+0.61*qv)*temp
            thv = tempv*np.power(1.e5/press,Rd/Cp)
            del(tempv)
            thv_pert = thv - thv_mean
            del(thv)
            buo = g*(thv_pert/thv_mean)
            bu_rst_mat[j,i,:] = buo.mean(dim=['Time','south_north','west_east'])
            del(buo)
            thvanom = thv_pert/thv_mean
            th_rst_mat[j,i,:] = thvanom.mean(dim=['Time','south_north','west_east'])
            del(thv_pert);del(thv_mean);del(thvanom)
            tanom = (temp-tem_mean)/(tem_mean)
            t_rst_mat[j,i,:] = tanom.mean(dim=['Time','south_north','west_east'])
            del(tanom)
            qvanom = (0.61*(qv-qv_mean))/(1+0.61*qv_mean)
            qv_rst_mat[j,i,:] = qvanom.mean(dim=['Time','south_north','west_east'])
            del(qv);del(qvanom)
    return(bu_rst_mat, th_rst_mat, t_rst_mat, qv_rst_mat, press)

### Buo components for a period of time
def buo_components(temp_t,press,Tv,qv,Rd,Cp):
    import xarray as xr
    import numpy as np
    thv = Tv*np.power(1.e5/(press*100),Rd/Cp)
    thvanom = (thv - thv.mean())/thv.mean()
    tanom = (temp_t - temp_t.mean())/temp_t.mean()
    qvanom = (0.61*(qv-qv.mean()))/(1+0.61*qv.mean())
    del(thv)
    return(tanom,thvanom,qvanom)

#################### Nature Paper
def cons_count_rever(df,var,count):
    df_copy = pd.concat([df,(df[var].notnull().astype(int)
            .groupby(df[var].isnull().astype(int).cumsum())
            .cumsum().to_frame(count))], axis=1)
    
    df_m24 = df_copy[df_copy['consec_count']>=24] 
    index_24 = df_copy[df_copy['consec_count']==24].index
    
    df_p24 = pd.DataFrame()
    for i in range(len(index_24)):
        z = df.loc[(index_24 - pd.Timedelta(23,unit='h')).strftime('%Y-%m-%d %X')[i] : 
                    index_24.strftime('%Y-%m-%d %X')[i]]
        df_p24 = df_p24.append(z, ignore_index=False)
    dates = pd.concat([df_m24, df_p24], axis=1).index
    
    df_fin = df[df.index.isin(dates)]

    return(df_fin)

def var_data(filename,path,var,varname):
    ds = xr.open_dataset(path+filename)
    if var == 'tpwGrid':
        value = np.array(ds[var]).flatten()
        fecha = [d.strftime('%Y-%m-%d %H:%M:%S') for d in pd.date_range('2016-10-01','2020-01-01', freq='1D')]
    elif var == 'tope':
        value = np.array(ds[var]).flatten()
        fecha = [d.strftime('%Y-%m-%d %H:%M:%S') for d in pd.date_range('2016-01-01','2019-12-31 23:00:00', 
                                                                        freq='1H')]
    else:
        value = np.array(ds[var][0:1461]).flatten()
        fecha = [d.strftime('%Y-%m-%d %H:%M:%S') for d in pd.date_range('2016-01-01','2019-12-31', freq='1D')]
    temp = pd.DataFrame(np.column_stack((fecha,value)),columns = ['Dates',varname])
    temp.to_csv(path+var+'_data.csv', index = False)
    del(temp)
    varia = pd.read_csv(path+var+'_data.csv', index_col=0, parse_dates = True)
    return(varia)

def vert_means(data, spe):
    lista = np.array([])
    for i,var in enumerate(spe):
        lista = np.append(lista, data[var].mean())
    return(lista)

def daily_states(latitud,longitud,fac,daily):
    ### This is used to input the lat, lon and fac, to calculate the reversals and organised days!
    ### This also select the days in the meteorological variables! and calculate the daymean!
    ### C fraction is in other file, as would probably happen also for the cloud properties!
    ### This means that those variables are in another variables, not in the common met one!
    ### The True statement is to select if we want the daily values, or hourly, but careful with other variables
    ### Such as CTOP or cloud properties!
    path = '/home/tompkins-archive/acasallas/RRTMG_data/Area_'+latitud+'_'+longitud+'_ERA5/'
    #### Read Data
    slopes = pd.read_csv(path+'Slopes_data_ERA5.csv', parse_dates = True, index_col = 0)
    met_data= pd.read_csv(path+'Complete_met_data.csv', parse_dates = True, index_col = 0)
    met_data['CIN'] = met_data['CIN'].fillna(value=met_data['CIN'].mean())
    utc = 0
    slopes.index = slopes.index - pd.Timedelta(utc,unit='h')
    met_data.index = met_data.index - pd.Timedelta(utc,unit='h')
    ## Extra vars, they have to be added at the end also
    ctop = var_data('cloud_area.nc', path, 'tope', 'CTop') 
    ctop['CTop'] = ctop['CTop'].fillna(value=ctop['CTop'].mean())
    ### Select the reversals and organize states
    slopes_crit_lo=-fac*slopes.std()[0] # slopes_sd # critical threshold for an event
    slopes_crit_hi=fac*slopes.std()[0] # slopes_sd # critical threshold for an event
    slopes_mask_lo=slopes[slopes<slopes_crit_lo]
    slopes_mask_hi=slopes[slopes>slopes_crit_hi]
    ### Here we produce the counts and select the days with reversals and with organized states!!!
    # First the states in the slopes
    df_slop_lo = cons_count_rever(slopes_mask_lo, 'Slope', 'consec_count')
    df_slop_hi = cons_count_rever(slopes_mask_hi, 'Slope', 'consec_count')
    df_slop_ne = pd.concat([slopes, df_slop_lo]).drop_duplicates(keep=False)
    df_slop_ne = pd.concat([df_slop_ne, df_slop_hi]).drop_duplicates(keep=False)
    #Second select those days in the meteorological data!
    met_lo_data = met_data[met_data.index.isin(df_slop_lo.index)]
    met_hi_data = met_data[met_data.index.isin(df_slop_hi.index)]
    met_ne_data = pd.concat([met_data, met_lo_data]).drop_duplicates(keep=False)
    met_ne_data = pd.concat([met_ne_data, met_hi_data]).drop_duplicates(keep=False)
    if daily == True:
        ### Daymean and selection of variables that are not on the met_data variable
        ### If some variable is missing, it is better to put it here, at the end!
        ### If a variable is added here, it need to be also added to the return! 
        ### CTOP is only daily resolution! so it can not be included in the total met data as hourly!
        daymean_t = met_data.resample('D').mean().dropna()
        daymean_lo = met_lo_data.resample('D').mean().dropna()
        daymean_hi = met_hi_data.resample('D').mean().dropna()
        daymean_ne = met_ne_data.resample('D').mean().dropna()

        daymean_st = slopes.resample('D').mean().dropna()
        daymean_slo = df_slop_lo.resample('D').mean().dropna()
        daymean_shi = df_slop_hi.resample('D').mean().dropna()
        daymean_sne = df_slop_ne.resample('D').mean().dropna()

        #Clod TOP
        daymean_ctt = ctop.resample('D').mean().dropna()
        daymean_ctlo = daymean_ctt[daymean_ctt.index.isin(daymean_lo.index)]
        daymean_cthi = daymean_ctt[daymean_ctt.index.isin(daymean_hi.index)]
        daymean_ctne = daymean_ctt[daymean_ctt.index.isin(daymean_ne.index)]
        return(daymean_t,daymean_lo,daymean_hi,daymean_ne,
               daymean_st,daymean_slo,daymean_shi,daymean_sne,
               daymean_ctt,daymean_ctlo,daymean_cthi,daymean_ctne)
    else: 
        return(slopes, df_slop_lo, df_slop_hi, df_slop_ne,
               met_data, met_lo_data, met_hi_data, met_ne_data)

def RRTMG_nat(sst,moist,ctop,air_T,r_ice,ciwp,r_liq,clwp,press,himawari):
    import climlab
    import pdb
    alb = 0.25
    state = climlab.column_state(lev=np.array(press))
    state.Ts[0] = sst.mean()
    state.Tatm[:] = air_T
    lev = state.Tatm.domain.axes['lev'].points
    if himawari == True:
        cldfrac = ctop.mean()
    else:
        cldfrac = ctop
    
    #CRE_LW = {}
    #CRE_SW = {}
    for thickness in ciwp:
        for i in range(lev.size):
            # Whole-column cloud characteristics
            #  The cloud fraction is a Gaussian bump centered at the current level, when a exp is in there! 
            if himawari == False:
                mycloud = {'cldfrac': cldfrac,
                           'ciwp': np.zeros_like(state.Tatm) + ciwp[thickness],
                           'r_ice': np.zeros_like(state.Tatm) + r_ice,
                           'clwp': np.zeros_like(state.Tatm) + clwp[thickness],
                           'r_liq': np.zeros_like(state.Tatm) + r_liq,}
            else:
                mycloud = {'cldfrac': cldfrac*np.exp(-(lev-lev[i])*2/(2*25.)*2),
                           'ciwp': np.zeros_like(state.Tatm) + ciwp[thickness],
                           'r_ice': np.zeros_like(state.Tatm) + r_ice,
                           'clwp': np.zeros_like(state.Tatm) + clwp[thickness],
                           'r_liq': np.zeros_like(state.Tatm) + r_liq,}
            
            rad_diag = climlab.radiation.RRTMG(state=state, 
                        albedo=alb,
                        specific_humidity=moist,
                        verbose=False,
                        **mycloud)
            #pdb.set_trace()
            rad_diag.compute_diagnostics()
        #CRE_LW[thickness] = (rad_diag.TdotLW - rad_diag.TdotLW_clr)
        #CRE_SW[thickness] = (rad_diag.TdotSW - rad_diag.TdotSW_clr)    
    return(rad_diag)

def cambio_list(data,days):
    chan = []
    for day in days:
        day = int(day)
        chan.append(data[day])
    return(chan)

def exp_rrtmg(sst,hum,tem,cfr,ciwp,clwp,pres_long):
    r_ice = 25.
    r_liq = 14.
    rad = RRTMG_nat(sst,hum,cfr,tem,
                    r_ice,{'med': ciwp*1000,},
                    r_liq,{'med': clwp*1000,},
                    np.array(pres_long), False)
    return(rad.LW_flux_up[0]-1,rad.LW_flux_up[-1]-1,rad.LW_flux_up_clr[0]-1,rad.LW_flux_up_clr[-1]-1)

def battery_exp(varis,rads,days,exp_inp,area,tipo,pres_long,sta,end):
    ### To save the output
    exps_bat = {
               varis[0]+'_OLR':[],varis[0]+'_OLRclr':[], varis[0]+'_SFC':[],varis[0]+'_SFCclr':[]
               } 
    for var in varis:
        if var == 'Clouds':
            chci = cambio_list(exp_inp[var][4],days)
            chcl = cambio_list(exp_inp[var][5],days)
            change = cambio_list(exp_inp[var][3],days)
        elif var == 'SST':
            change = cambio_list(exp_inp[var][0],days)
        elif var == 'Hum':
            change = cambio_list(exp_inp[var][1],days)
        elif var == 'Temp':
            change = cambio_list(exp_inp[var][2],days)
        for day in days:
            day = int(day)
            print(var+' day: '+str(day)+' of '+str(int(days[-1])))
            for j,chan in enumerate(change):
                if var == 'SST':
                    olr,sfc,olrclr,sfcclr = exp_rrtmg(chan,
                                                      exp_inp[var][1][day],
                                                      exp_inp[var][2][day],
                                                      exp_inp[var][3][day],
                                                      exp_inp[var][4][day],
                                                      exp_inp[var][5][day],
                                                      pres_long)
                elif var == 'Hum':
                    olr,sfc,olrclr,sfcclr = exp_rrtmg(exp_inp[var][0][day],
                                                      chan,
                                                      exp_inp[var][2][day],
                                                      exp_inp[var][3][day],
                                                      exp_inp[var][4][day],
                                                      exp_inp[var][5][day],
                                                      pres_long)
                elif var == 'Temp':
                    olr,sfc,olrclr,sfcclr = exp_rrtmg(exp_inp[var][0][day],
                                                      exp_inp[var][1][day],
                                                      chan,
                                                      exp_inp[var][3][day],
                                                      exp_inp[var][4][day],
                                                      exp_inp[var][5][day],
                                                      pres_long)
                elif var == 'Clouds':
                    olr,sfc,olrclr,sfcclr = exp_rrtmg(exp_inp[var][0][day],
                                                      exp_inp[var][1][day],
                                                      exp_inp[var][2][day],
                                                      chan,chci[j],
                                                      chcl[j],pres_long)
                ### Adding results to a dictionary
                exps_bat[var+'_OLR'].append(olr)
                exps_bat[var+'_SFC'].append(sfc)
                exps_bat[var+'_OLRclr'].append(olrclr)
                exps_bat[var+'_SFCclr'].append(sfcclr)
                
    path = '/home/tompkins-archive/acasallas/RRTMG_data/'
    df_exp = pd.DataFrame.from_dict(exps_bat)
    df_exp.to_csv(path+'RRTMG_exp_'+tipo+'_'+area+'_'+varis[0]+'_'+sta+'-'+end+'.csv')
    return(df_exp)

def real_rrtmg(days,exp_inp,area,name,pres_long):
    exps_bat = {'OLR':[],'OLRclr':[], 'SFC':[],'SFCclr':[]}
    for day in days:
        day = int(day)
        print(name+' day: '+str(day)+' of '+str(int(days[-1])))
        olr,sfc,olrclr,sfcclr = exp_rrtmg(exp_inp['Real'][0][day],
                                          exp_inp['Real'][1][day],
                                          exp_inp['Real'][2][day],
                                          exp_inp['Real'][3][day],
                                          exp_inp['Real'][4][day]*1000,
                                          exp_inp['Real'][5][day]*1000,
                                          pres_long)
        exps_bat['OLR'].append(olr)
        exps_bat['SFC'].append(sfc)
        exps_bat['OLRclr'].append(olrclr)
        exps_bat['SFCclr'].append(sfcclr)

    path = '/home/tompkins-archive/acasallas/RRTMG_data/' 
    df_exp = pd.DataFrame.from_dict(exps_bat)
    df_exp.to_csv(path+'RRTMG_'+name+'_'+area+'.csv')
    return(df_exp)

