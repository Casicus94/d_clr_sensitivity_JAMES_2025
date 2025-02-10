import pdb
import os
import xarray as xr
import numpy as np
import pandas as pd

### Here (following link) is how we can mask a region depending on a shapefile
## https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html
# Path to the input and output files!

path  = '/home/netapp-clima/scratch/acasallas/Rosario_data/ParamoFireAtlas-20231211T193915Z-001/ParamoFireAtlas/'
scr = '/home/netapp-clima/scratch/acasallas/Rosario_data/ParamoFireAtlas-20231211T193915Z-001/'

# Create a list with the names of all the paramos, depend on the directory!
paramos = os.listdir(path)
   
for i,paramo in enumerate(paramos):
    print('--------------> Paramo: '+paramo)
    try:
        data = pd.read_csv(scr+'Data_Fire_Num_Extension/'+paramo+'fire_num_ext.csv')
    except:
        ### This for is to make a list of all the files with netcdf extension
        files = []
        for filename in os.listdir(path+paramo):
            if filename.endswith('.nc'):
                files.append(filename)
        ### Dictionary for fire number and fire extension
        dict_par = {'Dates':[], 'Fire_number':np.array([]), 'Fire_extension':np.array([])}
        for filename in files:
            ds = xr.open_dataset(path+paramo+'/'+filename)
            dict_par['Dates'].append(int(filename[-11:-3]))
            ### Number of fires
            tmp = ds.Objects.where(ds.Objects>0)
            # The last value of the labels, aka the max, imply the number of labels
            # Which is the same as the number of fires
            dict_par['Fire_number'] = np.append(dict_par['Fire_number'],tmp.max())
            ### Size in Hectares, 30*30m2 is the pixel size, we converted to km2, 1000*1000
            ### Then we convert it to hectares, dividing by 100
            tmp = ds.test.where(ds.Objects>0).sum()*30*30/(1000*1000*100) 
            dict_par['Fire_extension'] = np.append(dict_par['Fire_extension'],tmp)
        ### Here we save the data in a csv file for plotting, this makes the process faster
        ### It also makes everything less expensive in terms of computing power
        print('Data being save into a csv file')
        # Including the data into a DataFrame for simplicity
        data = pd.DataFrame.from_dict(dict_par)
        ### Since the for is not sorted, we need to sort the dataframe (which is also easier)
        ### It is important to have the dates sorted, to then make calculations for tendencies
        sort = data.sort_values('Dates',ascending=True)
        sort.to_csv(scr+'Data_Fire_Num_Extension/'+paramo+'fire_num_ext.csv', index=False)


