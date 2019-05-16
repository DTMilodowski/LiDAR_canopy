"""
A set of scripts to generate basic plots of the SEOS scottish forests LiDAR
surveys
--------------------------------------------------------------------------------
D.T.Milodowski, 16/05/2019
"""
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt

# Open the datasets
file_id = 'carbomap_site1_5m'
path2data = '/exports/csce/datastore/geos/users/dmilodow/FOREST2020/LiDAR/carbomap_highlands/canopy_metrics/'
path2fig = '/exports/csce/datastore/geos/users/dmilodow/FOREST2020/LiDAR/carbomap_highlands/figures/'
pad_file = '%s%s_pad.tif' % (path2data,file_id)
pai_file = '%s%s_pai.tif' % (path2data,file_id)
pd_file = '%s%s_pulse_density.tif' % (path2data,file_id)

pad = xr.open_rasterio(pad_file)
pai = xr.open_rasterio(pai_file).sel(band=1)
pd = xr.open_rasterio(pd_file).sel(band=1)
pd.values[pd.values==0]=np.nan

# Plot facet_grid for PAD at each height interfal
#fig1,ax = plt.subplots(1,1,figsize=(10,10))
p = pad.plot(x='x',y='y',col='band',col_wrap=8,vmin=0,vmax=1,cmap='plasma',cbar_kwargs={'label':'PAI for each 1m slice in canopy'})
for a in p.axes.flat:
    a.set_aspect('equal', 'box-forced')

plt.savefig('%s%s_pai_by_canopy_layer.png' % (path2fig,file_id))
plt.show()

p = pd.plot(x='x',y='y',vmin=0,vmax=1000,cmap='plasma',cbar_kwargs={'label':'Pulse density / pulses m$^{-2}$'})
p.axes.set_aspect('equal', 'box-forced')
plt.savefig('%s%s_pulse_density.png' % (path2fig,file_id))
plt.show()
