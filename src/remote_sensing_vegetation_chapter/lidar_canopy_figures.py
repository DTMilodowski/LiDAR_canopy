###############################################################################################################
# This driver function analyses both LiDAR data and field inventory data to produce independent estimates of
# canopy structure.  These are compared against each other and their integrated LAD is compared against LAI
# estimates from hemispherical photographs.
###############################################################################################################
import numpy as np
from scipy import stats
import xarray as xr
#import auxilliary_functions as aux
#import canopy_structure_plots as csp
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

#------------------------------------------------------------------------------------
# DIRECTORIES
# start by defining input files
las_file = '/home/dmilodow/DataStore_DTM/OTHER/RemoteSensingofVegetation_chapter/LiDAR/LNF_FR.las'
dem_file = '/home/dmilodow/DataStore_DTM/OTHER/RemoteSensingofVegetation_chapter/LiDAR/LNF_FR_dem.tif'
dsm_file = '/home/dmilodow/DataStore_DTM/OTHER/RemoteSensingofVegetation_chapter/LiDAR/LNF_FR_dsm.tif'
output_dir = '/home/dmilodow/DataStore_DTM/OTHER/RemoteSensingofVegetation_chapter/figures/'

#-----------------------------------------------------------------------------------
# Plot rasters (Figure 2 from LiDAR work)
dsm = xr.open_rasterio(dsm_file)
dsm.values[dsm.values>905]=np.nan
dsm.values[dsm.values < 0] = np.nan
x = dsm.coords['x'].values.copy()
y = dsm.coords['y'].values.copy()
origin = [x[0],y[0]]
dsm.coords['x']=x-x[0]
dsm.coords['y']=y-y[-1]

dem = xr.open_rasterio(dem_file)
dem.values[dem.values < 0] = np.nan
dem.coords['x']=x-x[0]
dem.coords['y']=y-y[-1]


chm = dsm-dem
mask = np.any((chm.values<0,chm.values>100),axis=0)
chm.values[mask] = np.nan

elev_range =[np.nanmin(dem.values),np.nanmax(dsm.values)]
height_range =[np.nanmin(chm.values),60]
fig,axes= plt.subplots(1,3,figsize = [9,5])
dsm.plot(vmin=elev_range[0],vmax=elev_range[1],cmap = 'viridis',ax=axes[0],
        add_colorbar=True,cbar_kwargs={'label': 'Elevation A.S.L. / m',
                    'orientation':'horizontal'})
dem.plot(vmin=elev_range[0],vmax=elev_range[1],cmap = 'viridis',ax=axes[1],
        add_colorbar=True,cbar_kwargs={'label': 'Elevation A.S.L. / m',
                    'orientation':'horizontal'})
chm.plot(vmin=height_range[0],vmax=height_range[1],cmap = 'viridis',ax=axes[2],
        add_colorbar=True,cbar_kwargs={'label': 'Height above ground / m',
                    'orientation':'horizontal'})

axes[0].set_title('Digital Surface Model',fontsize=12)
axes[1].set_title('Digital Elevation Model',fontsize=12)
axes[2].set_title('Canopy Height Model',fontsize=12)

yticklabels = axes[1].get_yticklabels() + axes[2].get_yticklabels()
axes[0].set_ylabel('Distance / m')
axes[1].set_ylabel('')
axes[2].set_ylabel('')
axes[1].set_xlabel('Distance / m')
axes[0].set_xlabel('')
axes[2].set_xlabel('')
plt.setp(yticklabels,visible=False)
for ax in axes:
    ax.yaxis.set_ticks_position('both')
    ax.set_aspect("equal")
fig.savefig('%s/canopy_height_model.png' % output_dir)
fig.show()

#-----------------------------------------------------------------------------------
# Plot point cloud samples (classification, height above ground)
N=4397200.
S=4396800.
W=647400.
E=647410.
plot_bbox = np.asarray([[W,N],[E,N],[E,S],[W,S]])
pts, starting_ids, trees = io.load_lidar_file_by_polygon(las_file,plot_bbox,max_pts_per_tree = 5*10**5)
pts[:,0]=pts[:,0]-x[0]
pts[:,1]=pts[:,1]-y[-1]

XX,YY = np.meshgrid(x-x[0],y-y[-1])
XX=XX.ravel()
YY=YY.ravel()
ZZ=dem.values.ravel()
mask = np.all((XX<=E-x[0],XX>=W-x[0],YY>=S-y[-1],YY<=N-y[-1]),axis=0)
temp_ids, dem_trees = io.create_KDTree(np.array([XX[mask],YY[mask]]).T)
heights = np.zeros(pts.shape[0])
for pp in range(0,heights.size):
    dist,id = dem_trees[0].query(pts[pp,:2],k=1)
    heights[pp] = pts[pp,2]-ZZ[mask][id]

fig,axes= plt.subplots(2,1,figsize = [8,7])
axes[0].scatter(pts[:,1],pts[:,2],marker='.',s=0.1,c=pts[:,4],cmap = 'viridis_r',
                vmin=0.6)
axes[1].scatter(pts[:,1],pts[:,2],marker='.',s=0.1,c=heights)
axes[0].annotate('a - Point Classifications', fontsize=12, xy=(0.95,0.95), xycoords='axes fraction',horizontalalignment='right', verticalalignment='top')
axes[1].annotate('b - Height Above Ground', fontsize=12, xy=(0.95,0.95), xycoords='axes fraction',horizontalalignment='right', verticalalignment='top')
for ax in axes:
    ax.set_aspect('equal')
fig.savefig('%s/point_cloud.png' % output_dir)
fig.show()

#------------------------------------------------------------------------------------
# PARAMETERS
# define important parameters for canopy profile estimation
max_height = 60
layer_thickness = 1.
heights = np.arange(0,max_height,layer_thickness)+layer_thickness
n_layers = heights.size

N=4396930.
S=4396920.
W=647400.
E=647410.
plot_bbox = np.asarray([[W,N],[E,N],[E,S],[W,S]])
pts, starting_ids, trees = io.load_lidar_file_by_polygon(las_file,plot_bbox,max_pts_per_tree = 5*10**5)
pts[:,0]=pts[:,0]-x[0]
pts[:,1]=pts[:,1]-y[-1]
pts=pts[pts[:,3]==1]


XX,YY = np.meshgrid(x-x[0],y-y[-1])
XX=XX.ravel()
YY=YY.ravel()
ZZ=dem.values.ravel()
mask = np.all((XX<=E-x[0],XX>=W-x[0],YY>=S-y[-1],YY<=N-y[-1]),axis=0)
temp_ids, dem_trees = io.create_KDTree(np.array([XX[mask],YY[mask]]).T)
heights = np.zeros(pts.shape[0])
for pp in range(0,heights.size):
    dist,id = dem_trees[0].query(pts[pp,:2],k=1)
    heights[pp] = pts[pp,2]-ZZ[mask][id]
