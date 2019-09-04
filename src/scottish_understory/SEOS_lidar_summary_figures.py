"""
SEOS_lidar_summary_figures.py
--------------------------------------------------------------------------------
Code to generate a suite of figures summarising the results from the SEOS
UAV LiDAR surveys.
Analysis focussed on site 1
Encompasses:
- Canopy density profiles, inverted from point clouds
    - a random selection of 5m x 5m resolution profiles
- Sensitivity analysis to spatial resolution (1 ha)
- Sensitivity analysis to pulse density (1 ha)
- Sensitivity analysis to pulse density (5m x 5m individual profiles)
- 3D mapping of canopy density
--------------------------------------------------------------------------------
"""
# import standard libraries
import os
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import seaborn as sns
sns.set()

# import general LiDAR libraries
import sys
sys.path.append('../')
import LiDAR_io as io
import LiDAR_tools as lidar

# import SEOS specific libraries
import SEOS_sensitivity_analysis_plots as splt

"""
PART A: LOAD DATA
--------------------------------------------------------------------------------
Load in:
LiDAR data
Sensitivity analyses results
--------------------------------------------------------------------------------
"""
# Define paths etc
file_id = 'carbomap_site1_5m'
path2data = '/exports/csce/datastore/geos/users/dmilodow/FOREST2020/LiDAR/carbomap_highlands/canopy_metrics/'
path2fig = '/exports/csce/datastore/geos/users/dmilodow/FOREST2020/LiDAR/carbomap_highlands/figures/'
os.mkdir('%s/sensitivity_analysis' % path2fig)
las_list = '../las_lists/carbomap_site1_lastiles.txt'
pad_file = '%s%s_pad.tif' % (path2data,file_id)

# Load LiDAR data by bounding box
N=6206000.
S=6205900.
W=344520.
E=344620.
plot_bbox = np.asarray([[W,N],[E,N],[E,S],[W,S]])
pts, starting_ids, trees = io.load_lidar_data_by_polygon(las_list,plot_bbox,laz_files=False,max_pts_per_tree = 5*10**5)
N_trees = len(trees)
# filter LiDAR data by return classification
pts[np.any((pts[:,4]==3,pts[:,4]==4,pts[:,4]==5),axis=0),4]=1
pts[pts[:,2]<0,2]=0

# Load sensitivity analysis to spatial resolution
resolution_sensitivity = np.load('SEOS_MH_sensitivity_resolution.npy')[()]
density_sensitivity = np.load('SEOS_MH_sensitivity_pulse_density.npy')[()]

# Load 3D maps of canopy density
pad = xr.open_rasterio(pad_file)

# other parameters
max_height = 40.
layer_thickness = 0.5
heights = np.arange(0.,max_height,layer_thickness)+0.5
n_layers = heights.size
plot_width = 100.
sample_res = 5.
kappa = 1.
n_iter = 100

"""
PART B: Plot canopy profiles
--------------------------------------------------------------------------------
A couple of plots here
1) point cloud transect; 1 ha vertical point cloud distributions; 1 ha average
    canopy profile plus bootstrap 95% CI
2) Random sample of individual 5m x 5m profiles from within the 1 ha plot
--------------------------------------------------------------------------------
"""


"""
PART C: Sensitivity analysis
--------------------------------------------------------------------------------
Three plots here
1) Sensitivity of 1ha average to spatial resolution
2) Sensitivity of 1ha average to pulse density
3) Sensitivity of individual profiles to pulse density
--------------------------------------------------------------------------------
"""
# 1 ha average profile sensitivity to spatial resolution
splt.plot_profile_sensitivity_density('%s/sensitivity_analysis/SEOS_sensitivity_analysis_resolution.png' % path2fig,
    heights,resolution_sensitivity)

# 1 ha average profile sensitivity to point density
splt.plot_profile_sensitivity_density('%s/sensitivity_analysis/SEOS_sensitivity_analysis_density.png' % path2fig,
    heights,density_sensitivity)

# 5m x 5m profiles sensitivity to point density
N = 10
for ii in range(0, N):
    profile_idx = np.random.randint(density_sensitivity['20'].shape[1])
    splt.plot_profile_sensitivity_density_individual_profile('%s/sensitivity_analysis/SEOS_sensitivity_analysis_density_individual_5m_profile_%i.png' % (path2fig,profile_idx),
        heights,density_sensitivity,profile_idx=profile_idx)

"""
PART D: 3D map of plant area density
--------------------------------------------------------------------------------
"""
p = pad.plot(x='x',y='y',col='band',col_wrap=8,vmin=0,vmax=1,cmap='plasma',
            cbar_kwargs={'label':'PAI for each 1m slice in canopy'})
for a in p.axes.flat:
    a.set_aspect('equal')
    a.set_xlabel(None)
    a.set_ylabel(None)
    a.set_title(None)
    a.set_xticklabels(a.get_xticklabels(),visible=False)
    a.set_yticklabels(a.get_yticklabels(),visible=False)

p.axes.flat[-1].plot([344672.5,345172.5],[6205442.5,6205442.5],'-',color='black',linewidth='2')
p.axes.flat[-1].annotate('500 m',xy=(344922.5,6205392.5), xycoords='data',
                backgroundcolor='none',ha='center', va='top')
plt.savefig('%s%s_pai_by_canopy_layer.png' % (path2fig,file_id))
plt.show()
