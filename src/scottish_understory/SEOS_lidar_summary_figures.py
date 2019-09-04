"""
SEOS_lidar_summary_figures.py
--------------------------------------------------------------------------------
Code to generate a suite of figures summarising the results from the SEOS
UAV LiDAR surveys.
Analysis focussed on site 1
Encompasses:
- 3D mapping of canopy density
- Canopy density profiles, inverted from point clouds
    - a random selection of 5m x 5m resolution profiles
- Sensitivity analysis to spatial resolution (1 ha)
- Sensitivity analysis to point density (1 ha)
- Sensitivity analysis to point density (5m x 5m individual profiles)
--------------------------------------------------------------------------------
"""
# import standard libraries
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
# Load LiDAR data by bounding box
las_list = '../las_lists/carbomap_site1_lastiles.txt'
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






max_height = 40.
layer_thickness = 0.5
heights = np.arange(0.,max_height,layer_thickness)+0.5
heights_rad = np.arange(0,max_height+0.5,layer_thickness)
n_layers = heights.size
plot_width = 100.
sample_res = 5.
kappa = 1.
n_iter = 100
